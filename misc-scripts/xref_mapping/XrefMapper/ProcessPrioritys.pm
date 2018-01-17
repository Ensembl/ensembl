=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package XrefMapper::ProcessPrioritys;
use strict;

use vars '@ISA';
@ISA = qw{ XrefMapper::BasicMapper };

use warnings;
use XrefMapper::BasicMapper;

use Cwd;
use DBI;
use File::Basename;
use IPC::Open3;

# Process the priority xrefs.

#
# 1) create a list of source "names" that are priority xrefs
#
# 2) Just to be sure set all ox_status in object_xref to 'DUMP_OUT'
#    set dumped in xref to NULL
# 
# 3) for each of the source names 
#    set ox_status to 'FAILED_PRIORITY' for those not the best match
#        Also do this for its depenedents
#

sub new {
  my($class, $mapper) = @_;

  my $self ={};
  bless $self,$class;
#  $self->core($mapper->core);
  $self->xref($mapper->xref);
  $self->verbose($mapper->verbose);
  return $self;
}

sub get_priority_names{
  my ($self, $dbi) = @_;


  my $psth = $dbi->prepare("select s.priority_description, s.name from source s, xref x where x.source_id = s.source_id group by s.priority_description, s.name order by s.name") || die "prepare failed";
  $psth->execute() || die "execute failed";

  my @names;
  my %seen;

  my $last_name = "rubbish";
  my ($desc,$name);
  $psth->bind_columns(\$desc,\$name);
  while($psth->fetch()){
    if($name eq $last_name and !defined($seen{$name})){
      push @names, $name;
      $seen{$name} = 1;
    }
    $last_name = $name;
  }

  return @names;
}


sub process {
  my ($self) = @_;

  my $dbi = $self->xref->dbc;
  my @names = $self->get_priority_names($dbi);

  print "The following will be processed as priority xrefs\n" if($self->verbose);
  foreach my $name (@names){
    print "\t$name\n" if($self->verbose);
  }

  my $update_ox_sth = $dbi->prepare('update object_xref set ox_status = "FAILED_PRIORITY" where object_xref_id = ?');
  my $update_x_sth  = $dbi->prepare("update xref set dumped = 'NO_DUMP_ANOTHER_PRIORITY' where xref_id = ?");

  # 1) Set to failed all those that have no object xrefs.

  my $f_sql =(<<FSQL);
   SELECT x.xref_id
     FROM  source s, xref x 
     LEFT JOIN object_xref ox ON  ox.xref_id = x.xref_id 
       WHERE x.source_id = s.source_id
         AND s.name = ? 
         AND ox.object_xref_id is null
FSQL

  my $f_sth =  $dbi->prepare($f_sql);
  foreach my $name (@names){
    $f_sth->execute($name);
    my ($xref_id);
    $f_sth->bind_columns(\$xref_id);
    while($f_sth->fetch()){
      $update_x_sth->execute($xref_id);
    }
  }
  $f_sth->finish;


  # 
  # Now ALL object_xrefs have an identity_xref :-)
  # So we can do a straight join and treat all info_types the same way.
  # 
  my $new_sql =(<<NEWS);
   SELECT ox.object_xref_id, x.accession, x.xref_id, (ix.query_identity + ix.target_identity) as identity, ox.ox_status, ox.ensembl_object_type, ensembl_id, info_type
      FROM object_xref ox, xref x, source s, identity_xref ix
        WHERE ox.object_xref_id = ix.object_xref_id 
          AND ox.xref_id = x.xref_id
          AND s.source_id = x.source_id
          AND s.name = ?
         ORDER BY x.accession DESC, s.priority ASC , identity DESC, x.xref_id DESC
NEWS
  my $new_sth = $dbi->prepare($new_sql);
  #
  # Query to copy identity_xref values from one xref to another
  # This is to keep alignment information event if alignment was not the best match
  #

  my $idx_copy_sql = (<<IDXCP);
   UPDATE identity_xref SET query_identity = ?, target_identity = ?, hit_start = ?, hit_end = ?, translation_start = ?, translation_end = ?, cigar_line = ?, score = ?, evalue = ?
      WHERE object_xref_id = ?;
IDXCP

  my $idx_copy_sth = $dbi->prepare($idx_copy_sql);

  #
  # Query to copy synonyms from one xref to another
  #

  my $syn_copy_sql = (<<SYNCP);
    INSERT IGNORE INTO synonym (SELECT ?, synonym FROM synonym
       WHERE xref_id = ?);
SYNCP

  my $syn_copy_sth = $dbi->prepare($syn_copy_sql);

  my $best_ox_sth = $dbi->prepare("SELECT object_xref_id FROM object_xref WHERE xref_id = ? and ensembl_object_type = ? and ensembl_id = ?");

  my $seq_score_sql = (<<SEQCP);
    SELECT query_identity, target_identity, hit_start, hit_end, translation_start, translation_end, cigar_line, score, evalue
      FROM identity_xref WHERE object_xref_id = ?
SEQCP
  my $seq_score_sth = $dbi->prepare($seq_score_sql);


  foreach my $name (@names){
    $new_sth->execute($name);
    my ($object_xref_id, $acc, $xref_id, $identity, $status, $object_type, $ensembl_id, $info_type);
    $new_sth->bind_columns(\$object_xref_id, \$acc, \$xref_id, \$identity, \$status, \$object_type, \$ensembl_id, \$info_type);
    my $last_acc = "";
    my $last_name = "";
    my $best_xref_id = undef;
    my @best_ensembl_id = undef;
    my $last_xref_id = 0;
    my $seen = 0;
    my @gone; # list of xref_ids that we've already seen for this accession
    while($new_sth->fetch){
      if($last_acc eq $acc){
        if($xref_id != $best_xref_id){
          # We've already seen this accession before, and this xref_id is not the best one

          $seen = ($xref_id == $last_xref_id);
          
          $last_xref_id = $xref_id;
# If xref is a sequence_match, we want to copy the alignment identity_xref to prioritised mappings of the same ensembl_id
          if ($info_type eq 'SEQUENCE_MATCH') {
            my ($query_identity, $target_identity, $hit_start, $hit_end, $translation_start, $translation_end, $cigar_line, $score, $evalue, $best_object_xref_id);
            $seq_score_sth->execute($object_xref_id);
            $seq_score_sth->bind_columns(\$query_identity, \$target_identity, \$hit_start, \$hit_end, \$translation_start, \$translation_end, \$cigar_line, \$score, \$evalue);
            $seq_score_sth->fetch();
            $best_ox_sth->execute($best_xref_id, $object_type, $ensembl_id);
            $best_ox_sth->bind_columns(\$best_object_xref_id);
            $best_ox_sth->fetch();
            $idx_copy_sth->execute($query_identity, $target_identity, $hit_start, $hit_end, $translation_start, $translation_end, $cigar_line, $score, $evalue, $best_object_xref_id);
          }
          # If the xref is marked DUMP_OUT, set it to FAILED_PRIORITY
          if($status eq "DUMP_OUT"){
            $update_ox_sth->execute($object_xref_id);
## If it is the first time processing this xref_id, also process dependents and update status
            if (!$seen) {
              $update_x_sth->execute($xref_id);
# Copy synonyms across if they are missing
              $syn_copy_sth->execute($best_xref_id, $xref_id);
              $self->process_dependents($xref_id, $best_xref_id, $dbi);
            }
          }
          else{ # not DUMP_OUT
            $update_x_sth->execute($xref_id);
          }
        } else {
# Alignment did not pass, dismiss
          if ($status eq 'FAILED_CUTOFF') {
            next;
          }
          ## There might be several mappings for the best priority
          push @best_ensembl_id, $ensembl_id;
        }
        if(@gone){ #best priority failed so another one now found so set dumped;
          if($last_name eq $acc){
            foreach my $d (@gone){
              $update_x_sth->execute($d);
            }
          }
        }
      }
      else{ # NEW xref_id
        if($status eq "DUMP_OUT"){
          $last_acc = $acc;
          $best_xref_id = $xref_id;
          @best_ensembl_id = ($ensembl_id);
          if(@gone and $last_name eq $acc){
            foreach my $d (@gone){
              $update_x_sth->execute($d);
            }
            @gone=();
          }
        }
        else{ # new xref_id not DUMP_OUT
          if ($last_name ne $acc) { @gone = () } # new accession
          push @gone, $xref_id;	
          $last_name = $acc;
        }
      }
    }
  }
  $new_sth->finish;

  $update_ox_sth->finish;
  $update_x_sth->finish;
  $seq_score_sth->finish;
  $best_ox_sth->finish;
  $idx_copy_sth->finish;
  $syn_copy_sth->finish;

  my $sth = $dbi->prepare("insert into process_status (status, date) values('prioritys_flagged',now())");
  $sth->execute();
  $sth->finish;
}

sub process_dependents{
# master xref IDs are entries for the current accession via various methods. We take dependent xrefs from the old and add to the new
  my ($self, $old_master_xref_id, $new_master_xref_id, $dbi) = @_;


  my $matching_ens_sth  = $dbi->prepare("select distinct ensembl_object_type, ensembl_id from object_xref where ox_status not in ('FAILED_CUTOFF') and xref_id = ? order by ensembl_object_type");
  my $dep_sth           = $dbi->prepare("select distinct dx.dependent_xref_id, dx.linkage_annotation, dx.linkage_source_id, ox.ensembl_object_type from dependent_xref dx, object_xref ox where ox.xref_id = dx.dependent_xref_id and ox.master_xref_id = dx.master_xref_id and dx.master_xref_id = ? order by ox.ensembl_object_type");
  my $insert_dep_x_sth  = $dbi->prepare("insert into dependent_xref(master_xref_id, dependent_xref_id, linkage_annotation, linkage_source_id) values(?, ?, ?, ?)");
  my $insert_dep_ox_sth = $dbi->prepare("insert ignore into object_xref(master_xref_id, ensembl_object_type, ensembl_id, linkage_type, ox_status, xref_id) values(?, ?, ?, 'DEPENDENT', 'DUMP_OUT', ?)");
  my $dep_ox_sth        = $dbi->prepare("select object_xref_id from object_xref where master_xref_id = ? and ensembl_object_type = ? and ensembl_id = ? and linkage_type = 'DEPENDENT' AND ox_status = 'DUMP_OUT' and xref_id = ?");
  my $insert_dep_go_sth = $dbi->prepare("insert ignore into go_xref values(?, ?, ?)");
  my $insert_ix_sth     = $dbi->prepare("insert ignore into identity_xref(object_xref_id, query_identity, target_identity) values(?, 100, 100)");

  my @master_xrefs = ($old_master_xref_id);
  my $recursive = 0;

  my ($new_object_type, $new_ensembl_id, $old_object_type, $old_ensembl_id);
  my ($dep_xref_id, $linkage_annotation, $new_object_xref_id, $linkage_source_id, $object_type);


  # Create a hash of all possible mappings for this accession
  my %ensembl_ids;
  $matching_ens_sth->execute($new_master_xref_id);
  $matching_ens_sth->bind_columns(\$new_object_type, \$new_ensembl_id);
  while ($matching_ens_sth->fetch()) {
    push @{ $ensembl_ids{$new_object_type} }, $new_ensembl_id;
  }
  my %old_ensembl_ids;
  $matching_ens_sth->execute($old_master_xref_id);
  $matching_ens_sth->bind_columns(\$old_object_type, \$old_ensembl_id);
  while ($matching_ens_sth->fetch()) {
    push @{ $old_ensembl_ids{$old_object_type} }, $old_ensembl_id;
  }


  ## Loop through all dependent xrefs of old master xref, and recurse
  while(my $xref_id = pop(@master_xrefs)){ 
    
    # Get dependent xrefs, be they gene, transcript or translation
    $dep_sth->execute($xref_id);
    $dep_sth->bind_columns(\$dep_xref_id, \$linkage_annotation, \$linkage_source_id, \$object_type);
    if ($recursive) {
      $new_master_xref_id = $xref_id;
    }
    while($dep_sth->fetch()){


      # Remove all mappings to low priority xrefs
      # Then delete any leftover identity or go xrefs of it
      foreach my $ensembl_id (@{ $old_ensembl_ids{$object_type}} ) {
        $self->_detach_object_xref($xref_id, $dep_xref_id, $object_type, $ensembl_id, $dbi);
      }

      # Duplicate each dependent for the new master xref if it is the first in the chain
      unless ($recursive) {
        $insert_dep_x_sth->execute($new_master_xref_id, $dep_xref_id, $linkage_annotation, $linkage_source_id);
      }

      # Loop through all chosen (best) ensembl ids mapped to priority xref, and connect them with object_xrefs
      foreach my $ensembl_id (@{ $ensembl_ids{$object_type} }) {
        # Add new object_xref for each best_ensembl_id. 
        $insert_dep_ox_sth->execute($new_master_xref_id, $object_type, $ensembl_id, $dep_xref_id);
        ## If there is a linkage_annotation, it is a go xref
        if ($linkage_annotation) {
          ## Fetch the newly created object_xref to add them to go_xref
          $dep_ox_sth->execute($new_master_xref_id, $object_type, $ensembl_id, $dep_xref_id);
          $dep_ox_sth->bind_columns(\$new_object_xref_id);
          while ($dep_ox_sth->fetch()) {
            $insert_dep_go_sth->execute($new_object_xref_id, $linkage_annotation, $new_master_xref_id);
            $insert_ix_sth->execute($new_object_xref_id);
          }
        }
      }
      unless ($dep_xref_id == $xref_id) {
        push @master_xrefs, $dep_xref_id; # remember chained dependent xrefs
      }
    }
    $recursive = 1;
  }

  $matching_ens_sth->finish();
  $dep_sth->finish();
  $insert_dep_x_sth->finish();
  $insert_dep_ox_sth->finish();
  $dep_ox_sth->finish();
  $insert_dep_go_sth->finish();
  $insert_ix_sth->finish();
}

# Delete identity xrefs, go_xrefs for a given object xref
# Set unimportant object_xrefs to FAILED_PRIORITY, and delete all those that remain
sub _detach_object_xref {
  my $self = shift;
  my ($xref_id, $dep_xref_id, $object_type, $ensembl_id, $dbi) = @_;
  # Drop all the identity and go xrefs for the dependents of an xref
  my $remove_dep_ox_sth = $dbi->prepare(
    "DELETE ix, g FROM object_xref ox \
     LEFT JOIN identity_xref ix ON ix.object_xref_id = ox.object_xref_id \
     LEFT JOIN go_xref g ON g.object_xref_id = ox.object_xref_id \
     WHERE master_xref_id = ? AND ensembl_object_type = ? AND xref_id = ? AND ensembl_id = ?"
  );
  # Fail the object_xrefs that did link to the deleted identity/go xrefs.
  # This only updates one of potentially many, due to table contraints.
  my $update_dep_ox_sth = $dbi->prepare(
    "UPDATE IGNORE object_xref SET ox_status = 'FAILED_PRIORITY' \
    WHERE master_xref_id = ? AND ensembl_object_type = ? AND xref_id = ? AND ox_status = 'DUMP_OUT' AND ensembl_id = ?"
  );
  # This deletes everything left behind by the previous query.
  my $clean_dep_ox_sth  = $dbi->prepare(
    "DELETE FROM object_xref \
    WHERE master_xref_id = ? AND ensembl_object_type = ? AND xref_id = ? AND ox_status = 'DUMP_OUT' AND ensembl_id = ?"
  );

  $remove_dep_ox_sth->execute($xref_id, $object_type, $dep_xref_id, $ensembl_id);
  # change status of object_xref to FAILED_PRIORITY for record keeping
  $update_dep_ox_sth->execute($xref_id, $object_type, $dep_xref_id, $ensembl_id);
  # delete the duplicates.
  $clean_dep_ox_sth->execute($xref_id, $object_type, $dep_xref_id, $ensembl_id);

  $remove_dep_ox_sth->finish();
  $update_dep_ox_sth->finish();
  $clean_dep_ox_sth->finish();

}


1;
