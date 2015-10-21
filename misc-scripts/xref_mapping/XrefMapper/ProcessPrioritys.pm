=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
#        Also do this fro its depenedents
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
  my ($self) = @_;


  my $psth = $self->xref->dbc->prepare("select s.priority_description, s.name from source s, xref x where x.source_id = s.source_id group by s.priority_description, s.name order by s.name") || die "prepare failed";
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

  my @names = $self->get_priority_names();

  print "The foillowing will be processed as priority xrefs\n" if($self->verbose);
  foreach my $name (@names){
    print "\t$name\n" if($self->verbose);
  }

  my $update_ox_sth = $self->xref->dbc->prepare('update object_xref set ox_status = "FAILED_PRIORITY" where object_xref_id = ?');
  my $update_x_sth  = $self->xref->dbc->prepare("update xref set dumped = 'NO_DUMP_ANOTHER_PRIORITY' where xref_id = ?");

  #
  # Change of tact here to make the sql easier...
  #

  # 1) Set to failed all those that have no object xrefs.

  my $f_sql =(<<FSQL);
   SELECT x.xref_id
     FROM  source s, xref x 
     LEFT JOIN object_xref ox ON  ox.xref_id = x.xref_id 
       WHERE x.source_id = s.source_id
         AND s.name = ? 
         AND ox.object_xref_id is null
FSQL

  my $f_sth =  $self->xref->dbc->prepare($f_sql);
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
          AND ox_status not in ('FAILED_CUTOFF')
          AND s.name = ?
         ORDER BY x.accession DESC, s.priority ASC , identity DESC, x.xref_id DESC
NEWS

  #
  # Query to copy identity_xref values from one xref to another
  # This is to keep alignment information event if alignment was not the best match
  #

  my $idx_copy_sql = (<<IDXCP);
   UPDATE identity_xref SET query_identity = ?, target_identity = ?, hit_start = ?, hit_end = ?, translation_start = ?, translation_end = ?, cigar_line = ?, score = ?, evalue = ?
      WHERE object_xref_id = ?;
IDXCP

  my $idx_copy_sth = $self->xref->dbc->prepare($idx_copy_sql);

  #
  # Query to copy synonyms from one xref to another
  #

  my $syn_copy_sql = (<<SYNCP);
    INSERT IGNORE INTO synonym (SELECT ?, synonym FROM synonym
       WHERE xref_id = ?);
SYNCP

  my $syn_copy_sth = $self->xref->dbc->prepare($syn_copy_sql);

  my $best_ox_sth = $self->xref->dbc->prepare("SELECT object_xref_id FROM object_xref WHERE xref_id = ? and ensembl_object_type = ? and ensembl_id = ?");

  my $seq_score_sql = (<<SEQCP);
    SELECT query_identity, target_identity, hit_start, hit_end, translation_start, translation_end, cigar_line, score, evalue
      FROM identity_xref WHERE object_xref_id = ?
SEQCP
  my $seq_score_sth = $self->xref->dbc->prepare($seq_score_sql);

  my $sth = $self->xref->dbc->prepare($new_sql);
  foreach my $name (@names){
    $sth->execute($name);
    my ($object_xref_id, $acc, $xref_id, $identity, $status, $object_type, $ensembl_id, $info_type);
    $sth->bind_columns(\$object_xref_id, \$acc, \$xref_id, \$identity, \$status, \$object_type, \$ensembl_id, \$info_type);
    my $last_acc = "";
    my $last_name = "";
    my $best_xref_id = undef;
    my @best_ensembl_id = undef;
    my $last_xref_id = 0;
    my $seen = 0;
    my @gone;
    while($sth->fetch){
      if($last_acc eq $acc){
	if($xref_id != $best_xref_id){
          if ($xref_id == $last_xref_id) {
            $seen = 1;
          } else {
            $seen = 0;
          }
          $last_xref_id = $xref_id;
# If it is a sequence_match, we want to copy the alignment identity_xref to prioritised mappings to the same ensembl_id
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
	  if($status eq "DUMP_OUT"){
            $update_ox_sth->execute($object_xref_id);
## If it is the first time processing this xref_id, also process dependents and update status
            if (!$seen) {
              $update_x_sth->execute($xref_id);
# Copy synonyms across if they are missing
              $syn_copy_sth->execute($best_xref_id, $xref_id);
              $self->process_dependents($xref_id, $best_xref_id, $object_type, \@best_ensembl_id);
            }
          }
	  else{
	    $update_x_sth->execute($xref_id);
	  }
	} else {
          ## There might be several mappings for the best priority
          push @best_ensembl_id, $ensembl_id;
        }
	if(@gone){ #best priority failed so anothre one now found so set dumped;
	  if($last_name eq $acc){
	    foreach my $d (@gone){
	      $update_x_sth->execute($d);
	    }
	  }
	}
      }
      else{ # NEW
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
	else{
	  if($last_name eq $acc){
	    push @gone, $xref_id;	  }
	  else{
	    @gone = ();
	    push @gone, $xref_id;
	  }	
	  $last_name = $acc;
	}
      }
    }
  }
  $sth->finish;

  $update_ox_sth->finish;
  $update_x_sth->finish;
  $seq_score_sth->finish;
  $best_ox_sth->finish;
  $idx_copy_sth->finish;
  $syn_copy_sth->finish;


# We want to make sure that if a priority xref is NOT MAPPEd then we only



  $sth = $self->xref->dbc->prepare("insert into process_status (status, date) values('prioritys_flagged',now())");
  $sth->execute();
  $sth->finish;
}

sub process_dependents{
  my ($self, $old_master_xref_id, $new_master_xref_id, $object_type, $best_ensembl_id) = @_;

  my $dep_sth           = $self->xref->dbc->prepare("select distinct dependent_xref_id, dx.linkage_annotation, dx.linkage_source_id from dependent_xref dx where dx.master_xref_id = ?");
  my $insert_dep_x_sth  = $self->xref->dbc->prepare("insert into dependent_xref(master_xref_id, dependent_xref_id, linkage_annotation, linkage_source_id) values(?, ?, ?, ?)");
  my $remove_dep_ox_sth = $self->xref->dbc->prepare("delete ix, g from object_xref ox left join identity_xref ix on ix.object_xref_id = ox.object_xref_id left join go_xref g on g.object_xref_id = ox.object_xref_id where master_xref_id = ? and ensembl_object_type = ? and xref_id = ? and ensembl_id = ?");
  my $old_ens_id_sth    = $self->xref->dbc->prepare("select distinct ensembl_id from object_xref where xref_id = ?");
  my $update_dep_ox_sth = $self->xref->dbc->prepare("update ignore object_xref set ox_status = 'FAILED_PRIORITY' where master_xref_id = ? and ensembl_object_type = ? and xref_id = ? and ensembl_id = ? and ox_status = 'DUMP_OUT'");
  my $clean_dep_ox_sth  = $self->xref->dbc->prepare("delete from object_xref where master_xref_id = ? and ensembl_object_type = ? and xref_id = ? and ensembl_id = ? and ox_status = 'DUMP_OUT'");
  my $insert_dep_ox_sth = $self->xref->dbc->prepare("insert ignore into object_xref(master_xref_id, ensembl_object_type, ensembl_id, linkage_type, ox_status, xref_id) values(?, ?, ?, 'DEPENDENT', 'DUMP_OUT', ?)");
  my $dep_ox_sth        = $self->xref->dbc->prepare("select object_xref_id from object_xref where master_xref_id = ? and ensembl_object_type = ? and ensembl_id = ? and linkage_type = 'DEPENDENT' AND ox_status = 'DUMP_OUT' and xref_id = ?");
  my $insert_dep_go_sth = $self->xref->dbc->prepare("insert ignore into go_xref values(?, ?, ?)");
  my $insert_ix_sth     = $self->xref->dbc->prepare("insert ignore into identity_xref(object_xref_id, query_identity, target_identity) values(?, 100, 100)");

  my @master_xrefs;
  my @old_ensembl_ids;
  my $old_ensembl_id;
  my $recursive = 0;

  $old_ens_id_sth->execute($old_master_xref_id);
  $old_ens_id_sth->bind_columns(\$old_ensembl_id);
  my $skip = 0;
  while ($old_ens_id_sth->fetch()) {
    $skip = 0;
    foreach my $best_ensembl_id (@$best_ensembl_id) {
      if ($old_ensembl_id eq $best_ensembl_id) {
        $skip = 1;
        last;
      }
    }
    if ($skip) {
      next;
    }
    push @old_ensembl_ids, $old_ensembl_id;
  }
  $old_ens_id_sth->finish();

  push @master_xrefs, $old_master_xref_id;

  while(my $xref_id = pop(@master_xrefs)){
    my ($dep_xref_id, $linkage_type, $new_object_xref_id, $linkage_source_id);
    $dep_sth->execute($xref_id);
    $dep_sth->bind_columns(\$dep_xref_id, \$linkage_type, \$linkage_source_id);
    ## Loop through all dependent xrefs
    while($dep_sth->fetch()){
      ## Add dependent to priority xref if it is the first master
      if (!$recursive) {
        $insert_dep_x_sth->execute($new_master_xref_id, $dep_xref_id, $linkage_type, $linkage_source_id);
        foreach my $best_ensembl_id (@$best_ensembl_id) {
          $remove_dep_ox_sth->execute($xref_id, $object_type, $dep_xref_id, $best_ensembl_id);
          $update_dep_ox_sth->execute($xref_id, $object_type, $dep_xref_id, $best_ensembl_id);
          $clean_dep_ox_sth->execute($xref_id, $object_type, $dep_xref_id, $best_ensembl_id);
        }
      }
      ## object_xref for dismissed xref set to FAILED_PRIORITY
      ## also delete any leftover identity or go xrefs
      foreach my $old_id (@old_ensembl_ids) {
        $remove_dep_ox_sth->execute($xref_id, $object_type, $dep_xref_id, $old_id);
        $update_dep_ox_sth->execute($xref_id, $object_type, $dep_xref_id, $old_id);
        $clean_dep_ox_sth->execute($xref_id, $object_type, $dep_xref_id, $old_id);
      }
      ## Loop through all ensembl ids mapped to priority xref
      foreach my $best_ensembl_id (@$best_ensembl_id) {
        ## new object_xref for each transfered mapping
        $insert_dep_ox_sth->execute($new_master_xref_id, $object_type, $best_ensembl_id, $dep_xref_id);
        if ($linkage_type) {
          ## If there is a linkage_type, it is a go xref
          ## Fetch the newly created object_xref to add them to go_xref
          $dep_ox_sth->execute($new_master_xref_id, $object_type, $best_ensembl_id, $dep_xref_id);
          $dep_ox_sth->bind_columns(\$new_object_xref_id);
          while ($dep_ox_sth->fetch()) {
            $insert_dep_go_sth->execute($new_object_xref_id, $linkage_type, $new_master_xref_id); 
            $insert_ix_sth->execute($new_object_xref_id);
          }
        }
      }
      push @master_xrefs, $dep_xref_id; # remember dependents dependents
    }
    $recursive = 1;
    $new_master_xref_id = $dep_xref_id;
  }

  $dep_sth->finish();
  $insert_dep_x_sth->finish();
  $remove_dep_ox_sth->finish();
  $clean_dep_ox_sth->finish();
  $insert_dep_ox_sth->finish();
  $dep_ox_sth->finish();
  $insert_dep_go_sth->finish();

}

1;
