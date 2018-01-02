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

use strict;
use warnings;

use SeqStoreConverter::BasicConverter;

package SeqStoreConverter::CaenorhabditisBriggsae;

use vars qw(@ISA);

@ISA = qw(SeqStoreConverter::BasicConverter);


sub create_coord_systems {
  my $self = shift;

  $self->debug("CaenorhabditisBriggsae Specific: creating clone, scaffold," .
              " and contig coordinate systems");

  my $target = $self->target();
  my $dbh    = $self->dbh();

  my $ass_def = $self->get_default_assembly();

  my @coords = 
    (["scaffold" , $ass_def,   "default_version", 1     ],
     ['clone'      , undef   , 'default_version', 2     ],
     ["contig"     , undef   , "default_version,sequence_level", 3]);

  my @assembly_mappings =  ("scaffold:$ass_def|contig",
                            "clone|contig",
                            "scaffold:$ass_def|contig|clone");

  $self->debug("Building coord_system table");

  my $sth = $dbh->prepare("INSERT INTO $target.coord_system " .
                           "(name, version, attrib, rank) VALUES (?,?,?,?)");

  my %coord_system_ids;

  foreach my $cs (@coords) {
    $sth->execute(@$cs);
    $coord_system_ids{$cs->[0]} = $sth->{'mysql_insertid'};
  }
  $sth->finish();

  $self->debug("Adding assembly.mapping entries to meta table");

  $sth = $dbh->prepare("INSERT INTO $target.meta(meta_key, meta_value) " .
                       "VALUES ('assembly.mapping', ?)");

  foreach my $mapping (@assembly_mappings) {
    $sth->execute($mapping);
  }

  $sth->finish();

  return;
}



sub create_seq_regions {
  my $self = shift;

  $self->debug("CaenorhabditisBriggsae Specific: creating contig, " .
               "clone, contig and scaffold seq_regions");

  $self->contig_to_seq_region();
  $self->clone_to_seq_region();
  $self->chromosome_to_seq_region('scaffold');
}


sub chromosome_to_seq_region {
  my $self = shift;
  my $target_cs_name = shift;

  my $target = $self->target();
  my $source = $self->source();
  my $dbh    = $self->dbh();

  $target_cs_name ||= "chromosome";
  my $cs_id = $self->get_coord_system_id($target_cs_name);

  $self->debug("CaenorhabditisBriggsae Specific: Transforming " .
               "chromosomes into $target_cs_name seq_regions");


  ## For consistancy with mart and v19 we need to keep chr name the same for
  ## now, so the following section is commented out and replaced:
  ##strip off the leading 'cb25.' from the chromosome name
  #my $select_sth = $dbh->prepare
  #  ("SELECT chromosome_id,substring(name,6),length FROM $source.chromosome");

  my $select_sth = $dbh->prepare
    ("SELECT chromosome_id,name,length FROM $source.chromosome");


  my $insert_sth = $dbh->prepare
    ("INSERT INTO $target.seq_region (name, coord_system_id, length) " .
     "VALUES (?,?,?)");

  my $tmp_insert_sth = $dbh->prepare
    ("INSERT INTO $target.tmp_chr_map (old_id, new_id) VALUES (?, ?)");

  $select_sth->execute();

  my ($chrom_id, $name, $length);
  $select_sth->bind_columns(\$chrom_id, \$name, \$length);

  while ($select_sth->fetch()) {
    #insert into seq_region table
    $insert_sth->execute($name, $cs_id, $length);
    #copy old/new mapping into temporary table
    $tmp_insert_sth->execute($chrom_id, $insert_sth->{'mysql_insertid'});
  }

  $select_sth->finish();
  $insert_sth->finish();
  $tmp_insert_sth->finish();

  return;
}


sub create_assembly {
  my $self = shift;

  $self->debug("CaenorhabditisBriggsae Specific: loading assembly data");

  $self->assembly_contig_chromosome();
  $self->assembly_contig_clone();
}




#
# Override the assembly contig clone method because the briggsae database
# does not have any embl_offsets
#
sub assembly_contig_clone {
  my $self = shift;

  my $target = $self->target();
  my $source = $self->source();
  my $dbh    = $self->dbh();


  $self->debug("CaenorhabditisBriggsae Specific: loading contig/clone " .
               "assembly relationship");

  my $asm_sth = $dbh->prepare
    ("INSERT INTO $target.assembly " .
     "set asm_seq_region_id = ?, ".
     "    asm_start = ?, " .
     "    asm_end   = ?, " .
     "    cmp_seq_region_id = ?, ".
     "    cmp_start = ?, " .
     "    cmp_end   = ?, " .
     "    ori       = ?");

  # get a list of the contigs that have clones, their ids, and the
  # corresponding clone ids
  my $ctg_sth = $dbh->prepare
    ("SELECT ctg.name, ctg.contig_id, ctg.length, cln.new_id " .
     "FROM   $source.contig ctg, $target.tmp_cln_map cln " .
     "WHERE  ctg.name not like 'c%' " .  # only contigs w/ proper accessions
     "AND    ctg.clone_id = cln.old_id");

  $ctg_sth->execute();

  my ($ctg_name, $ctg_id, $ctg_len, $cln_id);

  $ctg_sth->bind_columns(\$ctg_name, \$ctg_id, \$ctg_len, \$cln_id);

  while($ctg_sth->fetch()) {
    my (undef,$cln_start, $cln_end) = split(/\./, $ctg_name);
    my $cln_len = $cln_end - $cln_start + 1;
    if($cln_len != $ctg_len) {
      die("Contig len $ctg_len != Clone len $cln_len");
    }

    $asm_sth->execute($cln_id, $cln_start, $cln_end,
                      $ctg_id, 1, $ctg_len, 1);
  }

  $ctg_sth->finish();
  $asm_sth->finish();

  return;
}



#
# Override contig_to_seq_region and clone_to_seq_region to provide 
# briggsae specific behaviour
#

# sub contig_to_seq_region {
#   my $self = shift;
#   my $target_cs_name = shift;

#   my $target = $self->target();
#   my $source = $self->source();
#   my $dbh     = $self->dbh();

#   $target_cs_name ||= 'contig';

#   $self->debug("CaenorhabditisBriggsae Specific: Transforming contigs into " .
#                "$target_cs_name seq_regions");

#   my $cs_id = $self->get_coord_system_id($target_cs_name);

#   #There are two types of contigs in briggsae:

#   #
#   # cosmids/clones
#   #
#   my $sth = $dbh->prepare
#     ("INSERT INTO $target.seq_region " .
#      "SELECT contig_id, name, $cs_id, length " .
#      "FROM $source.contig " .
#      "WHERE  name not like 'c%'");

#   $sth->execute();
#   $sth->finish();

#   #
#   # WGS contigs
#   #
#   $sth = $dbh->prepare
#     ("INSERT INTO $target.seq_region " .
#      "SELECT ctg.contig_id, cln.name, $cs_id, length " .
#      "FROM   $source.contig ctg, $source.clone cln " .
#      "WHERE  ctg.clone_id = cln.clone_id " .
#      "AND    ctg.name like 'c%'");

#   $sth->execute();
#   $sth->finish();

#   return;
# }



sub clone_to_seq_region {
  my $self = shift;
  my $target_cs_name = shift;

  my $target = $self->target();
  my $source = $self->source();
  my $dbh    = $self->dbh();

  # target coord_system will have a different ID
  $target_cs_name ||= "clone";
  my $cs_id = $self->get_coord_system_id($target_cs_name);

  $self->debug("CaenorhabditisBriggsae Specific:Transforming clones " .
               "into $target_cs_name seq_regions");

  #
  # We don't want to make clones out of the WGS contigs, only out of
  # the clones with proper embl accessions.  Also for some reason the embl_offset
  # is not set in the briggsae 17/18/19 databases, which means we have to deduce the
  # length from the name of the contigs!
  #
  my $select_sth = $dbh->prepare
    ("SELECT cl.clone_id,
             CONCAT(cl.embl_acc, '.', cl.embl_version),
             ctg.name
     FROM   $source.clone cl, $source.contig ctg
		 WHERE  cl.clone_id = ctg.clone_id
     AND    cl.embl_acc not like 'c%'
     ORDER BY cl.clone_id");

  $select_sth->execute();

  my ($clone_id, $embl_acc, $ctg_name);
  $select_sth->bind_columns(\$clone_id, \$embl_acc, \$ctg_name);

  my $highest_end = undef;
  my $current_clone = undef;
  my $current_clone_id = undef;
  my $length;

  my $insert_sth = $dbh->prepare
    ("INSERT INTO $target.seq_region (name, coord_system_id, length) " .
     "VALUES(?,?,?)");

  my $tmp_insert_sth = $dbh->prepare
    ("INSERT INTO $target.tmp_cln_map (old_id, new_id) VALUES (?, ?)");

  while ($select_sth->fetch()) {
    #extract the end position of the contig
    my $ctg_end;
    (undef,undef,$ctg_end) = split(/\./, $ctg_name);

    if(!defined($current_clone)) {
      $current_clone = $embl_acc;
      $current_clone_id = $clone_id;
      $highest_end   = $ctg_end;
    }

    if($current_clone ne $embl_acc) {
      #started new clone, store last one

      $insert_sth->execute($current_clone, $cs_id, $highest_end);
      #store mapping of old -> new ids in temp table
      $tmp_insert_sth->execute($current_clone_id, $insert_sth->{'mysql_insertid'});

      $current_clone = $embl_acc;
      $current_clone_id = $clone_id;
      $highest_end = $ctg_end;
    } elsif($ctg_end > $highest_end) {
      #same clone, adjust end if end of contig is highest yet seen
      $highest_end = $ctg_end;
    }
  }

  #insert the last clone
  $insert_sth->execute($current_clone, $cs_id, $highest_end);
  $tmp_insert_sth->execute($current_clone_id, $insert_sth->{'mysql_insertid'});


  $select_sth->finish();
  $insert_sth->finish();
  $tmp_insert_sth->finish();

  return;
}


1;
