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

package SeqStoreConverter::MusMusculus;

use vars qw(@ISA);

@ISA = qw(SeqStoreConverter::BasicConverter);


sub create_coord_systems {
  my $self = shift;

  $self->debug("MusMusculus Specific: loading assembly data");

  my $target = $self->target();
  my $dbh    = $self->dbh();

  my $ass_def = $self->get_default_assembly();

  my @coords = 
    (["chromosome" , $ass_def, "default_version"          ,1    ],
     ["supercontig", undef   , "default_version"          ,2     ],
     ['clone'      , undef   , 'default_version'          ,3     ],
     ["contig"     , undef   , "default_version,sequence_level",4]);

  my @assembly_mappings =  ("chromosome:$ass_def|contig",
                            "supercontig|contig",
                            "clone|contig",
                            "chromosome:$ass_def|contig|clone",
                            "chromosome:$ass_def|contig|supercontig",
                            "supercontig|contig|clone");

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

  $self->debug("MusMusculus Specific: creating contig, " .
               "clone, chromosome and supercontig seq_regions");

  $self->contig_to_seq_region();
  $self->clone_to_seq_region();
  $self->chromosome_to_seq_region();
  $self->supercontig_to_seq_region();
}

sub create_assembly {
  my $self = shift;

  $self->debug("MusMusculus Specific: loading assembly data");

  $self->assembly_contig_chromosome();
  $self->assembly_contig_clone();
  $self->assembly_contig_supercontig();
}

#
# Override contig_to_seq_region and clone_to_seq_region to provide 
# mouse specific behaviour
#

# sub contig_to_seq_region {
#   my $self = shift;
#   my $target_cs_name = shift;

#   my $target = $self->target();
#   my $source = $self->source();
#   my $dbh     = $self->dbh();

#   $target_cs_name ||= 'contig';

#   $self->debug("MusMusculus Specific: Transforming contigs into " .
#                "$target_cs_name seq_regions");

#   my $cs_id = $self->get_coord_system_id($target_cs_name);

#   #There are two types of contigs in mouse:

#   #
#   # Contigs which form BAC clones
#   #
#   my $sth = $dbh->prepare
#     ("INSERT INTO $target.seq_region " .
#      "SELECT contig_id, name, $cs_id, length " .
#      "FROM $source.contig " .
#      "WHERE  name not like 'C%'");

#   $sth->execute();
#   $sth->finish();

#   #
#   # Contigs which were created from whole genome shotgun
#   #
#   $sth = $dbh->prepare
#     ("INSERT INTO $target.seq_region " .
#      "SELECT ctg.contig_id, cln.name, $cs_id, length " .
#      "FROM   $source.contig ctg, $source.clone cln " .
#      "WHERE  ctg.clone_id = cln.clone_id " .
#      "AND    ctg.name like 'C%'");

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

  $self->debug("MusMusculus Specific:Transforming clones " .
               "into $target_cs_name seq_regions");

  #
  # We don't want to make clones out of the WGS contigs, only out of
  # the actual BACs with proper embl accessions
  #
  my $select_sth = $dbh->prepare
    ("SELECT cl.clone_id,
             CONCAT(cl.embl_acc, '.', cl.embl_version),
             MAX(ctg.embl_offset+ctg.length-1)
     FROM   $source.clone cl, $source.contig ctg
		 WHERE  cl.clone_id = ctg.clone_id
     AND    cl.embl_acc not like 'C%'
     GROUP BY ctg.clone_id");
 
  $select_sth->execute();

  my ($clone_id, $embl_acc, $length);
  $select_sth->bind_columns(\$clone_id, \$embl_acc, \$length);

  my $insert_sth = $dbh->prepare
    ("INSERT INTO $target.seq_region (name, coord_system_id, length) " .
     "VALUES(?,?,?)");

  my $tmp_insert_sth = $dbh->prepare
    ("INSERT INTO $target.tmp_cln_map (old_id, new_id) VALUES (?, ?)");

  while ($select_sth->fetch()) {
    $insert_sth->execute("$embl_acc", $cs_id, $length);

    #store mapping of old -> new ids in temp table
    $tmp_insert_sth->execute($clone_id, $insert_sth->{'mysql_insertid'});
  }

  $select_sth->finish();
  $insert_sth->finish();
  $tmp_insert_sth->finish();

  return;
}


1;
