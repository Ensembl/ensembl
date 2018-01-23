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

package SeqStoreConverter::GallusGallus;

use vars qw(@ISA);

@ISA = qw(SeqStoreConverter::BasicConverter);


sub create_coord_systems {
  my $self = shift;

  $self->debug("GallusGallus Specific: loading assembly data");

  my $target = $self->target();
  my $dbh    = $self->dbh();

  my $ass_def = $self->get_default_assembly();

  my @coords =
    (["chromosome" ,  $ass_def, "default_version", 1 ],
     ["supercontig",  $ass_def, "default_version", 2 ],
     ["contig"     ,  undef   , "default_version,sequence_level", 3]);

  my @assembly_mappings =  ("chromosome:$ass_def|contig",
                            "supercontig:$ass_def|contig",
                            "chromosome:$ass_def|contig|supercontig:$ass_def");

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

  $self->debug("GallusGallus Specific: creating contig, " .
               "clone, chromosome and supercontig seq_regions");

  $self->contig_to_seq_region();
  $self->chromosome_to_seq_region();
  $self->supercontig_to_seq_region();
}

#
# overridden to do trimming of contig names
#
sub contig_to_seq_region {
  my $self = shift;
  my $target_cs_name = shift;

  my $target = $self->target();
  my $source = $self->source();
  my $dbh     = $self->dbh();

  $target_cs_name ||= 'contig';

  $self->debug("GallusGallus Specific: Transforming contigs into " .
               "$target_cs_name seq_regions");

  my $cs_id = $self->get_coord_system_id($target_cs_name);

  # this ugly SQL simply takes the first part of the contig name
  # but trims everything after and including the second dot
  my $sth = $dbh->prepare
    ("INSERT INTO $target.seq_region " .
     "SELECT contig_id, SUBSTRING(name,1, LOCATE('.',name) + LOCATE('.',SUBSTRING(name,LOCATE('.',name)+1)) -1), $cs_id, length FROM $source.contig");

  $sth->execute();
  $sth->finish();
}

#
# overridden so that left over garbage in chromosome table is not used
#
sub chromosome_to_seq_region {
  my $self = shift;
  my $target_cs_name = shift;

  my $target = $self->target();
  my $source = $self->source();
  my $dbh    = $self->dbh();

  $target_cs_name ||= "chromosome";
  my $cs_id = $self->get_coord_system_id($target_cs_name);

  $self->debug("GallusGallus Specific: Transforming chromosomes into $target_cs_name seq_regions");

  # only take chromosomes which are actually in the assembly table
  my $select_sth = $dbh->prepare
    ("SELECT c.chromosome_id, c.name, c.length " .
     "FROM $source.chromosome c, $source.assembly a " .
     "WHERE c.chromosome_id = a.chromosome_id group by c.chromosome_id");

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

  $self->debug("GallusGallus Specific: loading assembly data");
  $self->assembly_contig_chromosome();
  $self->assembly_contig_supercontig();
}




1;
