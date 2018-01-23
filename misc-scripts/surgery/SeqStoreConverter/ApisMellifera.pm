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

package SeqStoreConverter::ApisMellifera;

use vars qw(@ISA);

@ISA = qw(SeqStoreConverter::BasicConverter);


sub create_coord_systems {
  my $self = shift;

  $self->debug("ApisMellifera Specific: creating " .
               "contig and scaffold coord systems");

  my $target = $self->target();
  my $dbh    = $self->dbh();

  my $ass_def = $self->get_default_assembly();

  my @coords = 
    (["scaffold" , $ass_def, "default_version", 1     ],
     ["contig",      undef   , "default_version,sequence_level", 2]);

  my @assembly_mappings =  ("scaffold:$ass_def|contig");


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

  $self->debug( "ApisMellifera Specific: creating contig and " .
               "scaffold seq_regions");

  $self->contig_to_seq_region('contig');
  $self->chromosome_to_seq_region("scaffold");
}


sub create_assembly {
  my $self = shift;

  $self->debug("ApisMellifera Specific: loading assembly data");

  $self->assembly_contig_chromosome();
}


1;
