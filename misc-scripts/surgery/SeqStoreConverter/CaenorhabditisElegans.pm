=head1 LICENSE

See the NOTICE file distributed with this work for additional information
regarding copyright ownership.

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

package SeqStoreConverter::CaenorhabditisElegans;

use vars qw(@ISA);

@ISA = qw(SeqStoreConverter::BasicConverter);


sub create_coord_systems {
  my $self = shift;

  $self->debug("Caenorhabditis Specific: creating " .
               "clone and chromosome coord systems");

  my $target = $self->target();
  my $dbh    = $self->dbh();

  my $ass_def = $self->get_default_assembly();

  my @coords = 
    (["chromosome" , $ass_def, "default_version"             , 1],
     ["clone",      undef   , "default_version,sequence_level", 2]);

  my @assembly_mappings =  ("chromosome:$ass_def|clone");

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

  $self->debug("CaenorhabditisElegans Specific: creating clone and " .
               "chromosome seq_regions");

  $self->contig_to_seq_region('clone');
  $self->chromosome_to_seq_region();
}


sub create_assembly {
  my $self = shift;

  $self->debug("CaenorhabditisElegans Specific: loading assembly data");

  $self->assembly_contig_chromosome();
}


#
# override the contig_to_seqregion method so that contigs are given clone
# names instead
#
sub contig_to_seq_region {
  my $self = shift;
  my $target_cs_name = shift;

  my $target = $self->target();
  my $source = $self->source();
  my $dbh     = $self->dbh();

  $target_cs_name ||= 'contig';

  $self->debug("CaenorhabditisElegans Specific: Transforming contigs " .
               "into $target_cs_name seq_regions");

  my $cs_id = $self->get_coord_system_id($target_cs_name);

  my $sth = $dbh->prepare
    ("INSERT INTO $target.seq_region " .
     "SELECT ctg.contig_id, CONCAT(cln.embl_acc, '.', cln.embl_version), " .
     "       $cs_id, ctg.length " .
     "FROM   $source.contig ctg, $source.clone cln " .
     "WHERE  ctg.clone_id = cln.clone_id");

  $sth->execute();
  $sth->finish();

  return;
}



1;
