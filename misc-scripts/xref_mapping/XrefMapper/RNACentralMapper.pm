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

package XrefMapper::RNACentralMapper;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw);

use base qw(XrefMapper::ChecksumMapper);

# Target file to use for checksum mapping
# RNACentral uses mRNA file
sub target {
  my ($self) = @_;
  return $self->mapper->core->dna_file();
}

sub external_db_name {
  my ($self) = @_;
  return 'RNAcentral';
}

# RNACentral xrefs are mapped on a transcript level
sub object_type {
  my ($self) = @_;
  return 'Transcript';
}

1;

