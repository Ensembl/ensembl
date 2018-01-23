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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

package Bio::EnsEMBL::IdMapping::StableIdGenerator::AedesAegypti;

# Package that implements incrementing and verification of Aedes aegypti
# stable IDs as used by the VectorBase project.

use strict;
use warnings;

use base qw(Bio::EnsEMBL::IdMapping::StableIdGenerator::EnsemblGeneric);

sub increment_stable_id {

  # This method will increment a stable ID.  For Aedes aegypti, it will
  # pick out the numerical part of the stable ID (no matter what type of
  # stable ID it is) and increment it by one.  It will then replace the
  # numerical part by the incremented value and return the new stable
  # ID.  The parsing of the stable ID is very naive.

  my ( $self, $stable_id ) = @_;

  if ( !$self->is_valid($stable_id) ) {
    throw("Unknown or missing stable ID: $stable_id.");
  }

  $stable_id =~ /^(\D*)(\d+)(\D*)/;

  my $number_as_string = "$2";
  my $number           = $2 + 1;
  $stable_id = sprintf(
    "%s" . sprintf( "%%0%dd", length($number_as_string) ) . "%s",
    $1, $number, $3 );

  return $stable_id;
}

sub is_valid {

  # A stable ID is a valid Aedes aegypti stable ID if it begins with the
  # character string "AAEL".

  my ( $self, $stable_id ) = @_;

  if ( !( defined($stable_id) && $stable_id =~ /^AAEL/ ) ) { return 0 }

  return 1;
}

1;
