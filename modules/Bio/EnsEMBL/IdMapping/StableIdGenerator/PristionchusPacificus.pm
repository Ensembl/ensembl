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

  Please email comments or questions to the public WormBase
  help list at <help@wormbase.org>.

=cut

package Bio::EnsEMBL::IdMapping::StableIdGenerator::PristionchusPacificus;

# Class to generate WormBase conform Pristionchus IDs to be injected into the mapper

use strict;
use warnings;

use base qw(Bio::EnsEMBL::IdMapping::StableIdGenerator::EnsemblGeneric);

# PPAxxxxx
# ups the id by one
sub increment_stable_id {
  my ( $self, $lastId ) = @_;

  throw("invalid stable ID: $lastId.") unless ($lastId=~/PPA/);

  $lastId =~ /^PPA(\d+)/;

  my $number = $1+1;
  my $stable_id = sprintf("PPA%05d",$number);

  return $stable_id;
}

# just in case it is actually used somewhere
sub is_valid {
  my ( $self, $stableId ) = @_;

  if ( defined($stableId) && $stableId =~ /^PPA\d+/)
     {return 1} else {return undef}
}

1;
