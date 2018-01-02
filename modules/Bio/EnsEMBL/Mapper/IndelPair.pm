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

=head1 NAME

Bio::EnsEMBL::Mapper::IndelPair

=head1 SYNOPSIS

=head1 DESCRIPTION

Two regions mapped between different coordinate systems are each
represented by a Bio::EnsEMBL::Mapper::Unit and joined together as a
Bio::EnsEMBL::Mapper::Pair, when one of the regions is an indel.

=head1 METHODS

=cut

package Bio::EnsEMBL::Mapper::IndelPair;

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::Mapper::Pair);

sub new {
  my ($proto, @args) = @_;

  my $class = ref($proto) || $proto;

  my $self = $class->SUPER::new(@args);    # create the Pair object
  $self->{'indel'} = 1;                    # and add the Indel flag

  return $self;
}

1;
