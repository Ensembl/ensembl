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

Bio::EnsEMBL::PepDnaAlignFeature - Ensembl specific pep-dna pairwise
alignment feature

=head1 SYNOPSIS

  See BaseAlignFeature

=head1 METHODS

=cut 

package Bio::EnsEMBL::PepDnaAlignFeature;

use Bio::EnsEMBL::BaseAlignFeature;

use vars qw(@ISA);
use strict;

@ISA = qw( Bio::EnsEMBL::BaseAlignFeature );

=head2 transform 

  Arg [1]    : none
  Example    : none
  Description: Overwrites Bio:EnsEMBL:Feature->transform as
               to give error message
  Status     : Stable

=cut

sub transform {
  my $self = shift;

  $self->throw( "PepDnaAlignFeatures cant be transformed as".
		" they are not on EnsEMBL coord system" );
}

=head2 _hit_unit

  Arg [1]    : none
  Example    : none
  Description: PRIVATE implementation of abstract superclass method.  Returns
               3 as the 'unit' used for the hit sequence.
  Returntype : int
  Exceptions : none
  Caller     : Bio::EnsEMBL::BaseAlignFeature
  Status     : Stable


=cut

sub _hit_unit {
  return 3;
}

=head2 _query_unit

  Arg [1]    : none
  Example    : none
  Description: PRIVATE implementation of abstract superclass method.  Returns
               1 as the 'unit' used for the query sequence.
  Returntype : int
  Exceptions : none
  Caller     : Bio::EnsEMBL::BaseAlignFeature
  Status     : Stable


=cut

sub _query_unit {
  return 1;
}

1;
