package Bio::EnsEMBL::PepDnaAlignFeature;

# EnsEMBL module for storing dna-dna pairwise alignments
#
# Cared for by Michele Clamp <michele@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself
#

=head1 NAME

  Bio::EnsEMBL::PepDnaAlignFeature - Ensembl specific pep-dna pairwise alignment feature

=head1 SYNOPSIS

  See BaseAlignFeature

=cut 

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
