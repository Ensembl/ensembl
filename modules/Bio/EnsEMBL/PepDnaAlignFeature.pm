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


sub transform {
  my $self = shift;

  $self->throw( "PepDnaAlignFeatures cant be transformed as".
		" they are not on EnsEMBL coord system" );
}

sub _hit_unit {
  return 3;
}

sub _query_unit {
  return 1;
}

1;
