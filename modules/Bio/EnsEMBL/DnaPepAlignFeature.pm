package Bio::EnsEMBL::DnaPepAlignFeature;

# EnsEMBL module for storing dna-dna pairwise alignments
#
# Cared for by Michele Clamp <michele@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself
#

=head1 NAME

  Bio::EnsEMBL::DnaPepAlignFeature - Ensembl specific dna-pep pairwise alignment feature

=head1 SYNOPSIS

  See BaseAlignFeature

=cut 


use Bio::EnsEMBL::BaseAlignFeature;

use vars qw(@ISA);
use strict;

@ISA = qw( Bio::EnsEMBL::BaseAlignFeature );


sub _hit_unit {
  return 1;
}

sub _query_unit {
  return 3;
}




1;
