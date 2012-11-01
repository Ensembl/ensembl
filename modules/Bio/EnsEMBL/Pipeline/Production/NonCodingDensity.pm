package Bio::EnsEMBL::Pipeline::Production::NonCodingDensity;

use base qw/Bio::EnsEMBL::Pipeline::Production::DensityGenerator/;


use strict;
use warnings;


sub get_option {
  my ($self) = @_;
  my @biotypes = $self->get_biotype_group("noncoding");
  return \@biotypes;
}

1;


