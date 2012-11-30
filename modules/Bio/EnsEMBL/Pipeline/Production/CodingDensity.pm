package Bio::EnsEMBL::Pipeline::Production::CodingDensity;

use base qw/Bio::EnsEMBL::Pipeline::Production::DensityGenerator/;


use strict;
use warnings;


sub get_option {
  my ($self) = @_;
  return $self->get_biotype_group("coding");
}

1;


