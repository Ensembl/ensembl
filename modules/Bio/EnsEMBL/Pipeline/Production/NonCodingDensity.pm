package Bio::EnsEMBL::Pipeline::Production::NonCodingDensity;

use base qw/Bio::EnsEMBL::Pipeline::Production::DensityGenerator/;


use strict;
use warnings;


sub get_option {
  my ($self) = @_;
  return $self->get_biotype_group("noncoding");
}

1;


