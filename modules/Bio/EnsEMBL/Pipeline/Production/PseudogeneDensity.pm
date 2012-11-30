package Bio::EnsEMBL::Pipeline::Production::PseudogeneDensity;

use base qw/Bio::EnsEMBL::Pipeline::Production::DensityGenerator/;


use strict;
use warnings;


sub get_option {
  my ($self) = @_;
  return $self->get_biotype_group("pseudogene");
}


1;


