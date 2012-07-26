package Bio::EnsEMBL::Pipeline::Production::PercentGC;

use base qw/Bio::EnsEMBL::Pipeline::Production::DensityGenerator/;


use strict;
use warnings;



sub get_density {
  my ($self, $block) = @_;
  my $gc = $block->get_base_count->{'%gc'};
  return $gc;
}

return 1;

