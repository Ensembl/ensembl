package Bio::EnsEMBL::Pipeline::Production::PercentGC;

use base qw/Bio::EnsEMBL::Pipeline::Production::DensityGenerator/;


use strict;
use warnings;



sub get_density {
  my ($self, $block) = @_;
  my $gc = $block->get_base_count->{'%gc'};
  return $gc;
}

sub get_total {
  my ($self) = @_;
  my $species = $self->param('species');
  my $slices = scalar(@{  Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'slice')->fetch_all('toplevel') });
  return $slices*$self->param('bin_count')*100;
}

return 1;

