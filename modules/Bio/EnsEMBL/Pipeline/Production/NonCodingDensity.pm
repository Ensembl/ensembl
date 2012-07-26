package Bio::EnsEMBL::Pipeline::Production::NonCodingDensity;

use base qw/Bio::EnsEMBL::Pipeline::Production::DensityGenerator/;


use strict;
use warnings;


sub get_density {
  my ($self, $block) = @_;
  my @biotypes = $self->get_biotype_group("non-coding");
  my $count = 0;
  foreach my $biotype (@biotypes) {
    $count += scalar(@{ $block->get_all_Genes_by_type($biotype) });
  }
  return $count;
}

1;


