package Bio::EnsEMBL::Pipeline::Production::NonCodingDensity;

use base qw/Bio::EnsEMBL::Pipeline::Production::DensityGenerator/;


use strict;
use warnings;


sub get_option {
  my ($self) = @_;
  my @biotypes = $self->get_biotype_group("noncoding");
  return \@biotypes;
}

sub get_density {
  my ($self, $block, $biotypes) = @_;
  my $count = 0;
  foreach my $biotype (@$biotypes) {
    $count += scalar(@{ $block->get_all_Genes_by_type($biotype) });
  }
  return $count;
}

sub get_total {
  my ($self, $option) = @_;
  my $species = $self->param('species');
  my $ga = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'gene');
  return scalar(@{ $ga->fetch_all_by_biotype($option) });
}

1;


