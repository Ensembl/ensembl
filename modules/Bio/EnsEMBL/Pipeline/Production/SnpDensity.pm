package Bio::EnsEMBL::Pipeline::Production::SnpDensity;

use base qw/Bio::EnsEMBL::Pipeline::Production::DensityGenerator/;


use strict;
use warnings;



sub get_density {
  my ($self, $block) = @_;
  my $variation_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($self->param('species'), 'variation');
  my $helper = $variation_adaptor->dbc()->sql_helper();
  my $sql = q{
     SELECT count(*) FROM variation_feature
     WHERE seq_region_id = ?
     AND seq_region_start <= ?
     AND seq_region_end >= ? };
  my @params = [$block->get_seq_region_id, $block->end, $block->start];
  my $count = $helper->execute_single_result(-SQL => $sql, -PARAMS => @params);
  return $count;
}

sub get_total {
  my ($self, $option) = @_;
  my $species = $self->param('species');
  my $variation_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'variation');
  my $helper = $variation_adaptor->dbc()->sql_helper();
  my $sql = "SELECT count(*) FROM variation_feature";
  return $helper->execute_single_result(-SQL => $sql);
}

return 1;

