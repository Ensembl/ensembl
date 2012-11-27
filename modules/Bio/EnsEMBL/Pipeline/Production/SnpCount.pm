package Bio::EnsEMBL::Pipeline::Production::SnpCount;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Pipeline::Production::StatsGenerator/;



sub get_feature_count {
  my ($self, $slice, $key) = @_;
  my $variation_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($self->param('species'), 'variation');
  my $helper = $variation_adaptor->dbc()->sql_helper();
  my $sql = q{
     SELECT count(*) FROM variation_feature
     WHERE seq_region_id = ? };
  my @params = [$slice->get_seq_region_id];
  my $count = $helper->execute_single_result(-SQL => $sql, -PARAMS => @params);
  return $count;
}


sub get_total {
  my ($self) = @_;
  my $variation_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($self->param('species'), 'variation');
  my $helper = $variation_adaptor->dbc()->sql_helper();
  my $sql = q{
     SELECT count(*) FROM variation_feature };
  my $count = $helper->execute_single_result(-SQL => $sql);
  return $count;
}


sub get_slices {
  my ($self, $species) = @_;
  my @slices;
  my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'variation');
  my $sa = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'slice');
  my $helper = $dba->dbc()->sql_helper();
  my $sql = q{
    SELECT DISTINCT seq_region_id FROM variation_feature };
  my @ids = @{ $helper->execute_simple(-SQL => $sql) };
  foreach my $id(@ids) {
    push @slices, $sa->fetch_by_seq_region_id($id);
  }
  return \@slices;
}


sub get_attrib_codes {
  my ($self) = @_;
  my $prod_dba   = $self->get_production_DBAdaptor();
  my $prod_helper     = $prod_dba->dbc()->sql_helper();
  my $sql = q{
    SELECT code, name
    FROM attrib_type
    WHERE name = 'SNP count' };
  my %attrib_codes = %{ $prod_helper->execute_into_hash(-SQL => $sql) };
  return %attrib_codes;
}

# Blank as I don't think there are any alt attrib codes for this ATMO
sub get_alt_attrib_codes {
  my ($self) = @_;
  return;
}


1;

