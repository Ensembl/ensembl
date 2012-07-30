package Bio::EnsEMBL::Pipeline::Production::GeneCount;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Pipeline::Production::StatsGenerator/;


sub get_attrib_codes {
  my ($self) = @_;
  my @attrib_codes = ('coding_cnt', 'pseudogene_cnt', 'noncoding_cnt');
  return @attrib_codes;
}


sub get_feature_count {
  my ($self, $slice, $key) = @_;
  my $prod_dba    = $self->get_production_DBAdaptor();
  my $prod_helper = $prod_dba->dbc()->sql_helper();
  my ($group) = $key =~ /(\w+)\_cnt/;
  my $sql = q{
    SELECT name
    FROM biotype
    WHERE biotype_group = ?
    AND object_type = 'gene'
    AND is_current = 1
    AND db_type like '%core%' };
  my @biotypes = @{ $prod_helper->execute_simple(-SQL => $sql, -PARAMS => [$group]) };
  my $count = 0;
  foreach my $biotype (@biotypes) {
    $count += scalar(@{ $slice->get_all_Genes_by_type($biotype) });
  }
  return $count;
}


1;

