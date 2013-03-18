package Bio::EnsEMBL::Pipeline::Production::GeneGCBatch;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Pipeline::Production::GeneGC/;

# modified version of GeneGC that uses the faster store_batch_on_Gene method from AttributeAdaptor
sub run {
  my ($self) = @_;
  my $species = $self->param('species');
  my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');

  my $attrib_code = 'GeneGC';
  $self->delete_old_attrib($dba, $attrib_code);

  my $genes = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'gene')->fetch_all();
  my $aa = Bio::EnsEMBL::Registry->get_adaptor($self->param('species'), 'core', 'Attribute');

  my $prod_helper = $self->get_production_DBAdaptor()->dbc()->sql_helper();
  my ($name, $description) = @{
	$prod_helper->execute(
	  -SQL => q{
    SELECT code, name, description
    FROM attrib_type
    WHERE code = ? },
	  -PARAMS => [$attrib_code])->[0]};

  my $attributes = {};
  while (my $gene = shift @$genes) {
	my $count = $gene->feature_Slice()->get_base_count->{'%gc'};
	if ($count > 0) {
	  push @{$attributes->{$gene->dbID()}},
		Bio::EnsEMBL::Attribute->new(-NAME        => $name,
									 -CODE        => $attrib_code,
									 -VALUE       => $count,
									 -DESCRIPTION => $description);
	}
  }
  $aa->store_batch_on_Gene($attributes);

} ## end sub run

sub delete_old_attrib {
  my ($self, $dba, $attrib_code) = @_;
  my $helper = $dba->dbc()->sql_helper();
  my $sql    = q{
    DELETE ga
    FROM gene_attrib ga, attrib_type at, gene g, seq_region s, coord_system cs
    WHERE s.seq_region_id = g.seq_region_id
    AND g.gene_id = ga.gene_id
    AND cs.coord_system_id = s.coord_system_id
    AND at.attrib_type_id = ga.attrib_type_id
    AND cs.species_id = ?
    AND at.code = ? };
  $helper->execute_update(-SQL => $sql, -PARAMS => [$dba->species_id(), $attrib_code]);
}

1;

