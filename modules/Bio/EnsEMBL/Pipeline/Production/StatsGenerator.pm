package Bio::EnsEMBL::Pipeline::Production::StatsGenerator;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Pipeline::Base/;


use Bio::EnsEMBL::Attribute;

sub run {
  my ($self) = @_;
  my $species    = $self->param('species');
  my $dba        = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');

  my %attrib_codes = $self->get_attrib_codes();
  $self->delete_old_attrib($dba, %attrib_codes);

  my $slices = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'slice')->fetch_all('toplevel');
  while (my $slice = shift @$slices) {
    foreach my $key (keys %attrib_codes) {
      my $count = $self->get_feature_count($slice, $key);
      if ($count > 0) {
        $self->store_attrib($slice, $count, $attrib_codes{$key});
      }
    }
  }
}


sub delete_old_attrib {
  my ($self, $dba, %attrib_codes) = @_;
  my $helper = $dba->dbc()->sql_helper();
  my $sql = q{
    DELETE sa
    FROM seq_region_attrib sa, attrib_type at, seq_region s, coord_system cs
    WHERE s.seq_region_id = sa.seq_region_id
    AND cs.coord_system_id = s.coord_system_id
    AND at.attrib_type_id = sa.attrib_type_id
    AND cs.species_id = ?
    AND at.code = ? };
  foreach my $key (keys %attrib_codes) {
    $helper->execute_update(-SQL => $sql, -PARAMS => [$dba->species_id(), $attrib_codes{$key}]);
  }
}


sub store_attrib {
  my ($self, $slice, $count, $code) = @_;
  my $aa          = Bio::EnsEMBL::Registry->get_adaptor($self->param('species'), 'core', 'Attribute');
  my $prod_dba    = $self->get_production_DBAdaptor();
  my $prod_helper = $prod_dba->dbc()->sql_helper();
  my $sql = q{
    SELECT name, description
    FROM attrib_type
    WHERE code = ? };
  my ($name, $description) = $prod_helper->execute(-SQL => $sql, -PARAMS => [$code]);
  my $attrib = Bio::EnsEMBL::Attribute->new(
    -NAME        => $name,
    -CODE        => $code,
    -VALUE       => $count,
    -DESCRIPTION => $description
  );
  my @attribs = ($attrib);
  $aa->store_on_Slice($slice, \@attribs);
}


sub get_production_DBAdaptor {
  my ($self) = @_;
  return Bio::EnsEMBL::Registry->get_DBAdaptor('multi', 'production');
}


1;

