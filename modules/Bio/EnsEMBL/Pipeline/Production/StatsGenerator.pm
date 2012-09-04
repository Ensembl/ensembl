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
  my $total = $self->get_total();
  my $sum = 0;

  my $slices = $self->get_slices($species);
  my @sorted_slices = 
     sort( { $a->coord_system()->rank() <=> $b->coord_system()->rank()
             || $b->seq_region_length() <=> $a->seq_region_length() } @$slices) ;

  while (my $slice = shift @sorted_slices) {
    foreach my $code (keys %attrib_codes) {
      my $count = $self->get_feature_count($slice, $code, $attrib_codes{$code});
      if ($count > 0) {
        $self->store_attrib($slice, $count, $code);
      }
      $sum += $count;
    }
    if ($sum >= $total) {
      last;
    }
  }
}

sub get_slices {
  my ($self, $species) = @_;
  my @slices = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'slice')->fetch_all('toplevel');
  return \@slices;
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
  foreach my $code (keys %attrib_codes) {
    $helper->execute_update(-SQL => $sql, -PARAMS => [$dba->species_id(), $code]);
  }
}


sub store_attrib {
  my ($self, $slice, $count, $code) = @_;
  my $aa          = Bio::EnsEMBL::Registry->get_adaptor($self->param('species'), 'core', 'Attribute');
  my $prod_dba    = $self->get_production_DBAdaptor();
  my $prod_helper = $prod_dba->dbc()->sql_helper();
  my ($name, $description);
  my $sql = q{
    SELECT name
    FROM attrib_type
    WHERE code = ? };
  my @names = @{ $prod_helper->execute_simple(-SQL => $sql, -PARAMS => [$code]) };
  foreach my $bit (@names) {
    $name .= $bit . " ";
  }
  $sql = q{
    SELECT description
    FROM attrib_type
    WHERE code = ? };
  my @descriptions = @{ $prod_helper->execute_simple(-SQL => $sql, -PARAMS => [$code]) };
  foreach my $bit (@descriptions) {
    $description .= $bit . " ";
  }
  my $attrib = Bio::EnsEMBL::Attribute->new(
    -NAME        => $name,
    -CODE        => $code,
    -VALUE       => $count,
    -DESCRIPTION => $description
  );
  my @attribs = ($attrib);
  $aa->store_on_Slice($slice, \@attribs);
}

sub get_biotype_group {
  my ($self, $biotype) = @_;
  my $prod_dba = $self->get_production_DBAdaptor();
  my $helper = $prod_dba->dbc()->sql_helper();
  my $sql = q{
     SELECT name
     FROM biotype
     WHERE object_type = 'gene'
     AND is_current = 1
     AND biotype_group = ?
     AND db_type like '%core%' };
  my @biotypes = @{ $helper->execute_simple(-SQL => $sql, -PARAMS => [$biotype]) };
  return \@biotypes;
}



1;

