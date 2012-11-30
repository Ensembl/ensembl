package Bio::EnsEMBL::Pipeline::Production::PepStats;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Pipeline::Base/;

use Bio::EnsEMBL::Attribute;

sub run {
  my ($self)  = @_;
  my $species = $self->param('species');
  my $dbtype  = $self->param('dbtype');
  my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, $dbtype);
  if ($dbtype =~ 'vega' || $dbtype =~ 'otherf') {
	my $core_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');
	$dba->dnadb($core_dba);
  }
  my $helper = $dba->dbc()->sql_helper();

  my @attrib_codes = $self->get_attrib_codes();
  $self->delete_old_attrib($dba, @attrib_codes);

  my $tmpfile = $self->param('tmpdir') . "/$$.pep";
  $self->dump_translation($dba, $tmpfile);
  my $results = $self->run_pepstats($tmpfile);

  my $attrib_types = $self->get_attrib_types();

  foreach my $translation (keys %$results) {
	$self->store_attrib($attrib_types, $translation, $results->{$translation});
  }
}

sub get_attrib_types {
  my ($self)       = @_;
  my $prod_helper  = $self->get_production_DBAdaptor()->dbc()->sql_helper();
  my $attrib_types = {};
  for my $row (
	@{$prod_helper->execute(
		-SQL => q{
    SELECT code,name,description
    FROM attrib_type
    WHERE code in ('NumResidues','MolecularWeight','AvgResWeight','Charge','IsoPoint') })})
  {
	$attrib_types->{$row->[0]} = {name        => $row->[1],
								  description => $row->[2]};
  }
  return $attrib_types;
}

my $key_names        = {};
my $key_descriptions = {};

sub store_attrib {
  my ($self, $attrib_types, $translation, $results) = @_;
  my $dbtype = $self->param('dbtype');
  my $aa = Bio::EnsEMBL::Registry->get_adaptor($self->param('species'), $dbtype, 'Attribute');
  my @attribs = ();
  foreach my $key (keys %$results) {
	my $value = $results->{$key};
	my $attrib = Bio::EnsEMBL::Attribute->new(-NAME        => $attrib_types->{$key}{name},
											  -CODE        => $key,
											  -VALUE       => $value,
											  -DESCRIPTION => $attrib_types->{$key}{description});
	push(@attribs, $attrib);
  } ## end foreach my $key (keys %$results)
  $aa->store_on_Translation($translation, \@attribs);
} ## end sub store_attrib

sub run_pepstats {
  my ($self, $tmpfile) = @_;
  my $PEPSTATS = $self->param('binpath') . '/bin/pepstats';
  open(OUT, "$PEPSTATS -filter < $tmpfile 2>&1 |");
  my @lines   = <OUT>;
  my $attribs = {};
  my $tid;
  close(OUT);
  foreach my $line (@lines) {

	if ($line =~ /PEPSTATS of ([^ ]+)/) {
	  $tid = $1;
	} elsif (defined $tid) {
	  if ($line =~ /^Molecular weight = (\S+)(\s+)Residues = (\d+).*/) {
		$attribs->{$tid}{'NumResidues'}     = $3;
		$attribs->{$tid}{'MolecularWeight'} = $1;
	  } elsif ($line =~ /^Average(\s+)(\S+)(\s+)(\S+)(\s+)=(\s+)(\S+)(\s+)(\S+)(\s+)=(\s+)(\S+)/) {
		$attribs->{$tid}{'AvgResWeight'} = $7;
		$attribs->{$tid}{'Charge'}       = $12;
	  } elsif ($line =~ /^Isoelectric(\s+)(\S+)(\s+)=(\s+)(\S+)/) {
		$attribs->{$tid}{'IsoPoint'} = $5;
	  }
	}
  }
  return $attribs;
} ## end sub run_pepstats

sub delete_old_attrib {
  my ($self, $dba, @attrib_codes) = @_;
  my $helper = $dba->dbc()->sql_helper();
  my $sql    = q{
     DELETE ta
     FROM translation_attrib ta, attrib_type at, translation tl, transcript tr, seq_region s, coord_system c
     WHERE at.attrib_type_id = ta.attrib_type_id
     AND ta.translation_id = tl.translation_id
     AND tl.transcript_id = tr.transcript_id
     AND tr.seq_region_id = s.seq_region_id
     AND s.coord_system_id = c.coord_system_id
     AND c.species_id = ?
     AND at.code = ? };
  foreach my $code (@attrib_codes) {
	$helper->execute_update(-SQL => $sql, -PARAMS => [$dba->species_id(), $code]);
  }
}

sub get_attrib_codes {
  my ($self)      = @_;
  my $prod_dba    = $self->get_production_DBAdaptor();
  my $prod_helper = $prod_dba->dbc()->sql_helper();
  my $sql         = q{
    SELECT code
    FROM attrib_type
    WHERE description = 'Pepstats attributes' };
  my @attrib_codes = @{$prod_helper->execute_simple(-SQL => $sql)};
  return @attrib_codes;
}

sub dump_translation {
  my ($self, $dba, $tmpfile) = @_;
  my $helper = $dba->dbc()->sql_helper();
  my $dbtype = $self->param('dbtype');
  my $ta     = Bio::EnsEMBL::Registry->get_adaptor($self->param('species'), $dbtype, 'translation');
  open(TMP, "> $tmpfile");
  my $sql = q{
    SELECT tl.translation_id
    FROM translation tl, transcript tr, seq_region s, coord_system cs
    WHERE tr.transcript_id = tl.transcript_id
    AND tr.seq_region_id = s.seq_region_id
    AND s.coord_system_id = cs.coord_system_id
    AND cs.species_id = ? };
  my @translation_ids = @{$helper->execute_simple(-SQL => $sql, -PARAMS => [$dba->species_id()])};

  for my $dbid (@translation_ids) {
	my $translation = $ta->fetch_by_dbID($dbid);
	my $peptide_seq = $translation->seq();
	if ($peptide_seq !~ /\n$/) {
	  $peptide_seq .= "\n";
	}
	print TMP ">$dbid\n$peptide_seq";
  }
  close(TMP);
} ## end sub dump_translation

1;
