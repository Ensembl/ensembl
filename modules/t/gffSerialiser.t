## no critic (RequireFilenameMatchesPackage)
package Test::SO::Term;
use strict;
use warnings;

sub new {
  my ($class) = @_;
  return bless({}, ref($class) || $class);
}

sub name {
  my ($self) = @_;
  return 'feature';
}

package Test::SO;

use base qw/Bio::EnsEMBL::DBSQL::OntologyTermAdaptor/;

sub new {
  my ($class) = @_;
  return bless({}, ref($class) || $class);
}

sub fetch_by_accession {
  my ($self) = @_;
  return Test::SO::Term->new();
}

package main;

use strict;
use warnings;
use Test::More;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Utils::IO::GFFSerializer;
use IO::String;

my $db = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $db->get_DBAdaptor('core');

my $id = 'ENSG00000131044';

my $ga = $dba->get_GeneAdaptor();

{
  my $gene = $ga->fetch_by_stable_id($id);
  delete $gene->{source};
  $gene->{description} = undef; #empty value means don't emit the key
  my $expected = <<'OUT';
##gff-version 3
##sequence-region   20 30274334 30300924
OUT
  #Have to do this outside of the HERETO thanks to tabs
  $expected .= join("\t", 
    qw/20  ensembl feature 30274334  30300924  . + ./,
    'ID=ENSG00000131044;biotype=protein_coding;external_name=C20orf125;logic_name=ensembl' 
  );
  $expected .= "\n";

  assert_gff3($gene, $expected, 'Gene with no source serialises to GFF3 as expected. Source is ensembl');
}

{
  my $gene = $ga->fetch_by_stable_id($id);
  $gene->source('wibble');
  my $expected = <<'OUT';
##gff-version 3
##sequence-region   20 30274334 30300924
OUT
  #Have to do this outside of the HERETO thanks to tabs
  $expected .= join("\t", 
    qw/20  wibble feature 30274334  30300924  . + ./,
    'ID=ENSG00000131044;biotype=protein_coding;description=DJ310O13.1.2 (NOVEL PROTEIN SIMILAR DROSOPHILA PROTEIN CG7474%2C ISOFORM 2 ) (FRAGMENT). [Source:SPTREMBL%3BAcc:Q9BR18];external_name=C20orf125;logic_name=ensembl' 
  );
  $expected .= "\n";

  assert_gff3($gene, $expected, 'Gene with custom source serialises to GFF3 as expected. Source is wibble');
}

sub assert_gff3 {
  my ($feature, $expected, $msg) = @_;
  my $ota = Test::SO->new();
  my $fh = IO::String->new();
  my $ser = Bio::EnsEMBL::Utils::IO::GFFSerializer->new($ota, $fh);
  $ser->print_main_header([$feature->feature_Slice()]);
  $ser->print_feature($feature);
  is(${$fh->string_ref()}, $expected, $msg);
}

done_testing();