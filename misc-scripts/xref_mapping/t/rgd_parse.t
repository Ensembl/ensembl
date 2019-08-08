use strict;
use warnings;
use Test::More;

use XrefParser::RGDParser;
use Xref::Test::TestDB;

my $db = Xref::Test::TestDB->new();
# Initialise two sources:
my $rgd_source = $db->create_db_row('Source',{
  name => 'RGD',
  status => 'KNOWN',
  priority_description => 'Test RGD source',
  priority => 10
});

my $refseq_source = $db->create_db_row('Source',{
  name => 'refseq',
  status => 'XREF',
  priority_description => 'Test RefSeq source for dependent xref creation',
  priority => 20
});

# Create a tame RefSeq xref to which the RGD test data refers
my $refseq_xref = $db->create_db_row('Xref',{
  accession => 'NM1234',
  version => 1,
  label => 'RefSeq innit?',
  source_id => $refseq_source->source_id,
  species_id => 34,
  info_type => 'COORDINATE_OVERLAP',
  info_text => '',
});


my $parser = XrefParser::RGDParser->new($db->dbh);

# Test without any prior RefSeq entries
$parser->run({ files => ['test-data/rgd.txt' ], verbose => 1, species_id => 34, source_id => $rgd_source->source_id });


ok ($db->schema->resultset('Xref')->check_direct_xref({
  accession => 1594427,
  label => '2331ex4-5',
  description => 'class I gene fragment 2331',
}), 'Sample rat Xref has been inserted');

cmp_ok($db->schema->resultset('Xref')->count(), '==', 4, 'Three xrefs in source file and one precursor are now in DB');

# Test that an RGD dependent xref has been assigned to $refseq_xref above

my $hit = $db->schema->resultset('DependentXref')->fetch_dependent_xref('NM1234','1500000');

ok($hit,'Found a dependent Xref');
is ($hit->dependent_xref->accession, '1500000', 'RGD xref is dependent on RefSeq accession');


done_testing();
