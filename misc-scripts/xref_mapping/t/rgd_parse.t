use strict;
use warnings;
use Test::More;

use XrefParser::RGDParser;
use Xref::Test::TestDB;

my $db = Xref::Test::TestDB->new();

my $parser = XrefParser::RGDParser->new($db->dbh);

# Test without any prior RefSeq entries
$parser->run({ files => ['test_data/rgd.txt' ], verbose => 1, species_id => 34, source_id => 2800 });

ok ($db->schema->resultset('Xref')->check_direct_xref({
  accession => 1594427,
  label => '2331ex4-5',
  description => 'class I gene fragment 2331',
}), 'Sample rat direct Xref has been inserted');


done_testing();