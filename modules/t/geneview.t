use strict;
use warnings;
use vars qw( $verbose );

BEGIN { $| = 1;  
	use Test;
	plan tests => 4;
}


use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

$verbose = 0;

ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

ok($multi);

my $db = $multi->get_DBAdaptor( "core" );

my $gene = $db->get_GeneAdaptor->fetch_by_transcript_stable_id( "ENST00000217347" );
my $geneid = $gene->stable_id;
ok( $geneid );

$gene = $db->get_GeneAdaptor->fetch_by_translation_stable_id( "ENSP00000278995" );
$geneid = $gene->stable_id;

$gene = $db->get_GeneAdaptor()->fetch_by_stable_id( "ENSG00000101321" );
ok( $gene );


