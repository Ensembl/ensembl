use strict;

BEGIN { $| = 1;  
	use Test;
	plan tests => 9;
}

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Map::Ditag;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor( 'core' );

my $ditag;
my $dbID      = 1;
my $name      = "101A01-2";
my $type      = "ZZ11";
my $tag_count = 2;
my $sequence  = "GAGAACTTGGACCGCAGAGAATACACACAAATCAAACC";
my $adaptor   = $db->get_DitagAdaptor;
my $ditag_id  = 3278337;

######
# 1  #
######

#test new

$ditag = Bio::EnsEMBL::Map::Ditag->new (
					   -dbID      => $dbID,
                                           -name      => $name, 
                                           -type      => $type,
					   -tag_count => $tag_count,
                                           -sequence  => $sequence, 
                                           -adaptor   => $adaptor,
                                        );
ok($ditag && $ditag->isa('Bio::EnsEMBL::Map::Ditag'));

#######
# 2-6 #
#######

#test dbID, name, type, sequence, tag-count

ok($ditag->dbID eq $dbID);
ok($ditag->name eq $name);
ok($ditag->type eq $type);
ok($ditag->sequence() eq $sequence);
ok($ditag->tag_count > 0);

######
# 7  #
######

#test adaptor

ok($ditag->adaptor->isa('Bio::EnsEMBL::Map::DBSQL::DitagAdaptor'));

#######
# 8-9 #
#######

#test get_ditagFeatures

my $ditagFeatures = $adaptor->fetch_by_dbID($ditag_id)->get_ditagFeatures();
ok(scalar @$ditagFeatures);
ok($ditagFeatures->[0]->isa('Bio::EnsEMBL::Map::DitagFeature'));

1;
