use strict;

use Test::More;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Map::MarkerFeature;
use Bio::EnsEMBL::Map::MarkerSynonym;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0; #set to 1 to turn on debug printouts

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $multi->get_DBAdaptor( 'core' );

######
# 1  #
######

#test constructor

my $adaptor = $db->get_MarkerFeatureAdaptor;
my $dbID = 111;
my $start = 100;
my $end   = 10;
my $slice = 
  $db->get_SliceAdaptor->fetch_by_region('contig','AL359765.6.1.13780');

my $analysis = Bio::EnsEMBL::Analysis->new;
my $marker_id = 1;
my $mapweight = 1;


#####
# 1 #
#####

#test construction
my $mf = Bio::EnsEMBL::Map::MarkerFeature->new
  ($dbID, $adaptor, $start, $end, $slice, $analysis, $marker_id, 
   $mapweight);

ok($mf && ref $mf && $mf->isa('Bio::EnsEMBL::Map::MarkerFeature'));


#######
# 2-3 #
#######

#test dbID
ok($dbID == $mf->dbID);
ok(&test_getter_setter($mf, 'dbID', undef));

#######
# 4-5 #
#######

#test adaptor
ok($adaptor == $mf->adaptor);
ok(&test_getter_setter($mf, 'adaptor', undef));

#######
# 6   #
#######

#test marker lazy-loading

my $marker = $mf->marker;
ok($marker->dbID == $marker_id);

#######
# 7-10#
#######

#test contig, start, end, strand (inherited)

ok($slice == $mf->slice);
ok($start == $mf->start);
ok($end == $mf->end);
ok($mf->strand == 0);

#######
#  11 #
#######

ok($mapweight == $mf->map_weight);



my $ms = Bio::EnsEMBL::Map::MarkerSynonym->new(1234, 'unists', 'a marker');

$mf->marker()->display_MarkerSynonym($ms);
ok($mf->display_id() eq $mf->marker()->display_MarkerSynonym()->name());

done_testing();
