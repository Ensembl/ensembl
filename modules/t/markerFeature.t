use lib 't';
use strict;

BEGIN { $| = 1;  
	use Test ;
	plan tests => 11
}

use MultiTestDB;
use Bio::EnsEMBL::Map::MarkerFeature;
use Bio::EnsEMBL::RawContig;
use Bio::EnsEMBL::Analysis;
use TestUtils qw(debug test_getter_setter);

our $verbose = 0; #set to 1 to turn on debug printouts

my $multi = MultiTestDB->new();
my $db = $multi->get_DBAdaptor( 'core' );

######
# 1  #
######

#test constructor

my $adaptor = $db->get_MarkerFeatureAdaptor;
my $dbID = 111;
my $start = 100;
my $end   = 10;
my $contig = Bio::EnsEMBL::RawContig->new;
my $analysis = Bio::EnsEMBL::Analysis->new;
my $marker_id = 1;
my $mapweight = 1;


#####
# 1 #
#####

#test construction
my $mf = Bio::EnsEMBL::Map::MarkerFeature->new
  ($dbID, $adaptor, $start, $end, $contig, $analysis, $marker_id, 
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

ok($contig == $mf->contig);
ok($start == $mf->start);
ok($end == $mf->end);
ok($mf->strand == 0);

#######
#  11 #
#######

ok($mapweight == $mf->map_weight);
