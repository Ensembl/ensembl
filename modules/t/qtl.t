use strict;
use warnings;
use Test::More;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Map::Marker;
use Bio::EnsEMBL::Map::MarkerSynonym;

use Bio::EnsEMBL::Map::DBSQL::QtlAdaptor;
use Bio::EnsEMBL::Map::DBSQL::QtlFeatureAdaptor;
use Bio::EnsEMBL::Slice;

use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0; #set to 1 to turn on debug printouts

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $multi->get_DBAdaptor( 'core' );
my $sa = $db->get_SliceAdaptor();



my $slice = $sa->fetch_by_region( 'chromosome', "20" );

my $qtl_features = $slice->get_all_QtlFeatures();

debug( "found ".scalar( @$qtl_features ) );
debug( join( "\n", map { "$_" } %{$qtl_features->[0]})."\n" );
debug( join( "\n", map { "$_" } %{$qtl_features->[1]})."\n" );

ok( scalar( @$qtl_features ) == 2 );

my $qtladptor = $db->get_QtlAdaptor();
my $traits = $qtladptor->list_traits();


debug( join( "\n", map {"Trait: $_"} @$traits ));
ok( scalar( @$traits ) == 2 );

my $qtls = $qtladptor->fetch_all_by_trait( $traits->[0] );
my $qf = $qtls->[0]->get_QtlFeature();

debug( join( "\n", %$qf ));
ok( $qf->isa( "Bio::EnsEMBL::Map::QtlFeature" ));

my $qtl = $qtls->[0];
ok($qtl->flank_marker_1->isa('Bio::EnsEMBL::Map::Marker'));
ok($qtl->flank_marker_2->isa('Bio::EnsEMBL::Map::Marker'));

my $synonyms = $qtl->get_synonyms;

ok($synonyms->{'rat genome database'} eq 'rqtl2');

done_testing();
