use strict;
use warnings;

use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::DBSQL::AffyFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::AffyProbeAdaptor;
use Bio::EnsEMBL::DBSQL::AffyArrayAdaptor;
use Bio::EnsEMBL::AffyProbe;
use Bio::EnsEMBL::AffyArray;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Attribute;

use Bio::EnsEMBL::Test::TestUtils qw( test_getter_setter debug count_rows);

BEGIN { $| = 1;
	use Test;
	plan tests => 21;
}


our $verbose = 0;


my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;

# get a core DBAdaptor
my $db = $multi->get_DBAdaptor("core");

$multi->hide( "core", "oligo_array", "oligo_probe", "oligo_feature" ); 

my $slice_adaptor = $db->get_SliceAdaptor();
my $affyFeatureAdaptor   = $db->get_AffyFeatureAdaptor();
my $slice = $slice_adaptor->fetch_by_region( "chromosome", 20, 30_000_000, 31_000_000 );

ok( ref( $affyFeatureAdaptor ) && 
    $affyFeatureAdaptor->isa( "Bio::EnsEMBL::DBSQL::AffyFeatureAdaptor" ));

my $affyProbeAdaptor = $db->get_AffyProbeAdaptor();

ok( ref( $affyProbeAdaptor ) && 
    $affyProbeAdaptor->isa( "Bio::EnsEMBL::DBSQL::AffyProbeAdaptor" ));

my $affyArrayAdaptor = $db->get_AffyArrayAdaptor();

ok( ref( $affyArrayAdaptor ) && 
    $affyArrayAdaptor->isa( "Bio::EnsEMBL::DBSQL::AffyArrayAdaptor" ));

# the full circle, make array and probe
# store them, make Feature(s) store them
# retrieve the lot

# array creation and storing
my $array1 = Bio::EnsEMBL::AffyArray->new( 
    '-name' => 'Affy-1',
    -setsize => 12
);

my $array2 = Bio::EnsEMBL::AffyArray->new( 
    '-name' => 'Affy-2',
    -setsize => 12,
    -included_in => $array1
);

$affyArrayAdaptor->store( $array1 );
$affyArrayAdaptor->store( $array2 );

ok( $array1->dbID() );
ok( $array2->dbID() );

# probe creation and storing
my $probe = Bio::EnsEMBL::AffyProbe->new
    ( -probenames => [ '125:130', '222:242' ],
      -arrays => [ $array1, $array2 ],
      -probeset => 'arnes_12112'
    );

$affyProbeAdaptor->store( $probe );

ok( $probe->dbID() );
ok( count_rows( $affyProbeAdaptor, "oligo_probe" ) == 2 );

#
# Now a feature on this probe
#

# first need to create an analysis or
#  use one already there :-)

my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name( "Vertrna" );

my $feature = Bio::EnsEMBL::AffyFeature->new
    ( -slice => $slice,
      -start => 100,
      -end => 124,
      -strand => -1,
      -mismatches => 0,
      -probe => $probe,
      -analysis => $analysis
    );

$affyFeatureAdaptor->store( $feature );

ok( $feature->dbID() );
ok( count_rows( $affyFeatureAdaptor, "oligo_feature" ) == 1 );

#
# Now getting the stuff from the database
#

# getting features

my $db_feat = $slice->get_all_AffyFeatures();

ok( scalar @$db_feat == 1);

$db_feat = $slice->get_all_AffyFeatures( "Somename" );

ok( scalar @$db_feat == 0);

$db_feat = $slice->get_all_AffyFeatures( "Affy-1" );

ok( scalar @$db_feat == 1);

$db_feat = $slice->get_all_AffyFeatures( "Affy-1", "Affy-2" );

ok( scalar @$db_feat == 1);

$db_feat = $affyFeatureAdaptor->fetch_all_by_probeset( 'arnes_12112' );

ok( scalar @$db_feat == 1);

$db_feat = $affyFeatureAdaptor->fetch_all_by_probeset( 'some probeset' );

ok( scalar @$db_feat == 0);

$db_feat = $slice->get_all_AffyFeatures( "Affy-2" );

ok( scalar @$db_feat == 1);


#
# getting probes
#

$probe = $db_feat->[0]->probe();

ok( $probe->isa( "Bio::EnsEMBL::AffyProbe" ));

my $names = $probe->get_all_complete_names();
debug( "Name 1".$names->[0] );
debug( "Name 2".$names->[1] );

ok( $names->[0] =~ /Affy-1/ || $names->[1] =~ /Affy-1/ );
ok( $names->[0] =~ /Affy-2/ || $names->[1] =~ /Affy-2/ );

my $arrays = $affyArrayAdaptor->fetch_all();

ok( @$arrays == 2 );

#
# And the Arrays for the probe
#

$arrays = $probe->get_all_AffyArrays();

ok( @$arrays == 2 );

$multi->restore();
