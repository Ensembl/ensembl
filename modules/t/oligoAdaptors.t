use strict;
use warnings;

use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::DBSQL::OligoFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::OligoProbeAdaptor;
use Bio::EnsEMBL::DBSQL::OligoArrayAdaptor;
use Bio::EnsEMBL::OligoProbe;
use Bio::EnsEMBL::OligoArray;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Attribute;

use Bio::EnsEMBL::Test::TestUtils qw( test_getter_setter debug count_rows);

BEGIN { $| = 1;
	use Test;
	plan tests => 22;
}


our $verbose = 0;


my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;

# get a core DBAdaptor
my $db = $multi->get_DBAdaptor("core");

$multi->hide( "core", "oligo_array", "oligo_probe", "oligo_feature" ); 

my $slice_adaptor = $db->get_SliceAdaptor();
my $oligoFeatureAdaptor   = $db->get_OligoFeatureAdaptor();
my $slice = $slice_adaptor->fetch_by_region( "chromosome", 20, 30_000_000, 31_000_000 );

ok( ref( $oligoFeatureAdaptor ) && 
    $oligoFeatureAdaptor->isa( "Bio::EnsEMBL::DBSQL::OligoFeatureAdaptor" ));

my $oligoProbeAdaptor = $db->get_OligoProbeAdaptor();

ok( ref( $oligoProbeAdaptor ) && 
    $oligoProbeAdaptor->isa( "Bio::EnsEMBL::DBSQL::OligoProbeAdaptor" ));

my $oligoArrayAdaptor = $db->get_OligoArrayAdaptor();

ok( ref( $oligoArrayAdaptor ) && 
    $oligoArrayAdaptor->isa( "Bio::EnsEMBL::DBSQL::OligoArrayAdaptor" ));

# the full circle, make array and probe
# store them, make Feature(s) store them
# retrieve the lot

# array creation and storing
my $array1 = Bio::EnsEMBL::OligoArray->new( 
    '-name' => 'Array-1',
    -setsize => 1,
	-type => 'OLIGO',
);

my $array2 = Bio::EnsEMBL::OligoArray->new( 
    '-name' => 'Array-2',
    -setsize => 1,
	-type => 'OLIGO',
    -included_in => $array1,
);

my $array3 = Bio::EnsEMBL::OligoArray->new( 
    '-name' => 'Affy-3',
    -setsize => 16,
	-type => 'AFFY',
);

$oligoArrayAdaptor->store( $array1 );
$oligoArrayAdaptor->store( $array2 );
$oligoArrayAdaptor->store( $array3 );

ok( $array1->dbID() );
ok( $array2->dbID() );
ok( $array3->dbID() );

# probe creation and storing
my $probe = Bio::EnsEMBL::OligoProbe->new
    ( -probename => 'Probe-1',
      -array => $array1,
      -description => 'Probe-1 description',
	  -probelength => 65,
    );

$oligoProbeAdaptor->store( $probe );

ok( $probe->dbID() );
ok( count_rows( $oligoProbeAdaptor, "oligo_probe" ) == 1 );

#
# Now a feature on this probe
#

# first need to create an analysis or
#  use one already there :-)

my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name( "Vertrna" );

my $feature = Bio::EnsEMBL::OligoFeature->new
    ( -slice => $slice,
      -start => 100,
      -end => 124,
      -strand => -1,
      -mismatches => 0,
      -probe => $probe,
      -analysis => $analysis
    );

$oligoFeatureAdaptor->store( $feature );

ok( $feature->dbID() );
ok( count_rows( $oligoFeatureAdaptor, "oligo_feature" ) == 1 );

#
# Now getting the stuff from the database
#

# getting features

my $db_feat = $slice->get_all_OligoFeatures();

ok( scalar @$db_feat == 1);

$db_feat = $slice->get_all_OligoFeatures_by_type( 'OLIGO' );

ok( scalar @$db_feat == 1);

$db_feat = $slice->get_all_OligoFeatures_by_type( 'OLIGO', 'some logic name' );

ok( scalar @$db_feat == 0);

$db_feat = $slice->get_all_OligoFeatures_by_type( 'AFFY' );

ok( scalar @$db_feat == 0);

$db_feat = $slice->get_all_OligoFeatures( "Somename" );

ok( scalar @$db_feat == 0);

$db_feat = $slice->get_all_OligoFeatures( "Array-1" );

ok( scalar @$db_feat == 1);


#
# getting probes
#

$probe = $db_feat->[0]->probe();

ok( $probe->isa( "Bio::EnsEMBL::OligoProbe" ));

ok( $probe->probelength() == 65 );

my $arrays = $oligoArrayAdaptor->fetch_all();

ok( @$arrays == 3 );

#
# And the Arrays for the probe
#

$arrays = $probe->get_all_Arrays();

ok( @$arrays == 1 );

#
# Get all the OLIGO arrays and then all the AFFY and OLIGO arrays
#

$arrays = $oligoArrayAdaptor->fetch_all_by_type('OLIGO');

ok( @$arrays == 2 );

$arrays = $oligoArrayAdaptor->fetch_all_by_type('AFFY','OLIGO');

ok( @$arrays == 3 );

$multi->restore();

