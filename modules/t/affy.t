use strict;
use warnings;

use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::AffyFeature;
use Bio::EnsEMBL::AffyProbe;
use Bio::EnsEMBL::AffyArray;

use Bio::EnsEMBL::Test::TestUtils qw( test_getter_setter debug );

BEGIN { $| = 1;
	use Test;
	plan tests => 16;
}


our $verbose = 0;


#
# Affy arrays dont have to have anythong more than a name, 
# but it should be possible to have some other information on them
#

my $another_array = undef;

my $affy_array = Bio::EnsEMBL::AffyArray->new 
  (
     -name => "Affy-1",
     -setsize => 11,
     -included_in => $another_array,
     -probecount => 4443
);

ok( ref( $affy_array ) && $affy_array->isa( "Bio::EnsEMBL::AffyArray" ));

#
# What Affy arryas should be able to do
#

ok( test_getter_setter( $affy_array, "dbID", 3 ));
ok( test_getter_setter( $affy_array, "adaptor", undef ));
ok( test_getter_setter( $affy_array, "name", "Some name" ));
ok( test_getter_setter( $affy_array, "setsize", 11 ));

# possibly do the following with database connection
my $probes = $affy_array->get_all_AffyProbes();

#
# Create an AffyProbe
#

my $affy_probe = Bio::EnsEMBL::AffyProbe->new
(
  -arrayname => "AFFY-1",
  -name => "123-145",
  -probeset => "affy_probeset"
);

$affy_probe = Bio::EnsEMBL::AffyProbe->new
(
  -array => $affy_array,
  -name => "123-145",
  -probeset => "affy_probeset"
);

ok( ref( $affy_probe ) && $affy_probe->isa( "Bio::EnsEMBL::AffyProbe" ));

my $merge_probe = Bio::EnsEMBL::AffyProbe->new
(
  -arraynames => [ "Affy-1", "Affy-2", "Affy-3" ],
  -names => [ "123-145", "23,24,56", "someplace" ],
  -probeset => "affy_probeset"
);

ok( ref( $merge_probe ) && $merge_probe->isa( "Bio::EnsEMBL::AffyProbe" ));




#
# What probes need to be able to do
#

ok( test_getter_setter( $merge_probe, "dbID", 1 ));
ok( test_getter_setter( $merge_probe, "adaptor", bless( {}, "Bio::EnsEMBL::DBSQL::BaseAdaptor" )));
ok( test_getter_setter( $merge_probe, "probeset", "Affy_probeset" ));


my $arrays = $affy_probe->get_all_AffyArrays();
ok( ref( $arrays ) eq "ARRAY" );
ok( $arrays->[0] && $arrays->[0]->isa( "Bio::EnsEMBL::AffyArray" ));


# expected to construct full probenames from Chipname, setname and probename
my $full_probenames = $merge_probe->get_all_complete_names();
ok( ref( $full_probenames ) eq 'ARRAY' );
ok( $full_probenames->[0] && ( $full_probenames->[0] =~ /affy_probeset/ ));
 

my $full_probename = $merge_probe->get_complete_name( "Affy-1" );
ok( $full_probename =~ /affy_probeset/ && $full_probename =~ /Affy-1/ );

my $probenames = $merge_probe->get_all_probenames();
my $probename = $merge_probe->get_probename( $affy_array->name() );


#
# When we implement storing this should be implemented as well
# for starters we are not creating these objects, but load the database with probes ...

# $affy_probe->add_array_name( $affy_array, $probename );


#
# Create an affy feature
#

my $coord_system = Bio::EnsEMBL::CoordSystem->new
   ( -name => "chromosome",
     -version => '',
     -rank => 1 );

my $slice = Bio::EnsEMBL::Slice->new
   ( 
     -seq_region_name => '1',
     -coord_system => $coord_system,
     -start => 1,
     -end => 50,
     -seq_region_length => 200_000_000
);
     
my $affy_feature = Bio::EnsEMBL::AffyFeature->new
    (
     -probe => $affy_probe,
     -mismatch_count => 1,
     -slice => $slice,
     -start => 1,
     -end => 10,
     -strand => -1,
     -probeset => "affy_probeset"	
);

ok( ref( $affy_feature ) && $affy_feature->isa( "Bio::EnsEMBL::AffyFeature"));

1;
