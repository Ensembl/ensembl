use strict;
use warnings;

use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::OligoFeature;
use Bio::EnsEMBL::OligoProbe;
use Bio::EnsEMBL::OligoArray;

use Bio::EnsEMBL::Test::TestUtils qw( test_getter_setter debug );

BEGIN { $| = 1;
	use Test;
	plan tests => 14;
}


our $verbose = 0;


#
# Oligo arrays don't have to have anything more than a name,
# but it should be possible to have some other information on them
#

my $another_array = undef;

my $oligo_array = Bio::EnsEMBL::OligoArray->new 
  (
     -name => "Array-1",
     -setsize => 1,
     -included_in => $another_array,
	 -type => 'OLIGO',
);

ok( ref( $oligo_array ) && $oligo_array->isa( "Bio::EnsEMBL::OligoArray" ));

#
# What oligo arrays should be able to do
#

ok( test_getter_setter( $oligo_array, "dbID", 3 ));
ok( test_getter_setter( $oligo_array, "adaptor", undef ));
ok( test_getter_setter( $oligo_array, "name", "Some name" ));
ok( test_getter_setter( $oligo_array, "setsize", 11 ));
ok( test_getter_setter( $oligo_array, "type", 'OLIGO' ));

# possibly do the following with database connection
#my $probes = $oligo_array->get_all_Probes();

#
# Create an OligoProbe
#

my $oligo_probe = Bio::EnsEMBL::OligoProbe->new
(
  -arrayname => "Array-1",
  -name => "123-145",
  -description => 'A jolly interesting description',
  -probelength => 65,
);

$oligo_probe = Bio::EnsEMBL::OligoProbe->new
(
  -array => $oligo_array,
  -name => "123-145",
  -description => 'A jolly interesting description',
  -probelength => 50,
);

ok( ref( $oligo_probe ) && $oligo_probe->isa( "Bio::EnsEMBL::OligoProbe" ));



#
# What probes need to be able to do
#

ok( test_getter_setter( $oligo_probe, "dbID", 1 ));
ok( test_getter_setter( $oligo_probe, "adaptor", bless( {}, "Bio::EnsEMBL::DBSQL::BaseAdaptor" )));
ok( test_getter_setter( $oligo_probe, "description", "Some description" ));
ok( test_getter_setter( $oligo_probe, "probelength", 60 ));


my $arrays = $oligo_probe->get_all_Arrays();
ok( ref( $arrays ) eq "ARRAY" );
ok( $arrays->[0] && $arrays->[0]->isa( "Bio::EnsEMBL::OligoArray" ));


my $probenames = $oligo_probe->get_all_probenames();
my $probename = $oligo_probe->get_probename( $oligo_array->name() );


#
# Create an oligo feature
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
     
my $oligo_feature = Bio::EnsEMBL::OligoFeature->new
    (
     -probe => $oligo_probe,
     -mismatch_count => 1,
     -slice => $slice,
     -start => 1,
     -end => 10,
     -strand => -1,
);

ok( ref( $oligo_feature ) && $oligo_feature->isa( "Bio::EnsEMBL::OligoFeature"));

1;
