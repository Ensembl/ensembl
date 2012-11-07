use strict;
use warnings;

use Test::More;

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose = 1;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $multi->get_DBAdaptor('core');

$multi->save('core', 'meta_coord');


#
# 1 - Can construct meta coord container
#

my $mcc = $db->get_MetaCoordContainer();
ok($mcc);


my @coord_systems = @{$mcc->fetch_all_CoordSystems_by_feature_type('exon')};

ok(@coord_systems == 1);

ok($coord_systems[0]->name eq 'chromosome');

my $cs = $db->get_CoordSystemAdaptor->fetch_by_name('contig');

my $count = count_rows($db, 'meta_coord');

my $max = -1;
for( my $i=0; $i<10; $i++ ) {
  my $length = int(rand( 1000) + 100);
  $mcc->add_feature_type($cs, 'exon', $length );
  $max = $length if ( $length > $max );
}

my $length = $mcc->fetch_max_length_by_CoordSystem_feature_type( $cs, "exon" ); 
#debug( "max = $max; length=$length ");
ok( $max == $length );
ok(count_rows($db, 'meta_coord') == $count + 1);

@coord_systems = @{$mcc->fetch_all_CoordSystems_by_feature_type('exon')};

ok(@coord_systems == 2);

ok($coord_systems[0]->name eq 'chromosome');
ok($coord_systems[1]->name eq 'contig');


$multi->restore('core', 'meta_coord');

done_testing();