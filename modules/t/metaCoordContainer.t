use strict;
use warnings;

use lib 't';

BEGIN { $| = 1;
	use Test;
	plan tests => 7;
}

use TestUtils qw( debug count_rows);

use MultiTestDB;


my $multi = MultiTestDB->new();
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

$mcc->add_feature_type($cs, 'exon');

ok(count_rows($db, 'meta_coord') == $count + 1);

@coord_systems = @{$mcc->fetch_all_CoordSystems_by_feature_type('exon')};

ok(@coord_systems == 2);

ok($coord_systems[0]->name eq 'chromosome');
ok($coord_systems[1]->name eq 'contig');


$multi->restore('core', 'meta_coord');

