use lib 't';
use strict;

BEGIN { $| = 1;
	use Test ;
	plan tests => 16;
}

use MultiTestDB;
use TestUtils qw(debug test_getter_setter);

our $verbose = 0; #set to 1 to turn on debug printouts

my $multi = MultiTestDB->new();
my $db = $multi->get_DBAdaptor( 'core' );


#
# Test constructor
#

my $csa = $db->get_CoordSystemAdaptor();

ok($csa && $csa->isa('Bio::EnsEMBL::DBSQL::CoordSystemAdaptor'));


#
# 2-6 Test fetch_by_name()
#

my $cs = $csa->fetch_by_name('chromosome');

ok($cs->name eq 'chromosome');
ok($cs->dbID());
ok($cs->version eq 'NCBI33');
ok($cs->is_top_level());
ok(!$cs->is_sequence_level());


#
# 7-8 Test fetch_all_by_name
#
my @cs_list = @{$csa->fetch_all_by_name('chromosome')};

ok(@cs_list == 1);

ok($cs_list[0]->equals($cs));


#
# 9-10 Test fetch_by_dbID()
#
$cs = $csa->fetch_by_dbID(3);

ok($cs->name() eq 'clone');
ok($cs->version() eq '');



#
# 11 Test fetch_top_level
#
$cs = $csa->fetch_top_level();

ok($cs->name eq 'chromosome');

#
# 12 Test fetch_all_top_level
#
($cs) = @{$csa->fetch_all_top_level()};

ok($cs->name eq 'chromosome');


#
# 13-14 Test fetch_sequence_level
#
my $cs = $csa->fetch_sequence_level();

ok($cs->name eq 'contig');
ok($cs->is_sequence_level());


#
# 15-16 Test get_mapping_path
#

my $ctg_cs = $csa->fetch_by_name('contig');
my $chr_cs = $csa->fetch_by_name('chromosome');
my $cln_cs = $csa->fetch_by_name('clone');

my $path = $csa->get_mapping_path($ctg_cs, $chr_cs);

ok(@$path == 2 &&
   $path->[0]->name() eq 'chromosome' &&
   $path->[1]->name() eq 'contig');

$path = $csa->get_mapping_path($chr_cs, $cln_cs);


ok(@$path == 3 &&
   (($path->[0]->name eq 'chromosome' &&  #there are 2 equally valid paths
     $path->[1]->name eq 'contig' &&
     $path->[2]->name eq 'clone')
    ||
    ($path->[0]->name eq 'clone' &&
     $path->[1]->name eq 'contig' &&
     $path->[2]->name eq 'chromosome'))
   );














