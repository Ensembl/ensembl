use lib 't';
use strict;

BEGIN { $| = 1;
	use Test ;
	plan tests => 2;
}

use MultiTestDB;
use TestUtils qw(debug test_getter_setter);

our $verbose = 0; #set to 1 to turn on debug printouts

my $multi = MultiTestDB->new();
my $db = $multi->get_DBAdaptor( 'core' );

#
# 1 Test AssemblyMapperAdaptor constructor
#
my $asma = $db->get_AssemblyMapperAdaptor();


ok($asma && $asma->isa('Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor'));


#
# 2 Test fetch_by_CoordSystems
#

my $csa = $db->get_CoordSystemAdaptor();
my $ctg_cs = $csa->fetch_by_name('contig');
my $chr_cs = $csa->fetch_by_name('chromosome');

my $asm_mapper = $asma->fetch_by_CoordSystems($ctg_cs, $chr_cs);

ok($asm_mapper && $asm_mapper->isa('Bio::EnsEMBL::AssemblyMapper'));


