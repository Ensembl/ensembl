use strict;
use warnings;

use lib 't';

BEGIN { $| = 1;  
	use Test;
	plan tests => 9;
}

use TestUtils qw( debug );

use MultiTestDB;


my $mdb = MultiTestDB->new();
my $db = $mdb->get_DBAdaptor('core');

$mdb->save('core', 'meta');


#
# 1 - Can construct meta container
#

my $mc = $db->get_MetaContainer();
ok($mc);


#
# 2 - list_value_by_key
#

my ($asm_default) = @{$mc->list_value_by_key('assembly.default')};
ok($asm_default eq 'NCBI_30');


#
# 3-4 store key value
#

$mc->store_key_value('testkey', 'testvalue1');
$mc->store_key_value('testkey', 'testvalue2');

my $listref = $mc->list_value_by_key('testkey');
ok($listref->[0] eq 'testvalue1');
ok($listref->[1] eq 'testvalue2');


#
# 5-6 - get_Species
#

my $species = $mc->get_Species();
ok($species->common_name eq 'Human');
my $bin = $species->binomial;
ok($bin eq 'Homo sapiens');

#
# 7 - get_taxon_id
#
my $taxid = $mc->get_taxonomy_id();
ok($taxid == 9606);


#
# 8 - get_default_assembly
#
$asm_default = $mc->get_default_assembly();
ok($asm_default eq 'NCBI_30');

#
# 9 - get_max_assembly_contig
#
my $maxac = $mc->get_max_assembly_contig();
ok($maxac == 1_000_000);

$mdb->restore('core', 'meta');



#
# 10 - get_all_coord_systems
#
my @cs = grep {$_ eq 'chromosome' || $_ eq 'clone' || $_ eq 'contig'}
         @{$mc->get_all_coord_systems()};

ok(@cs == 3);


#
# 11 - get_top_coord_system
#
ok($mc->get_top_coord_system eq 'chromosome');

#
# 12-13 - is_valid_coord_system
#
ok($mc->is_valid_coord_system('chromosome'));
ok(!$mc->is_valid_coord_system('non-existant'));