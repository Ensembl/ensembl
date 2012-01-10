use strict;
use warnings;


BEGIN { $| = 1;
	use Test;
	plan tests => 8;
}

use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Test::MultiTestDB;


my $mdb = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $mdb->get_DBAdaptor('core');

$mdb->save('core', 'meta');


#
# 1 - Can construct meta container
#

my $mc = $db->get_MetaContainer();
ok($mc);


#
# list_value_by_key
#

my ($asm_default) = @{$mc->list_value_by_key('assembly.default')};
ok($asm_default eq 'NCBI34');


#
#  store key value
#

$mc->store_key_value('testkey', 'testvalue1');
$mc->store_key_value('testkey', 'testvalue2');

my $listref = $mc->list_value_by_key('testkey');
ok($listref->[0] eq 'testvalue1');
ok($listref->[1] eq 'testvalue2');

$mc->delete_key('testkey');

$listref = $mc->list_value_by_key('testkey');
ok(@$listref == 0);

ok($mc->get_common_name() eq 'Human');
my $bin = $mc->get_scientific_name();
ok($bin eq 'Homo sapiens');

#
# 7 - get_taxon_id
#
my $taxid = $mc->get_taxonomy_id();
ok($taxid == 9606);

$mdb->restore('core', 'meta');

