use lib 't';
use strict;

BEGIN { $| = 1;
	use Test ;
	plan tests => 12
}

use MultiTestDB;
use Bio::EnsEMBL::MiscSet;
use TestUtils qw(debug test_getter_setter);

our $verbose = 0; #set to 1 to turn on debug printouts


my $multi_db = MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');


#
# Test constructor
#

my $msa = $db->get_MiscSetAdaptor();


ok($msa && ref($msa) && $msa->isa('Bio::EnsEMBL::DBSQL::MiscSetAdaptor'));


#
# Test fetch_all
#

ok(@{$msa->fetch_all()} == 4);

#
# Test fetch_by_dbID
#

my $ms = $msa->fetch_by_dbID(1);

ok($ms->dbID(1));
ok($ms->code() eq 'ntctgs');
ok($ms->name() eq 'NT contigs');
ok($ms->description eq '');
ok($ms->longest_feature == 7e7);

#
# Test fetch_by_code
#

$ms = $msa->fetch_by_code('cloneset');
ok($ms->dbID() == 4);
ok($ms->code() eq 'cloneset');
ok($ms->name() eq '1Mb cloneset');
ok($ms->description eq '');
ok($ms->longest_feature == 6e6);
