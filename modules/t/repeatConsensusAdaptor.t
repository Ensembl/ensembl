use strict;
use lib 't';
use TestUtils qw(test_getter_setter);

use MultiTestDB;


BEGIN { $| = 1;
	use Test;
	plan tests => 22;
}

my $multi_db = MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');

my $verbose = 0;

# Test Creation

my $rca = $db->get_RepeatConsensusAdaptor();

ok($rca && ref($rca) && $rca->isa('Bio::EnsEMBL::RepeatConsensusAdaptor'));

#
# Test fetch_by_dbID
#

my $rc = $rca->fetch_by_dbID(9);
ok($rc->name() eq 'MIR3');
ok($rc->dbID == 9);
ok($rc->repeat_consensus eq '');
ok($rc->length() == 0);
ok($rc->repeat_class eq 'Type I Transposons/SINE');

#
# Test fetch_by_name
#
my $rc = $rca->fetch_by_name('MIR');
ok($rc->name() eq 'MIR');
ok($rc->dbID() == 1);
ok($rc->repeat_consensus eq '');
ok($rc->length() == 0);
ok($rc->repeat_class eq 'Type I Transposons/SINE');

#
# Test fetch_by_name_class
#

$rc = $rca->fetch_by_name_class('MER65A', 'LTRs');
ok($rc->name() eq 'MER65A');
ok($rc->dbID() == 283);
ok($rc->repeat_class eq 'LTRs');
ok($rc->repeat_consensus eq '');
ok($rc->length() == 0);

#
# Test fetch_all_by_class_seq
#
ok(@{$rca->fetch_all_by_class_seq('LTRs', '')} == 38);

#
# Test store
#

$multi_db->save('core', 'repeat_consensus');

my $rc = Bio::EnsEMBL::RepeatConsensus->new
  (-REPEAT_CONSENSUS => 'ACTG',
   -NAME             => 'ACTG(n)',
   -LENGTH           => 4,
   -REPEAT_CLASS    => 'Simple_repeat');


$rca->store($rc);

ok($rc->dbID && $rc->adaptor());

$rc = $rca->fetch_by_dbID($rc->dbID);

ok($rc->repeat_consensus eq 'ACTG');
ok($rc->repeat_class eq  'Simple_repeat');
ok($rc->length() == 4);
ok($rc->name eq 'ACTG(n)');


$multi_db->restore('core', 'repeat_consensus');

