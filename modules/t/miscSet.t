use lib 't';
use strict;

BEGIN { $| = 1;
	use Test ;
	plan tests => 10
}

use Bio::EnsEMBL::MiscSet;
use TestUtils qw(debug test_getter_setter);

our $verbose = 0; #set to 1 to turn on debug printouts



my $misc_set_id = 3;
my $code = '1mbcloneset';
my $name = '1mb Clone Set';
my $desc = 'This is a 1MB cloneset';
my $max_len = 1e7;


#
# Test constructor
#
my $ms = Bio::EnsEMBL::MiscSet->new($misc_set_id,
                                    undef,
                                   $code,
                                   $name,
                                   $desc,
                                   $max_len);

ok($ms->dbID == $misc_set_id);
ok($ms->code eq $code);
ok($ms->description eq $desc);
ok($ms->longest_feature eq $max_len);
ok($ms->name eq $name);

ok(test_getter_setter($ms, 'dbID', 12));
ok(test_getter_setter($ms, 'code', 'testcode'));
ok(test_getter_setter($ms, 'description', 'new description'));
ok(test_getter_setter($ms, 'longest_feature', 1e8));
ok(test_getter_setter($ms, 'name', $name));


