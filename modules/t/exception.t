use strict;
use warnings;

use Test::More;
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Utils::Exception qw(warning verbose throw info
                                      deprecate stack_trace_dump stack_trace);

our $verbose= 0;


if(!$verbose) {
  verbose('NONE');
}

#
#1 - test throw
#

eval {
  throw('test exception');
};

ok($@);

debug($@);

#
#2-4 - Test verbosity, warnings
#

$verbose && verbose('EXCEPTION');

ok(verbose() == 1000 || (!$verbose && verbose() == 0));

warning('This warn should not appear');

ok(1);

warning('This warn should appear', 1000);

ok(1);

info("This info should not appear");

ok(1);

$verbose && verbose('ALL');

info("This info should appear");

ok(1);

#
# 5-6 Test stack trace
#

my $std = stack_trace_dump();
ok($std =~ /[A-Z]/);

debug(stack_trace_dump);

ok(stack_trace());


#
# 7 Test deprecate
#

test_deprecate();

ok(7);

sub test_deprecate {
  deprecate('This deprecate warning should appear');
}

verbose('DEPRECATE');

done_testing();