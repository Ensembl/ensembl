use strict;
use warnings;

use lib 't';

BEGIN { $| = 1;  
	use Test;
	plan tests => 7;
}

use TestUtils qw( debug );

use Bio::EnsEMBL::Utils::Exception qw(warning verbose throw 
                                      deprecate stack_trace_dump stack_trace);

our $verbose= 0;


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

verbose(-1);

ok(verbose() == -1);

warning('This warn should not appear');

ok(1);

warning('This warn should appear');

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


verbose(0);