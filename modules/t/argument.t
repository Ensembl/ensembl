use strict;
use warnings;


BEGIN { $| = 1;  
	use Test;
	plan tests => 4;
}

use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);

our $verbose= 0;

my @args = ('-TWO' => 2,
            '-ONE' => 1,
            '-THREE' => 3);

my ($one, $two, $three) = rearrange(['ONE','TWO','THREE'],@args);

ok($one == 1);
ok($two == 2);
ok($three == 3);

#
#regression test args with 0 were being set to undef instead of 0
#
@args = ('-ZERO' => 0, '-ONE' => 1);
my $zero;
($one, $zero) = rearrange(['ONE', 'ZERO'], @args);
ok(defined($zero) && $zero == 0 && $one == 1);
