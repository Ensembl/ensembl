use strict;
use warnings;

use lib 't';

BEGIN { $| = 1;  
	use Test;
	plan tests => 3;
}

use TestUtils qw( debug );

use Bio::EnsEMBL::Utils::Argument qw(rearrange);

our $verbose= 0;

my @args = ('-TWO' => 2,
            '-ONE' => 1,
            '-THREE' => 3);

my ($one, $two, $three) = rearrange(['ONE','TWO','THREE'],@args);

ok($one == 1);
ok($two == 2);
ok($three == 3);
