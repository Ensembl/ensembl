use strict;
use warnings;

use lib 't';

BEGIN { $| = 1;
	use Test;
	plan tests => 40;
}

use TestUtils qw( debug );

use Bio::EnsEMBL::Mapper::RangeRegistry;

our $verbose= 1;

my $rr = Bio::EnsEMBL::Mapper::RangeRegistry->new();

my $id = 'ID1';

#expect [100,400] back
my $range = $rr->check_and_register($id, 200,300, 100,400);
ok($range && $range->[0]->[0] == 100 && $range->[0]->[1] == 400);
print_ranges($range);

#expect undef back
$range = $rr->check_and_register($id, 150,190, 100,200);
print_ranges($range);
ok(!defined($range));

#expect [401,500] back
$range = $rr->check_and_register($id, 200, 500);
ok($range && $range->[0]->[0] == 401 && $range->[0]->[1] == 500);
print_ranges($range);


sub print_ranges {
  my $rangelist = shift;
  if(!defined($rangelist)) {
    debug("UNDEF");
    return;
  }

  foreach my $range (@$rangelist) {
    debug('['.$range->[0].'-'.$range->[1].']');
  }
}
