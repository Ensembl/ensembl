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

#expect undef back
$range = $rr->check_and_register($id, 300, 500);
ok(!defined($range));
print_ranges($range);

#expect 700-900 back
$range = $rr->check_and_register($id, 700, 900);
ok($range && $range->[0]->[0] == 700 && $range->[0]->[1] == 900);
print_ranges($range);

# expect 1000-1200 back
$range = $rr->check_and_register($id, 1050, 1150, 1000, 1200);
ok($range && $range->[0]->[0] == 1000 && $range->[0]->[1] == 1200);
print_ranges($range);

#expect 50-99, 501-699, 901-950 back
$range = $rr->check_and_register($id, 50, 200, 50, 950);
ok($range && $range->[0]->[0] == 50 && $range->[0]->[1] == 99);
ok($range && $range->[1]->[0] == 501 && $range->[1]->[1] == 699);
ok($range && $range->[2]->[0] == 901 && $range->[2]->[1] == 950);
print_ranges($range);

#make sure that the interal list has been merged into 2 ranges
#we have to do this to make sure that it is efficient
my $internal_list = $rr->{'registry'}->{$id};
ok(@$internal_list == 2);

#check that creating adjacent regions merges the list correctly
$range = $rr->check_and_register($id, 40,45,30,49);
ok(@$internal_list == 2);
ok($range && $range->[0]->[0] == 40 && $range->[0]->[1] == 49);
print_ranges($range);


$range = $rr->check_and_register($id, 951, 999);
ok(@$internal_list == 1);
ok($range && $range->[0]->[0] == 951 && $range->[0]->[1] == 999);
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
