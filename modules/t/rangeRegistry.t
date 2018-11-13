# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;

use Test::More;
use Test::Warnings;

use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Mapper::RangeRegistry;

our $verbose= 0;

my $rr = Bio::EnsEMBL::Mapper::RangeRegistry->new();

my $id = 'ID1';

#expect [100,400] back
my $range = $rr->check_and_register($id, 200,300, 100,400);
ok(@$range==1 && $range->[0]->[0] == 100 && $range->[0]->[1] == 400);
print_ranges($range);

#expect undef back
$range = $rr->check_and_register($id, 150,190, 100,200);
ok(!defined($range));
print_ranges($range);

#expect [401,500] back
$range = $rr->check_and_register($id, 200, 500);
ok(@$range==1 && $range->[0]->[0] == 401 && $range->[0]->[1] == 500);
print_ranges($range);

#expect undef back
$range = $rr->check_and_register($id, 300, 500);
ok(!defined($range));
print_ranges($range);

#expect 700-900 back
$range = $rr->check_and_register($id, 700, 900);
ok(@$range==1 && $range->[0]->[0] == 700 && $range->[0]->[1] == 900);
print_ranges($range);

# expect 1000-1200 back
$range = $rr->check_and_register($id, 1050, 1150, 1000, 1200);
ok(@$range==1 && $range->[0]->[0] == 1000 && $range->[0]->[1] == 1200);
print_ranges($range);

#expect 50-99, 501-699, 901-950 back
$range = $rr->check_and_register($id, 50, 200, 50, 950);
ok(@$range==3 && $range->[0]->[0] == 50 && $range->[0]->[1] == 99);
ok(@$range==3 && $range->[1]->[0] == 501 && $range->[1]->[1] == 699);
ok(@$range==3 && $range->[2]->[0] == 901 && $range->[2]->[1] == 950);
print_ranges($range);

#make sure that the interal list has been merged into 2 ranges
#we have to do this to make sure that it is efficient
my $internal_list = $rr->{'registry'}->{$id};
ok(@$internal_list == 2);

#check that creating adjacent regions merges the list correctly
$range = $rr->check_and_register($id, 40,45,30,49);
ok(@$internal_list == 2);
ok(@$range==1 && $range->[0]->[0] == 30 && $range->[0]->[1] == 49);
print_ranges($range);

$range = $rr->check_and_register($id, 951, 999);
ok(@$internal_list == 1);
ok($range && $range->[0]->[0] == 951 && $range->[0]->[1] == 999);
print_ranges($range);


# Check that a single range can be added to the beginning
$range = $rr->check_and_register($id, 1, 10, 1,20);
ok(@$internal_list == 2);
ok(@$range==1 && $range->[0]->[0] == 1 && $range->[0]->[1] == 20);
print_ranges($range);

#check range that spans entire existing ranges
$range = $rr->check_and_register($id, 1, 1200);
ok(@$internal_list == 1);
ok(@$range==1 && $range->[0]->[0] == 21 && $range->[0]->[1] == 29);
print_ranges($range);

#check adding identical range to existing internal range
$range = $rr->check_and_register($id, 1, 1200);
ok(!defined($range));
print_ranges($range);


#check requesting small area of size 1
$range = $rr->check_and_register($id,10,10, 1, 1e6);
ok(!defined($range));
print_ranges($range);

#check that adding a range to a different id works
$range = $rr->check_and_register("ID2", 100,500, 1, 600);
ok($range && @$range==1 && $range->[0]->[0] == 1 && $range->[0]->[1] == 600);
print_ranges($range);


# I suspect there is a small bug in respect to handling the 
# second argument extended range

# setup some ranges
$range = $rr->check_and_register( "rr_bug", 2,10 );
$range = $rr->check_and_register( "rr_bug", 15,20 );
$range = $rr->check_and_register( "rr_bug", 25,30 );

my $overlap = $rr->overlap_size(  "rr_bug", 3, 40 );
ok( $overlap == 20 );

$range = $rr->check_and_register( "rr_bug", 28, 35, 3, 40 );
debug( "*** extended bug test ***" );
print_ranges( $range );
# this should result in 2,40 to be covered in the range registry
ok(@$range==3 && $range->[0]->[0] == 11 && $range->[0]->[1] == 14);
ok(@$range==3 && $range->[1]->[0] == 21 && $range->[1]->[1] == 24);
ok(@$range==3 && $range->[2]->[0] == 31 && $range->[2]->[1] == 40);



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

done_testing();
