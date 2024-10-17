# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
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

use Bio::EnsEMBL::CoordSystem;

use Bio::EnsEMBL::Test::TestUtils;

use Test::Exception;

our $verbose = 0;

use Test::More;
use Test::Warnings;

my $name    = 'chromosome';
my $version = 'NCBI33';
my $dbID    = 1;
my $top_level = 0;
my $sequence_level = 0;
my $default = 1;
my $rank    = 1;

#
# Test constructor
#
my $coord_system = Bio::EnsEMBL::CoordSystem->new
  (-NAME    => $name,
   -VERSION => $version,
   -DBID    => $dbID,
   -TOP_LEVEL => $top_level,
   -RANK    => $rank,
   -SEQUENCE_LEVEL => $sequence_level,
   -DEFAULT => 1);


ok($coord_system && $coord_system->isa('Bio::EnsEMBL::CoordSystem'));


#
# Test name()
#
ok($coord_system->name() eq $name);

#
# Test version()
#
ok($coord_system->version() eq $version);

#
# Test is_top_level()
#
ok(!$coord_system->is_top_level());

#
# Test is_sequence_level()
#
ok(!$coord_system->is_sequence_level());

#
# Test is_default()
#
ok($coord_system->is_default());

#
# Test rank()
#
ok($coord_system->rank() == $rank);

#
# Test equals()
#

my $coord_system2 = Bio::EnsEMBL::CoordSystem->new
  (-NAME    => $name,
   -VERSION => $version,
   -DBID    => 123,
   -RANK    => $rank,
   -TOP_LEVEL => $top_level);

ok($coord_system->equals($coord_system2));

$coord_system2 = Bio::EnsEMBL::CoordSystem->new
  (-NAME    => 'chromosome',
   -VERSION => 'NCBI34',
   -DBID    => 123,
   -RANK    => $rank,
   -TOP_LEVEL => $top_level);

ok(!$coord_system->equals($coord_system2));

#
# test creation of toplevel system
#
$name    = 'toplevel';
$top_level = 1;

#
# Test constructor
#
$coord_system = Bio::EnsEMBL::CoordSystem->new
  (-NAME    => $name,
   -TOP_LEVEL => $top_level);

ok($coord_system->name() eq $name);
ok($coord_system->is_top_level());
ok($coord_system->rank() == 0);

#
# Test constructor with alias set
#
{
  my $name = "primary_assembly";
  my $alias = "chromosome";
  my $coord_system = Bio::EnsEMBL::CoordSystem->new
  (-NAME     => $name,
   -RANK     => $rank,
   -ALIAS_TO => $alias);

  is($coord_system->alias_to(), $alias, "Correctly initialised alias variable in constructor to '$alias'");

  my $bad_alias = "badalias";
  throws_ok {
    my $cs_bad_alias = Bio::EnsEMBL::CoordSystem->new
    (-NAME     => $name,
    -RANK     => $rank,
    -ALIAS_TO => $bad_alias);
  } qr /The ALIAS_TO argument can only be defined as 'chromosome'/,
  "Checks that error is thrown if -ALIAS_TO in constructor is defined as something other than '$alias'";
}


#
# Test alias_to getter/setter
#
my $alias = 'chromosome';
ok(test_getter_setter($coord_system, 'alias_to', $alias)); 
is($coord_system->alias_to(), $alias, "Getter method correctly retrieved alias '$alias'");

throws_ok{ $coord_system->alias_to('somethingelse') } qr/The alias can only be set to/,
  'Checks that error is thrown if requested alias is not called chromosome';

done_testing();
