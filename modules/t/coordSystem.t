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

use Bio::EnsEMBL::CoordSystem;

use Bio::EnsEMBL::Test::TestUtils;

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

done_testing();
