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
use Test::More;
use Test::Warnings;

use Bio::EnsEMBL::Map::MapLocation;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0; #set to 1 to turn on debug printouts


my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $multi->get_DBAdaptor( 'core' );


######
# 1  #
######

#test constructor
my $mapname = 'genethon';
my $name = 'DS1234';
my $chr_name = 'X';
my $pos = '12.5';
my $lod = 0.23;

my $mloc = 
  Bio::EnsEMBL::Map::MapLocation->new($name, $mapname, $chr_name, $pos, $lod);

ok($mloc && ref $mloc && $mloc->isa('Bio::EnsEMBL::Map::MapLocation'));



#######
# 2-3 #
#######

#test map_name

ok($mapname eq $mloc->map_name);
ok(&test_getter_setter($mloc, 'map_name', 'marshfield'));

#######
# 4-5 #
#######

#test chromosome_name
ok($chr_name eq $mloc->chromosome_name);
ok(&test_getter_setter($mloc, 'chromosome_name', undef));

#######
# 6-7 #
#######

#test name
ok($name eq $mloc->name);
ok(&test_getter_setter($mloc, 'name', 'Z123213'));


#######
# 8-9 #
#######

#test position
ok($pos eq $mloc->position);
ok(&test_getter_setter($mloc, 'position', 'q13.1'));

########
# 10-11#
########

#test lod_score

ok($lod == $mloc->lod_score);
ok(&test_getter_setter($mloc, 'lod_score', 0.03 ));

done_testing();
