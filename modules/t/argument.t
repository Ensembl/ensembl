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

done_testing();
