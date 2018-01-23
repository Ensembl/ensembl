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

use Bio::EnsEMBL::Utils::Exception qw(warning verbose throw info
                                      deprecate stack_trace_dump stack_trace);

our $verbose= 0;


if(!$verbose) {
  verbose('NONE');
}

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

$verbose && verbose('EXCEPTION');

ok(verbose() == 1000 || (!$verbose && verbose() == 0));

warning('This warn should not appear');

ok(1);

warning('This warn should appear', 1000);

ok(1);

info("This info should not appear");

ok(1);

$verbose && verbose('ALL');

info("This info should appear");

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

verbose('DEPRECATE');

done_testing();
