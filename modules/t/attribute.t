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

use Bio::EnsEMBL::Attribute;

use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0; #set to 1 to turn on debug printouts

use Test::More;
use Test::Warnings;

use Bio::EnsEMBL::Test::TestUtils;

#
# test constructor
#

my $code = 'testcode';
my $name = 'testname';
my $desc = 'testdesc';
my $value = 'testval';

my $attrib = Bio::EnsEMBL::Attribute->new
  (-CODE => $code,
   -NAME => $name,
   -DESCRIPTION => $desc,
   -VALUE => $value);

ok($attrib->code()  eq $code);
ok($attrib->name()  eq $name);
ok($attrib->description()  eq $desc);
ok($attrib->value() eq $value);

#
# test getter/setters
#
ok(test_getter_setter($attrib, 'name', 'newname'));
ok(test_getter_setter($attrib, 'code', 'newcode'));
ok(test_getter_setter($attrib, 'description', 'newdesc'));
ok(test_getter_setter($attrib, 'value', 'newvalue'));

done_testing();
