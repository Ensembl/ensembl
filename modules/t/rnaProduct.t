# Copyright [2018] EMBL-European Bioinformatics Institute
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

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::RNAProduct;

use Test::More;
use Test::Warnings;

my $loaded = 0;
END { print "not ok 1\n" unless $loaded; }

#turn on/off debug prints:
our $verbose = 0;

use Bio::EnsEMBL::Test::MultiTestDB;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

$loaded = 1;

ok(1);

my $db = $multi->get_DBAdaptor('core');

my $rp = Bio::EnsEMBL::RNAProduct->new();

ok($rp, 'RNAProduct constructor works without arguments');


# TODO: More RNAProduct tests


#
# Tests for the mature-RNA adaptor
##################################

my $rp_a  = $db->get_RNAProductAdaptor();

ok($rp_a, 'Can get RNAProductAdaptor from core DBAdaptor');


# TODO: More RNAProductAdaptor tests


# Test generic_count(), inherited method from BaseAdaptor
is($rp_a->generic_count(), @{$rp_a->list_dbIDs()}, "Number of features from generic_count is equal to the number of dbIDs from list_dbIDs");

done_testing();
