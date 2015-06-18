# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

require_ok('Bio::EnsEMBL::UTR');

our $verbose = 0; #set to 1 to turn on debug printouts

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

ok(1);

my $db = $multi->get_DBAdaptor( 'core' );

ok($db);

my $stable_id = 'ENST00000217347';
my $transcript_adaptor = $db->get_TranscriptAdaptor();
my $transcript = 
  $transcript_adaptor->fetch_by_stable_id($stable_id);


my @five_utrs = @{ $transcript->get_all_five_prime_utrs() };  
my $five = scalar(@five_utrs);
my @three_utrs = @{ $transcript->get_all_three_prime_utrs() };
my $three = scalar(@three_utrs);

is($five_utrs[0]->start, $transcript->start, "Five prime starts at transcript start");
is($five_utrs[$five-1]->end + 1, $transcript->coding_region_start, "Five prime ends starts the coding region");
is($three_utrs[0]->start - 1, $transcript->coding_region_end, "Three prime starts at end of coding region");
is($three_utrs[$three-1]->end, $transcript->end, "Three prime ends at transcript end");

done_testing();
