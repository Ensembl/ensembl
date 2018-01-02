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

my @ets = @{ $transcript->get_all_ExonTranscripts() };
my @exons = @{ $transcript->get_all_Exons() };

my $n = scalar(@exons);

for (my $i = 0; $i < $n; $i++) {
  is($ets[$i]->start, $exons[$i]->start, "Starts match");
  is($ets[$i]->end, $exons[$i]->end, "Ends match");
  is($ets[$i]->rank, $i+1, "Rank is correct");
}


done_testing();
