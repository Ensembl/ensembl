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

print "Testing UTRs on transcript fetched by stable_id\n";
try_utrs($transcript);

my $slice_adaptor = $db->get_SliceAdaptor();
my $slice = $slice_adaptor->fetch_by_location('20:31166507-31166507', 'chromosome');
my $features = $transcript_adaptor->fetch_all_by_Slice($slice);
$transcript = shift @{$features};

print "Testing UTRs on transcript fetched by slice\n";
try_utrs($transcript);

$stable_id = 'ENST00000246229';
$transcript = 
  $transcript_adaptor->fetch_by_stable_id($stable_id);

print "Testing UTRs on reverse strand transcript fetched by stable_id\n";
try_utrs_reverse($transcript);

my $slice = $slice_adaptor->fetch_by_location('20:30583588-30583588', 'chromosome');
my $features = $transcript_adaptor->fetch_all_by_Slice($slice);
my $transcript = shift @{$features};

print "Testing UTRs on reverse strand transcript fetched by slice\n";
try_utrs_reverse($transcript);

sub try_utrs {
    my $transcript = shift;

    my @five_utrs = @{ $transcript->get_all_five_prime_UTRs() };  
    my $five = scalar(@five_utrs);
    my @three_utrs = @{ $transcript->get_all_three_prime_UTRs() };
    my $three = scalar(@three_utrs);

    is($five_utrs[0]->start, $transcript->start, "Five prime starts at transcript start");
    is($five_utrs[$five-1]->end, $transcript->coding_region_start - 1, "Five prime ends starts the coding region");
    is($three_utrs[0]->start, $transcript->coding_region_end + 1, "Three prime starts at end of coding region");
    is($three_utrs[$three-1]->end, $transcript->end, "Three prime ends at transcript end");
    is($five_utrs[0]->seq_region_start, $transcript->seq_region_start, "Five prime seq_region_starts at transcript seq_region_start");
    is($three_utrs[$three-1]->seq_region_end, $transcript->seq_region_end, "Three prime seq_region_ends at transcript seq_region_end");


}

sub try_utrs_reverse {
    my $transcript = shift;

    my @five_utrs = @{ $transcript->get_all_five_prime_UTRs() };  
    my $five = scalar(@five_utrs);
    my @three_utrs = @{ $transcript->get_all_three_prime_UTRs() };
    my $three = scalar(@three_utrs);

    is($five_utrs[$five-1]->start, $transcript->coding_region_end + 1, "Five prime starts at transcript start");
    is($five_utrs[0]->end, $transcript->end, "Five prime ends starts the coding region");
    is($three_utrs[0]->start, $transcript->start, "Three prime starts at end of coding region");
    is($three_utrs[$three-1]->end, $transcript->coding_region_start - 1, "Three prime ends at transcript end");

}

done_testing();
