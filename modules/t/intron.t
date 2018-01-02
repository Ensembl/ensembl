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

require_ok('Bio::EnsEMBL::Exon');
require_ok('Bio::EnsEMBL::Intron');

our $verbose = 0; #set to 1 to turn on debug printouts

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

ok(1);

my $db = $multi->get_DBAdaptor( 'core' );

ok($db);

my $stable_id = 'ENST00000217347';
my $transcript_adaptor = $db->get_TranscriptAdaptor();
my $transcript = 
  $transcript_adaptor->fetch_by_stable_id($stable_id);


my @exons = (@{$transcript->get_all_Exons()});  
my @introns = (@{$transcript->get_all_Introns()});  

my $rank=1;
foreach my $intron (@introns){
  is($intron->rank($transcript), $rank, "Checking intron rank");
  $rank++;
  ok($intron->prev_Exon->end == $intron->start-1);
  ok($intron->next_Exon->start == $intron->end+1);
  ok($intron->is_splice_canonical(), 'Checking Intron is canonical in its splicing');
}

#Make a fake exon pair
{
  my $slice = $transcript->slice();
  my $exon_one = Bio::EnsEMBL::Exon->new(-START => 10, -END => 20, -STRAND => 1, -SLICE => $slice);
  my $exon_two = Bio::EnsEMBL::Exon->new(-START => 30, -END => 40, -STRAND => 1, -SLICE => $slice);
  my $intron = Bio::EnsEMBL::Intron->new($exon_one, $exon_two);
  ok(!$intron->is_splice_canonical(), 'Checking fake Intron is not canonical');
}

$multi->restore();

done_testing();
