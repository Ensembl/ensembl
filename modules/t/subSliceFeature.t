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
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::SubSlicedFeature;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $multi->get_DBAdaptor('core');

my $ga = $dba->get_GeneAdaptor();
my $gene = $ga->fetch_by_stable_id('ENSG00000125964');

note($gene->stable_id);
note($gene->start);
note($gene->end);

my $transcript_list = $gene->get_all_Transcripts;

note(scalar(@$transcript_list));
foreach (@$transcript_list) {
    note($_->stable_id);
    note($_->start);
    note($_->end);
}
my $fake_gene = Bio::EnsEMBL::SubSlicedFeature->new(-feature => $gene, -start => 30840810, -end => 30859270);

$transcript_list = $fake_gene->get_all_Transcripts;
note (scalar(@$transcript_list));
foreach (@$transcript_list) {
    note($_->stable_id);
    note($_->start);
    note($_->end);
}

is ($transcript_list->[0]->stable_id,"ENST00000216932", "Only one transcript found in subsliced Gene");

my $exon_list = $fake_gene->get_all_Exons;
foreach (@$exon_list) {
    note($_->stable_id);
    note($_->start);
    note($_->end);
}
ok(scalar(@$exon_list) == 4, "Correct Exons for subsliced Gene");
is($exon_list->[0]->stable_id,'ENSE00001048819', "Correct Exon returned for subsliced Gene");

$fake_gene = Bio::EnsEMBL::SubSlicedFeature->new(-feature => $gene, -start => 1, -end => 2);
$transcript_list = $fake_gene->get_all_Transcripts;
ok (scalar(@$transcript_list) == 0, "Out of bounds search for features");

# Check normal feature functions are unaffected.

ok($fake_gene->start == 30822726, "Normal feature function behaves through proxy");

done_testing;
