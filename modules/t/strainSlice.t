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
#use Test::Warnings;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new('mus_musculus');

my $vdba = $multi_db->get_DBAdaptor('variation');
my $cdba = $multi_db->get_DBAdaptor('core');
$vdba->dnadb($cdba);

my $slice_adaptor = $cdba->get_SliceAdaptor;
my $slice = $slice_adaptor->fetch_by_region('chromosome', 19, 20380186, 20384187);

my $variation_adaptor = $vdba->get_VariationAdaptor;
my $vf_adaptor = $vdba->get_VariationFeatureAdaptor;

my $vfs = $vf_adaptor->fetch_all_by_Slice($slice);
my $vf = $vfs->[0];

my $A_J = $slice->get_by_strain('A/J');
my $FVB_NJ = $slice->get_by_strain('FVB/NJ');
my $with_coverage = 1;

my $display_slice_name = $A_J->display_Slice_name;
ok($display_slice_name eq 'A/J', 'display_slice_name');

my $strain_name = $A_J->strain_name;
ok($strain_name eq 'A/J', 'strain_name');

my $sample = $A_J->sample;
ok($sample && $sample->isa('Bio::EnsEMBL::Variation::Sample'), 'isa Bio::EnsEMBL::Variation::Sample');
ok($sample->name eq 'A/J', 'sample name');

my $seq = $A_J->seq;
ok(length($seq) == 4002, 'seq length');

$seq = $A_J->seq($with_coverage);
ok(length($seq) == 4002, 'seq length');

my $expanded_length = $A_J->expanded_length;
ok($expanded_length == 4002, 'expanded length');

my $a_j_vfs = $A_J->get_all_VariationFeatures();
ok(scalar @$a_j_vfs == 46, 'get_all_VariationFeatures');

my $variation = $variation_adaptor->fetch_by_name('rs46873854');
$vfs = $vf_adaptor->fetch_all_by_Variation($variation);
$vf = $vfs->[0];
my $af = $A_J->get_AlleleFeature($vf);
ok($vf->allele_string eq 'C/T', 'VF allele_string');
ok($af->allele_string eq 'T|T', 'AF allele_string');
ok($af->ref_allele_string eq 'C', 'AF ref_allele_string');

my $afs = $A_J->get_all_AlleleFeatures_Slice();
ok(scalar @$afs == 46, 'get_all_AlleleFeatures_Slice');

$afs = $A_J->get_all_AlleleFeatures_Slice($with_coverage);
ok(scalar @$afs == 46, 'get_all_AllelelFeatures_Slice with_coverage');

$afs = $A_J->get_all_differences_StrainSlice($FVB_NJ);
ok(scalar @$afs == 53, 'get_all_differences_StrainSlice');

my $sub_Slice = $A_J->sub_Slice(2, 12, 1);
ok(length($sub_Slice->seq) == 11, 'sub_Slice length');

my $ref_subseq = $A_J->ref_subseq(20380187, 20380197, 1);
ok(length($ref_subseq) == 11, 'ref_subseq length');

my $subseq = $A_J->subseq(1, 10, 1);
ok($subseq eq 'CACTGTTCCC', 'subseq');

my $position = $A_J->get_original_seq_region_position(20);
ok($position == 20, 'get_original_seq_region_position');

done_testing();
