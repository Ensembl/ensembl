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
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Variation::StrainSlice;
use Bio::EnsEMBL::DensityFeature;
my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new('mus_musculus');
my $vdba = $multi_db->get_DBAdaptor('variation');
my $cdba = $multi_db->get_DBAdaptor('core');
$vdba->dnadb($cdba);

my $slice_adaptor = $cdba->get_SliceAdaptor;
my $slice = $slice_adaptor->fetch_by_region('chromosome', 19, 20380186, 20384187);

my $strain_name = 'A/J';
my $strain_slice = Bio::EnsEMBL::Variation::StrainSlice->new(
  -START   => $slice->{'start'},
  -END     => $slice->{'end'},
  -STRAND  => $slice->{'strand'},
  -ADAPTOR => $slice->adaptor(),
  -SEQ     => $slice->{'seq'},
  -SEQ_REGION_NAME => $slice->{'seq_region_name'},
  -SEQ_REGION_LENGTH => $slice->{'seq_region_length'},
  -COORD_SYSTEM    => $slice->{'coord_system'},
  -STRAIN_NAME     => $strain_name);

my $vf_adaptor = $vdba->get_VariationFeatureAdaptor;

my $features = $vf_adaptor->fetch_all_by_Slice($strain_slice);
ok(scalar @$features == 189, 'Count VariationFeatures on StrainSlice');

my $feature = $features->[0];
ok($feature && $feature->isa('Bio::EnsEMBL::Variation::VariationFeature'), 'Feature is a VariationFeature');

my $feature_slice = $feature->feature_Slice;
ok($feature_slice && $feature_slice->isa('Bio::EnsEMBL::Variation::StrainSlice'), 'Feature slice is a StrainSlice');

$multi_db->save('core', 'analysis');
$multi_db->save('core', 'density_type');
$multi_db->save('core', 'density_feature');
my $aa  = $cdba->get_AnalysisAdaptor();
my $analysis = new Bio::EnsEMBL::Analysis(-database => "ensembl", -logic_name => "SNPDensity");
ok(!$analysis->is_stored($cdba));
$aa->store($analysis);
ok($analysis->is_stored($cdba));

my $dta = $cdba->get_DensityTypeAdaptor();
my $dt = Bio::EnsEMBL::DensityType->new(
          -analysis   => $analysis,
          -block_size => 4000,
          -value_type => 'sum');

ok(!$dt->is_stored($cdba));
$dta->store($dt); 
ok($dt->is_stored($cdba));

my $df = Bio::EnsEMBL::DensityFeature->new(
          -seq_region    => $slice,
          -start         => $slice->{'start'},
          -end           => $slice->{'end'},
          -density_type  => $dt,
          -density_value => 5 );
my $dfa = $cdba->get_DensityFeatureAdaptor();
$dfa->store(($df));

$features = $dfa->fetch_all_by_Slice($strain_slice, 'SNPDensity', 10, 1);
is( @$features, 10, "Number of stored SNP densities");

$multi_db->restore('core', 'analysis');
$multi_db->restore('core', 'density_type');
$multi_db->restore('core', 'density_feature');
done_testing();
