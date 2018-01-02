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
use Bio::EnsEMBL::Slice;

our $verbose = 0;

#
# TEST - Slice Compiles
#
ok(1);


my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');
my $vdb = $multi_db->get_DBAdaptor('variation');
#
# TEST - Slice creation from adaptor
#
my $slice_adaptor = $db->get_SliceAdaptor;
# tests for variation related methods
my $slice = $slice_adaptor->fetch_by_region('chromosome', 20, 30_252_000, 31_252_001);
my $study_adaptor = $vdb->get_StudyAdaptor;
my $variation_set_adaptor = $vdb->get_VariationSetAdaptor;
my $population_adaptor = $vdb->get_PopulationAdaptor;

my $vfs = $slice->get_all_VariationFeatures;
is(scalar @$vfs, 2, 'get_all_VariationFeatures');

$vfs = $slice->get_all_somatic_VariationFeatures;
is(scalar @$vfs, 1, 'get_all_somatic_VariationFeatures');

$vfs = $slice->get_all_somatic_VariationFeatures_by_source('dbSNP');
is(scalar @$vfs, 1, 'get_all_somatic_VariationFeatures_by_source');

$vfs = $slice->get_all_VariationFeatures_with_phenotype;
is(scalar @$vfs, 1, 'get_all_VariationFeatures_with_phenotype');

$vfs = $slice->get_all_somatic_VariationFeatures_with_phenotype;
is(scalar @$vfs, 1, 'get_all_somatic_VariationFeatures_with_phenotype');

my $svfs = $slice->get_all_StructuralVariationFeatures;
is(scalar @$svfs, 1, 'get_all_StructuralVariationFeatures');

$svfs = $slice->get_all_somatic_StructuralVariationFeatures;
is(scalar @$svfs, 1, 'get_all_somatic_StructuralVariationFeatures');

my $size_min         = 1;
my $size_max         = 10000;
my $include_evidence = 0;
$svfs = $slice->get_all_StructuralVariationFeatures_by_size_range($size_min, $size_max, $include_evidence);
is(scalar @$svfs, 1, 'get_all_StructuralVariationFeatures_by_size_range');

$svfs = $slice->get_all_somatic_StructuralVariationFeatures_by_size_range($size_min, $size_max, $include_evidence);
is(scalar @$svfs, 1, 'get_all_somatic_StructuralVariationFeatures_by_size_range');

my $study = $study_adaptor->fetch_by_name('study_name');
$svfs = $slice->get_all_StructuralVariationFeatures_by_Study($study);
is(scalar @$svfs, 3, 'get_all_StructuralVariationFeatures_by_Study');

$svfs = $slice->get_all_StructuralVariationFeatures_by_source('dbSNP');
is(scalar @$svfs, 1, 'get_all_StructuralVariationFeatures_by_source');

$svfs = $slice->get_all_somatic_StructuralVariationFeatures_by_source('dbSNP');
is(scalar @$svfs, 1, 'get_all_somatic_StructuralVariationFeatures_by_source');

my $variation_set = $variation_set_adaptor->fetch_by_name('variation_set');
$vfs = $slice->get_all_VariationFeatures_by_VariationSet($variation_set);
is(scalar @$vfs, 2, 'get_all_VariationFeatures_by_VariationSet');

$svfs = $slice->get_all_StructuralVariationFeatures_by_VariationSet($variation_set);
is(scalar @$svfs, 2, 'get_all_StructuralVariationFeatures_by_VariationSet');

my $pfs = $slice->get_all_PhenotypeFeatures;
is(scalar @$pfs, 2, 'get_all_PhenotypeFeatures');

$svfs = $slice->get_all_CopyNumberVariantProbeFeatures;
is(scalar @$svfs, 1, 'get_all_CopyNumberVariantProbeFeatures');

$vfs = $slice->get_all_genotyped_VariationFeatures;
is(scalar @$vfs, 1, 'get_all_genotyped_VariationFeatures');

my $population = $population_adaptor->fetch_by_name('population');
$vfs = $slice->get_all_VariationFeatures_by_Population($population);
is(scalar @$vfs, 1, 'get_all_VariationFeatures_by_Population');
done_testing();
