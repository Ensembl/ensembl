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
use Bio::EnsEMBL::MappedSliceContainer;
use Bio::EnsEMBL::Variation::DBSQL::StrainSliceAdaptor;
my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new('mus_musculus');
my $vdba = $multi_db->get_DBAdaptor('variation');
my $cdba = $multi_db->get_DBAdaptor('core');
$vdba->dnadb($cdba);

my $slice_adaptor = $cdba->get_SliceAdaptor;
my $slice = $slice_adaptor->fetch_by_region('chromosome', 19, 20380186, 20384187);

my $msc = Bio::EnsEMBL::MappedSliceContainer->new(-SLICE => $slice, -EXPANDED => 1, -ADAPTOR => $slice_adaptor);

$msc->set_StrainSliceAdaptor(Bio::EnsEMBL::Variation::DBSQL::StrainSliceAdaptor->new($vdba));

my $strain_name = 'A/J';
$msc->attach_StrainSlice($strain_name);

my @mapped_slices = @{$msc->get_all_MappedSlices};
ok(scalar(@mapped_slices) == 1, 'Return 1 MappedSlice'); 

my $mapped_slice = $mapped_slices[0];
my $container = $mapped_slice->container;
ok($container && $container->isa('Bio::EnsEMBL::MappedSliceContainer'), 'isa Bio::EnsEMBL::MappedSliceContainer');

my $seq = $mapped_slice->seq(1);
ok(length($seq) == 4002, 'sequence length');

my $coord_system_name = $mapped_slice->coord_system_name;
ok($coord_system_name eq 'container', 'coord_system_name');

my $coord_system = $mapped_slice->coord_system;
ok($coord_system_name eq 'container', 'coord_system');

my $adaptor = $mapped_slice->adaptor;
ok($adaptor && $adaptor->isa('Bio::EnsEMBL::Variation::DBSQL::StrainSliceAdaptor'), 'isa Bio::EnsEMBL::Variation::DBSQL::StrainSliceAdaptor');

my $slice_mapper_pairs = $mapped_slice->get_all_Slice_Mapper_pairs;
ok(scalar @{$slice_mapper_pairs} == 1, 'count slice mapper pairs');

done_testing();
