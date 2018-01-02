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

use Bio::EnsEMBL::DensityFeature;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Analysis;


our $verbose = 0; #set to 1 to turn on debug printouts
use Test::More;
use Test::Warnings;
use Bio::EnsEMBL::Test::TestUtils;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db = $multi->get_DBAdaptor('core');


my $analysis =  new Bio::EnsEMBL::Analysis (-program     => "densityFeature.t",
					   -database    => "ensembl",
					   -gff_source  => "densityFeature.t",
					   -gff_feature => "density",
					   -logic_name  => "GeneDensityTest");

my $dt = Bio::EnsEMBL::DensityType->new(-analysis   => $analysis,
					-block_size => 600,
					-value_type => 'sum');

my $slice_adaptor = $db->get_SliceAdaptor();
my $slice = $slice_adaptor->fetch_by_region('chromosome', '20', 1, 600);


#
#test the constructor
#
my $feat = Bio::EnsEMBL::DensityFeature->new(-seq_region    => $slice,
				             -start         => 1,
					     -end           => 300,
					     -density_type  => $dt,
					     -density_value => 123);

ok($feat && ref $feat && $feat->isa('Bio::EnsEMBL::DensityFeature'));


#
# Test the getter setter functions;
#

ok(&test_getter_setter($feat, 'start', 100));
ok(&test_getter_setter($feat, 'end', 500));
ok(&test_getter_setter($feat, 'density_value', 456));

done_testing();
