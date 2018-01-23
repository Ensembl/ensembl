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

use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Test::MultiTestDB;

use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0; #set to 1 to turn on debug printouts
use Test::More;
use Test::Warnings;
use Bio::EnsEMBL::Test::TestUtils;



my $analysis = Bio::EnsEMBL::Analysis->new
  (-program     => "test",
   -database    => "ensembl",
   -gff_source  => "densityFeature.t",
   -gff_feature => "density",
   -logic_name  => "GeneDensityTest");


#
# test constructor
#
my $dt = Bio::EnsEMBL::DensityType->new
  (-dbID       => 1200,
   -analysis   => $analysis,
   -block_size => 600,
   -value_type => 'sum');

ok($dt && ref($dt) && $dt->isa('Bio::EnsEMBL::DensityType'));
ok($dt->dbID == 1200);
ok($dt->analysis == $analysis);
ok($dt->block_size == 600);
ok($dt->value_type eq 'sum');


#
# test getter/setter methods
#
ok(&test_getter_setter($dt, 'dbID', 12));
ok(&test_getter_setter($dt, 'analysis', undef));
ok(&test_getter_setter($dt, 'block_size', 300));
ok(&test_getter_setter($dt, 'value_type', 'ratio'));

done_testing();
