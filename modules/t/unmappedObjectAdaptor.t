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
use Bio::EnsEMBL::UnmappedObject;
use Bio::EnsEMBL::DBSQL::UnmappedObjectAdaptor;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
ok( $multi );

my $db = $multi->get_DBAdaptor( "core" );

my $uma = $db->get_UnmappedObjectAdaptor();

ok(ref($uma));

#
# test fetch operations
#

my @objects = @{$uma->fetch_all()};

ok(scalar(@objects) == 4);

@objects = @{$uma->fetch_all_by_type('xref')};

ok(scalar(@objects) == 2);

my $analysis_adaptor = $db->get_AnalysisAdaptor();

my $analysis = $analysis_adaptor->fetch_by_logic_name("Unigene");

if(!defined($analysis)){
   die "ARSE\n";
}		
@objects = @{$uma->fetch_all_by_analysis($analysis)};

ok(scalar(@objects) == 2);

@objects = @{$uma->fetch_all_by_analysis($analysis,"UniGene")};

ok(scalar(@objects) == 2);

@objects = @{$uma->fetch_all_by_analysis($analysis,"RFAM")};

ok(scalar(@objects) == 0);

@objects = @{$uma->fetch_by_identifier("X5678")};

ok(scalar(@objects) == 1);

@objects = @{$uma->fetch_by_identifier("X5678","UniGene")};

ok(scalar(@objects) == 1);

@objects = @{$uma->fetch_by_identifier("X5678","RFAM")};

ok(scalar(@objects) == 0);

done_testing();
