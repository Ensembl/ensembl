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

use Config;
use Test::More;
use Test::Exception;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::ApiVersion;

my $version = software_version()-1;
Bio::EnsEMBL::Registry->load_registry_from_url('mysql://anonymous@ensembldb.ensembl.org:5306/'.$version);

my $dbas = Bio::EnsEMBL::Registry->get_all_DBAdaptors(-GROUP=>'core');
my $min = 10;
ok(defined $dbas && scalar(@$dbas)>$min,'More than '.$min.' DBAs found');
ok($dbas->[0]->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'),'First DBA is core');

my $species = 'homo_sapiens';
my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
ok(defined $dba && $dba->isa('Bio::EnsEMBL::DBSQL::DBAdaptor') && ($dba->species() eq $species), "$species adaptor retrieved");

done_testing();

