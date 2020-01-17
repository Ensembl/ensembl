# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

use Bio::EnsEMBL::Test::MultiTestDB;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $multi_config = $multi->db_conf();
if (lc($multi_config->{'driver'}) ne 'mysql') {
  plan skip_all => 'Registry only supports MySQL for now';
}

my $mysql_connect_string = "$multi_config->{'driver'}://";
if (exists $multi_config->{'user'}) {
  $mysql_connect_string .= "$multi_config->{'user'}";
  if (exists $multi_config->{'pass'}) {
    $mysql_connect_string .= ":$multi_config->{'pass'}";
  }
  $mysql_connect_string .= '@';
}
$mysql_connect_string .= $multi_config->{'host'};
if (exists $multi_config->{'port'}) {
  $mysql_connect_string .= ":$multi_config->{'port'}";
}
Bio::EnsEMBL::Registry->load_registry_from_url($mysql_connect_string);

my $dbas = Bio::EnsEMBL::Registry->get_all_DBAdaptors(-GROUP=>'core');
isnt($dbas, undef, 'Can retrieve list of DBAs from Registry');
cmp_ok(scalar @{$dbas}, '>', 0, 'DBA list is not empty');
ok($dbas->[0]->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'),'First DBA is core');

my $species = 'homo_sapiens';
my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
ok(defined $dba && $dba->isa('Bio::EnsEMBL::DBSQL::DBAdaptor') && ($dba->species() eq $species), "$species adaptor retrieved");

done_testing();

