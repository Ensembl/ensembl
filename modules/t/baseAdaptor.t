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
use Test::Exception;
use Bio::EnsEMBL::Test::MultiTestDB;
use DBI qw/:sql_types/;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $multi->get_DBAdaptor('core');
my $gene_adaptor = $dba->get_GeneAdaptor();

$gene_adaptor->bind_param_generic_fetch('protein_coding', SQL_VARCHAR);
my $count = $gene_adaptor->generic_count('g.biotype =?');
is($count, 21, 'Checking generic_count for protein_coding genes returns expected amounts');


# fetch_all_by_dbID_list tests
{
  my $gene_list = $gene_adaptor->_uncached_fetch_all_by_id_list([qw(ENSG00000101321 ENSG00000101346 ENSG00000101367)],undef,"stable_id");
  ok(scalar(@$gene_list) == 3, "Basic uncached fetch by list");
}
{
  my $gene_list = $gene_adaptor->_uncached_fetch_all_by_id_list([qw(ENSG00000101321 ENSG00000101346 ENSG00000101367)], undef, "stable_id", 0);
  ok(scalar(@$gene_list) == 3, "Basic uncached fetch by list forcing numeric type to false");
}

# CHecking non-numeric given to a numeric means no-go
dies_ok {
    $gene_adaptor->_uncached_fetch_all_by_id_list([qw(ENSG00000101321 ENSG00000101346 ENSG00000101367)], undef, "dbID")
} "Wrong data type for dbID";

# We want the proper exception telling us that we said it was a numeric but then just
# gave it the wrong data
throws_ok {
    $gene_adaptor->_uncached_fetch_all_by_id_list([qw(ENSG00000101321 ENSG00000101346 ENSG00000101367)], undef, "stable_id", 1)
} qr/You specified/, "Correct type given but forcing it to use numerics";

# Making sure we can execute an excessivly large list of IDs by splitting into multiple queries
my @ids = 1..256000;
my $id_derrived_gene_list = $gene_adaptor->_uncached_fetch_all_by_id_list(\@ids, undef, "dbID", 1);
is(scalar(@{$id_derrived_gene_list}), 21, 'Checking we get 21 genes back');

done_testing();
