use strict;
use warnings;
use Test::More;
use Test::Exception;
use Bio::EnsEMBL::Test::MultiTestDB;
use DBI qw/:sql_types/;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $multi->get_DBAdaptor('core');
my $gene_adaptor = $dba->get_GeneAdaptor();

$gene_adaptor->bind_param_generic_fetch('protein_coding', SQL_VARCHAR);
my $count = $gene_adaptor->generic_count('g.biotype =?');
is($count, 20, 'Checking generic_count for protein_coding genes returns expected amounts');


# fetch_all_by_dbID_list tests

my $gene_list = $gene_adaptor->_uncached_fetch_all_by_id_list([qw(ENSG00000101321 ENSG00000101346 ENSG00000101367)],undef,"stable_id");
ok(scalar(@$gene_list) == 3, "Basic uncached fetch by list");

dies_ok {
    $gene_list = $gene_adaptor->_uncached_fetch_all_by_id_list([qw(ENSG00000101321 ENSG00000101346 ENSG00000101367)],undef,"dbID")
} "Wrong data type for ID";

done_testing();