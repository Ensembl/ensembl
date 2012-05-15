use strict;
use warnings;
use Test::More;
use Bio::EnsEMBL::Test::MultiTestDB;
use DBI qw/:sql_types/;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $multi->get_DBAdaptor('core');
my $gene_adaptor = $dba->get_GeneAdaptor();

$gene_adaptor->bind_param_generic_fetch('protein_coding', SQL_VARCHAR);
my $count = $gene_adaptor->generic_count('g.biotype =?');
is($count, 20, 'Checking generic_count for protein_coding genes returns expected amounts');

done_testing();