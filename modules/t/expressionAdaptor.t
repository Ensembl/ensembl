use strict;

use Bio::EnsEMBL::Test::TestUtils;

use Test::More;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Expression;

our $verbose = 1;
our $clean = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;

# get a core DBAdaptor
my $db = $multi->get_DBAdaptor("core");

my $sa = $db->get_SliceAdaptor();
my $ea = $db->get_ExonAdaptor();
my $ga = $db->get_GeneAdaptor();


#
# Test get_ExpressionAdaptor works
#
my $ea = $db->get_ExpressionAdaptor();

ok($ea && ref($ea) && $ea->isa('Bio::EnsEMBL::DBSQL::ExpressionAdaptor'));


# hide the contents of the attrib_type, misc_attrib, seq_region_attrib tables
# so we can test storing etc. with a clean slate
$multi->hide('core', 'gene_expression', 'tissue');


#################
# Gene functionality tests
#


my $tissue = Bio::EnsEMBL::Expression->new(-NAME => 'test_name2',
                                           -ONTOLOGY => 'test_ontology2',
                                           -DESCRIPTION => 'test_desc2',
                                           -VALUE => 'test_value2');


my $gene = $ga->fetch_by_stable_id('ENSG00000131044');

$ea->store_on_Gene($gene, [$tissue]);

#
# make sure the gene_expression table was updated
#
my $count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT count(*) FROM gene_expression " .
   "WHERE gene_id = ".$gene->dbID())->[0]->[0];

is($count, 1, "Stored a gene expression for ENSG00000131044");

#
# make sure the tissue table was updated
#
$count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT count(*) FROM tissue " .
   "WHERE name = 'test_name2'")->[0]->[0];
is($count, 1, "Stored a tissue test_name2");

#
# test that we can now retrieve this expression
#
my @expressions = @{$ea->fetch_all_by_Gene($gene)};
is(@expressions, 1, "Fetched one expression for ENSG00000131044");

@expressions = @{$ea->fetch_all_by_Gene($gene,"rubbish")};
is(@expressions, 0, "Fetched no expression for ENSG00000131044 and tissue rubbish");

@expressions = @{$ea->fetch_all_by_Gene($gene,"test_name2")};
is(@expressions, 1, "Fetched one expression for ENSG00000131044 and tissue test_name2");

@expressions = @{$ea->fetch_all_by_Gene(undef,"test_name2")};
is(@expressions, 1, "Fetched one expression for tissue test_name2");


$tissue = $expressions[0];

is($tissue->name, 'test_name2', "Tissue name is test_name2");
is($tissue->ontology, 'test_ontology2', "Tissue ontology is test_ontology2");
is($tissue->description, 'test_desc2', "Tissue description is test_desc2");
is($tissue->value, 'test_value2', "Expression value is test_value2");


#
# test the removal of this expression with tissue code
#
note("Removing junk from gene");
$ea->remove_from_Gene($gene,"junk");
$count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT count(*) FROM gene_expression " .
   "WHERE gene_id = " . $gene->dbID())->[0]->[0];

is($count, 1, "One entry in gene_expression table for gene ENSG00000131044");




#
# test the removal of this expression
#

note("Removing test_name2 from gene");
$ea->remove_from_Gene($gene,"test_name2");
$count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT count(*) FROM gene_expression " .
   "WHERE gene_id = " . $gene->dbID())->[0]->[0];

is($count, 0, "No entry in gene_expression table for gene ENSG00000131044 and tissue test_name2");



#
# make sure the expression is no longer retrievable
#
@expressions = @{$ea->fetch_all_by_Gene($gene)};
is(@expressions, 0, "No expressions any more for ENSG00000131044");



#
# try to add an expression with an already existing code
#
$ea->store_on_Gene($gene, [$tissue]);
#
# make sure the gene_expression table was updated
#
$count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT count(*) FROM gene_expression " .
   "WHERE gene_id = ".$gene->dbID())->[0]->[0];

is($count, 1, "Found one entry in gene_expression for ENSG00000131044");

#
# make sure the tissue table was updated
#
$count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT count(*) FROM tissue " .
   "WHERE name = 'test_name2'")->[0]->[0];
is($count, 1, "Found one entry in tissue table for test_name2");

@expressions = @{$ea->fetch_all_by_Gene($gene)};
note "expressions: " . scalar(@expressions);
is(@expressions, 1, "One expression for ENSG00000131044");

@expressions = @{$ea->fetch_all_by_Gene(undef)};
is(@expressions, 1, "Fetched one gene expression");



#
# test the removal of this expression
#
$ea->remove_from_Gene($gene);
$count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT count(*) FROM gene_expression " .
   "WHERE gene_id = " . $gene->dbID())->[0]->[0];

is($count, 0, "No gene expressions fetched");



$multi->restore('core', 'gene_expression', 'tissue');

done_testing();
