use strict;

use Bio::EnsEMBL::Test::TestUtils;

use Test::More;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Expression;
use Bio::EnsEMBL::Analysis;

our $verbose = 1;
our $clean = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;

# get a core DBAdaptor
my $db = $multi->get_DBAdaptor("core");

my $sa = $db->get_SliceAdaptor();
my $exa = $db->get_ExonAdaptor();
my $ga = $db->get_GeneAdaptor();
my $ta = $db->get_TranscriptAdaptor();


#
# Test get_ExpressionAdaptor works
#
my $ea = $db->get_ExpressionAdaptor();

ok($ea && ref($ea) && $ea->isa('Bio::EnsEMBL::DBSQL::ExpressionAdaptor'));


# hide the contents of the attrib_type, misc_attrib, seq_region_attrib tables
# so we can test storing etc. with a clean slate
$multi->hide('core', 'gene_expression', 'transcript_expression', 'exon_expression', 'tissue');


#################
# Gene functionality tests
#

my $analysis = Bio::EnsEMBL::Analysis->new(-LOGIC_NAME => 'test_ana2');
my $aa = $db->get_AnalysisAdaptor();
$aa->store($analysis);

my $expression = Bio::EnsEMBL::Expression->new(-NAME    => 'test_name2',
                                           -ONTOLOGY    => 'test_ontology2',
                                           -DESCRIPTION => 'test_desc2',
                                           -ANALYSIS    =>  $analysis,
                                           -VALUE_TYPE  => 'count',
                                           -VALUE       => 'test_value2');


my $gene = $ga->fetch_by_stable_id('ENSG00000131044');

$ea->store_on_Gene($gene, [$expression]);

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

@expressions = @{$ea->fetch_all_by_Gene($gene, "rubbish")};
is(@expressions, 0, "Fetched no expression for ENSG00000131044 and tissue rubbish");

@expressions = @{$ea->fetch_all_by_Gene($gene, "test_name2")};
is(@expressions, 1, "Fetched one expression for ENSG00000131044 and tissue test_name2");

@expressions = @{$ea->fetch_all_by_Gene($gene, undef, "rubbish")};
is(@expressions, 0, "Fetched no expression for ENSG00000131044 and analysis rubbish");

@expressions = @{$ea->fetch_all_by_Gene($gene, undef, "test_ana2")};
is(@expressions, 1, "Fetched one expression for ENSG00000131044 and analysis test_ana2");

@expressions = @{$ea->fetch_all_by_Gene($gene, undef, undef, "RPKM")};
is(@expressions, 0, "Fetched no expression for ENSG00000131044 and value type RPKM");

@expressions = @{$ea->fetch_all_by_Gene($gene, undef, undef, "count")};
is(@expressions, 1, "Fetched one expression for ENSG00000131044 and value type count");

@expressions = @{$ea->fetch_all_by_Gene(undef,"test_name2")};
is(@expressions, 1, "Fetched one expression for tissue test_name2");


$expression = $expressions[0];

note("Fetching values from stored expression");
is($expression->name, 'test_name2', "Tissue name is test_name2");
is($expression->ontology, 'test_ontology2', "Tissue ontology is test_ontology2");
is($expression->description, 'test_desc2', "Tissue description is test_desc2");
is($expression->analysis->logic_name, 'test_ana2', "Tissue analysis is test_ana2");
is($expression->value, 'test_value2', "Expression value is test_value2");
is($expression->value_type, 'count', "Expression value type is count");


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
# try to add an expression with an already existing tissue
#
note("Storing new expression with already existing tissue");
$ea->store_on_Gene($gene, [$expression]);
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
is(@expressions, 1, "One expression for ENSG00000131044");

@expressions = @{$ea->fetch_all_by_Gene(undef)};
is(@expressions, 1, "Fetched one gene expression");



#
# test the removal of this expression
#
note("Removing all expressions");
$ea->remove_from_Gene($gene);
$count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT count(*) FROM gene_expression " .
   "WHERE gene_id = " . $gene->dbID())->[0]->[0];

is($count, 0, "No gene expressions fetched");


#################
# Transcript functionality tests
#

my $transcript = $ta->fetch_by_stable_id('ENST00000355555');

$ea->store_on_Transcript($transcript, [$expression]);

#
# make sure the transcript_expression table was updated
#
my $count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT count(*) FROM transcript_expression " .
   "WHERE transcript_id = ".$transcript->dbID())->[0]->[0];

is($count, 1, "Stored a transcript expression for ENST00000355555");

#
# test that we can now retrieve this expression
#
my @expressions = @{$ea->fetch_all_by_Transcript($transcript)};
is(@expressions, 1, "Fetched one expression for ENST00000355555");

#
# test the removal of this expression with tissue code
#
note("Removing test_name2 from transcript");
$ea->remove_from_Transcript($transcript, "test_name2");
$count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT count(*) FROM transcript_expression " .
   "WHERE transcript_id = " . $transcript->dbID())->[0]->[0];

is($count, 0, "No entry in transcript_expression table for transcript ENST00000355555 and tissue test_name2");


#################
# Exon functionality tests
#

my $exon = $exa->fetch_by_stable_id('ENSE00000859937');
print $exon . " fetched transcript\n" ;

$ea->store_on_Exon($exon, [$expression]);

#
# make sure the exon_expression table was updated
#
my $count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT count(*) FROM exon_expression " .
   "WHERE exon_id = ".$exon->dbID())->[0]->[0];

is($count, 1, "Stored a exon expression for ENSE00000859937");

#
# test that we can now retrieve this expression
#
my @expressions = @{$ea->fetch_all_by_Exon($exon)};
is(@expressions, 1, "Fetched one expression for ENSE00000859937");

#
# test the removal of this expression with tissue code
#
note("Removing test_name2 from exon");
$ea->remove_from_Exon($exon, "test_name2");
$count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT count(*) FROM exon_expression " .
   "WHERE exon_id = " . $exon->dbID())->[0]->[0];

is($count, 0, "No entry in exon_expression table for exon ENSE00000859937 and tissue test_name2");




$multi->restore('core', 'gene_expression', 'transcript_expression', 'exon_expression', 'tissue');

done_testing();
