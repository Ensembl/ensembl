use strict;
use warnings;

use Test::More;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::OntologyDBAdaptor;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::DnaDnaAlignFeature;

# switch on the note prints
our $verbose = 0;

note("Startup test");
ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('ontology');

my $odb = $multi->get_DBAdaptor("ontology");
note("Ontology database instatiated");
ok($odb);

my $human = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $human->get_DBAdaptor("core");
note("Test database instatiated");
ok($db);
my $go_adaptor = $odb->get_OntologyTermAdaptor();

my $accession = "GO:0000217";
my $term = $go_adaptor->fetch_by_accession($accession);
is($term, undef, "GO:0000217 does not exist in the non obsolete list");
$term = $go_adaptor->fetch_by_accession($accession, 1);
ok($term->is_obsolete, "GO:0000217 is obsolete");

$accession = "GO:0003698";
$term = $go_adaptor->fetch_by_accession($accession);
ok($term->name, "GO:0003698 alt_id was fetched");

$accession = "GO:0003677";
$term = $go_adaptor->fetch_by_accession($accession, 1);
ok(!$term->is_obsolete, "GO:0003677 is not obsolete");
ok(!$term->is_root, "GO:0003677 is not a root");

my $gene;
my $ga = $db->get_GeneAdaptor();

my $genes = $ga->fetch_all_by_GOTerm($term);
is(@{$genes}, 2, "Genes match the GO term");

my $pattern = '%binding%';
my $terms = $go_adaptor->fetch_all_by_name($pattern);
is(@{$terms}, 134, "Found binding terms");

$terms = $go_adaptor->fetch_all_by_name($pattern, undef, 1);
is(@{$terms}, 138, "Found binding terms, including obsolete ones");

my $roots = $go_adaptor->fetch_all_roots();
is(@{$roots}, 1, "Found roots");

my $go_roots = $go_adaptor->fetch_all_roots('GO');
is(@{$go_roots}, 1, "Found go roots");

my $efo_roots = $go_adaptor->fetch_all_roots('EFO');
is(@{$efo_roots}, 0, "Found no efo roots");

#Now go back to the OntologyXref & see if we can do the reverse lookup
{
  my $go_db_links = $genes->[0]->get_all_DBLinks('GO');
  foreach my $dbentry (@{$go_db_links}) {
    my $term = $dbentry->get_OntologyTerm();
    my $direct_term = $go_adaptor->fetch_by_accession($dbentry->primary_id());
    is_deeply($term, $direct_term, 'Fetching the OntologyTerm from the OntologyXref should match the same object from OntologyTermAdaptor');
  }
}
$term = $go_adaptor->fetch_by_accession('GO:0000182',1); # unintentionally picked an obsolete term for testing on.
my $term_list = $go_adaptor->fetch_all_by_descendant_term($term);
my $inclusive_term_list = $go_adaptor->fetch_all_by_descendant_term($term,undef,undef,1);
ok (scalar(@$term_list) == scalar(@$inclusive_term_list) - 1, "Zero_distance flag on fetch_all_by_descendant_term");

done_testing();
