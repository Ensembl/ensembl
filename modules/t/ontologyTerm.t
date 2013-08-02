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

# switch on the debug prints
our $verbose = 0;

debug("Startup test");
ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('ontology');

my $odb = $multi->get_DBAdaptor("ontology");
debug("Ontology database instatiated");
ok($odb);

my $human = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $human->get_DBAdaptor("core");
debug("Test database instatiated");
ok($db);

my $accession = "GO:0003677";
my $go_adaptor = $odb->get_OntologyTermAdaptor();
my $term = $go_adaptor->fetch_by_accession($accession);
ok(!$term->is_root, "Term is not a root");

my $gene;
my $ga = $db->get_GeneAdaptor();

my $genes = $ga->fetch_all_by_GOTerm($term);
is(@{$genes}, 2, "Genes match the GO term");

my $pattern = '%binding%';
my $terms = $go_adaptor->fetch_all_by_name($pattern);
is(@{$terms}, 138, "Found binding terms");

my $roots = $go_adaptor->fetch_all_roots();
is(@{$roots}, 1, "Found roots");

my $go_roots = $go_adaptor->fetch_all_roots('GO');
is(@{$go_roots}, 1, "Found go roots");

my $efo_roots = $go_adaptor->fetch_all_roots('EFO');
is(@{$efo_roots}, 0, "Found no efo roots");



done_testing();
