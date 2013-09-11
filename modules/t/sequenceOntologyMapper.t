#######################################################
#                                                     #
# Test Bio::EnsEMBL::Utils::SequenceOntologyMapper,   #
# the mapper from Ensembl features object to sequence #
# ontology terms                                      #
#######################################################

use strict;
use warnings;

use Data::Dumper;
use Test::More;
use Test::Exception;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Scalar qw ( assert_ref );
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Utils::SequenceOntologyMapper;

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon; 
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::SimpleFeature; 
use Bio::EnsEMBL::MiscFeature; 
use Bio::EnsEMBL::RepeatFeature;
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::StructuralVariationFeature;
use Bio::EnsEMBL::Compara::ConstrainedElement; 
use Bio::EnsEMBL::Funcgen::RegulatoryFeature;

# Module compiles
ok(1, 'Bio::EnsEMBL::Utils::SequenceOntologyMapper compiles');

#
# Test constructor
#
throws_ok { new Bio::EnsEMBL::Utils::SequenceOntologyMapper() } 
  qr /No ontology term adaptor/, "constructor with no argument raises exception";
throws_ok { new Bio::EnsEMBL::Utils::SequenceOntologyMapper(bless({}, 'Someclass')) } 
  qr /Argument is not/, "constructor with object not of class OntologyTermAdaptor raises exception";

# my $registry = 'Bio::EnsEMBL::Registry';
# Bio::EnsEMBL::Registry->no_version_check(1);
# Bio::EnsEMBL::Registry->no_cache_warnings(1);

# $registry->load_registry_from_db(-host => '127.0.0.1',
# 				 -port => 33061,
# 				 -db_version => 73,
# 				 -user => 'ensro',
# 				 -no_cache => 1);
				 
# my $oa = $registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );
# assert_ref($oa, 'Bio::EnsEMBL::DBSQL::OntologyTermAdaptor');

my $omulti = Bio::EnsEMBL::Test::MultiTestDB->new('ontology');
my $odb = $omulti->get_DBAdaptor('ontology');
my $oa = $odb->get_OntologyTermAdaptor();

my $mapper = Bio::EnsEMBL::Utils::SequenceOntologyMapper->new($oa);
isa_ok($mapper, 'Bio::EnsEMBL::Utils::SequenceOntologyMapper');

#
# Test translate method
#
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db = $multi->get_DBAdaptor('core');

my $mappings =
  [
   { obj => Bio::EnsEMBL::Feature->new, term => 'region' },
   { obj => Bio::EnsEMBL::Gene->new, term => 'gene' },
   { obj => Bio::EnsEMBL::Transcript->new, term => 'transcript' },
   { obj => Bio::EnsEMBL::Exon->new, term => 'exon' },
   { obj => $db->get_SliceAdaptor->fetch_by_region('chromosome', '20', 30_270_000, 31_200_000), term => 'region' },
   { obj => Bio::EnsEMBL::SimpleFeature->new(), term => 'biological_region' },
   { obj => Bio::EnsEMBL::MiscFeature->new(), term => 'biological_region' },
   { obj => Bio::EnsEMBL::RepeatFeature->new(), term => 'repeat_region' },
   { obj => Bio::EnsEMBL::Variation::VariationFeature->new(), term => 'sequence_variant' },
   { obj => Bio::EnsEMBL::Variation::StructuralVariationFeature->new(), term => 'structural_variant' },
   { obj => Bio::EnsEMBL::Compara::ConstrainedElement->new(), term => 'DNA_constraint_sequence' }

   # 'Bio::EnsEMBL::Funcgen::RegulatoryFeature' => 'regulatory_region'
  ];

throws_ok { $mapper->translate( bless({}, 'Someclass')) } 
  qr /not found/, "translate method raises exception";

map { 
  print ref($_->{obj}) and is($mapper->translate($_->{obj}), $_->{term}, sprintf "%s translates as assumed", ref($_->{obj})) }
  @{$mappings};

done_testing();

