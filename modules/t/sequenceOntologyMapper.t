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

# Module compiles
ok(1, 'Bio::EnsEMBL::Utils::SequenceOntologyMapper compiles');

#
# Test constructor
#
throws_ok { new Bio::EnsEMBL::Utils::SequenceOntologyMapper() } 
  qr /No ontology term adaptor/, "constructor with no arg raises exception";
throws_ok { new Bio::EnsEMBL::Utils::SequenceOntologyMapper(bless({}, 'Someclass')) } 
  qr /Argument is not/, "constructor with arg not of class OntologyTermAdaptor raises exception";

my $omulti = Bio::EnsEMBL::Test::MultiTestDB->new('ontology');
my $odb = $omulti->get_DBAdaptor('ontology');
my $oa = $odb->get_OntologyTermAdaptor();

my $mapper = Bio::EnsEMBL::Utils::SequenceOntologyMapper->new($oa);
isa_ok($mapper, 'Bio::EnsEMBL::Utils::SequenceOntologyMapper');

#
# Test to_accession/to_name methods
#
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db = $multi->get_DBAdaptor('core');

my $mappings =
  [
   # test generic feature
   { obj => Bio::EnsEMBL::Feature->new(), accession => 'SO:0000001', name => 'region' },
   # test various genes with/without biotypes
   { obj => Bio::EnsEMBL::Gene->new(), accession => 'SO:0000704', name => 'gene' },
   { obj => Bio::EnsEMBL::Gene->new(-biotype => 'protein_coding'), accession => 'SO:0001217', name => 'protein_coding_gene' },
   { obj => Bio::EnsEMBL::Gene->new(-biotype => 'tRNA'), accession => 'SO:0001272', name => 'tRNA_gene' },
   # test various transcripts with/without biotypes
   { obj => Bio::EnsEMBL::Transcript->new(), accession => 'SO:0000673', name => 'transcript' },
   { obj => Bio::EnsEMBL::Transcript->new(-biotype => 'processed_transcript'), accession => 'SO:0001503', name => 'processed_transcript' },
   { obj => Bio::EnsEMBL::Transcript->new(-biotype => 'retrotransposed'), accession => 'SO:0000569', name => 'retrotransposed' },
   # exons
   { obj => Bio::EnsEMBL::Exon->new, accession => 'SO:0000147', name => 'exon' },
   # slices
   { obj => $db->get_SliceAdaptor->fetch_by_region('chromosome', '20', 30_270_000, 31_200_000), accession => 'SO:0000001', name => 'region' },
   # simple features
   { obj => Bio::EnsEMBL::SimpleFeature->new(), accession => 'SO:0001411', name => 'biological_region' },
   # misc features
   { obj => Bio::EnsEMBL::MiscFeature->new(), accession => 'SO:0001411', name => 'biological_region' },
   # repeat feature
   { obj => Bio::EnsEMBL::RepeatFeature->new(), accession => 'SO:0000657', name => 'repeat_region' },
  ];

throws_ok { $mapper->translate( bless({}, 'Someclass')) } 
  qr /not found/, "translate method raises exception";

map { 
  is($mapper->to_accession($_->{obj}), $_->{accession}, sprintf "%s maps to correct accession", ref($_->{obj})) and
    is($mapper->to_name($_->{obj}), $_->{name}, sprintf "%s maps to correct", ref($_->{obj})) }
  @{$mappings};

done_testing();

