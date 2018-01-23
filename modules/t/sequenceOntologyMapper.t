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

#######################################################
#                                                     #
# Test Bio::EnsEMBL::Utils::SequenceOntologyMapper,   #
# the mapper from Ensembl features object to sequence #
# ontology terms                                      #
#######################################################

use strict;
use warnings;

use Test::More;
use Test::Warnings;
use Test::Exception;

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
   { obj => Bio::EnsEMBL::Gene->new(-biotype => 'dummy'), accession => 'SO:0000704', name => 'gene' }, # if we don't pass a non empty biotype it will default to protein_coding hence it will fail
   { obj => Bio::EnsEMBL::Gene->new(-biotype => 'protein_coding'), accession => 'SO:0001217', name => 'protein_coding_gene' },
   { obj => Bio::EnsEMBL::Gene->new(-biotype => 'tRNA'), accession => 'SO:0001263', name => 'ncRNA_gene' },
   # test various transcripts with/without biotypes
   { obj => Bio::EnsEMBL::Transcript->new(-biotype => 'dummy'), accession => 'SO:0000673', name => 'transcript' },
   { obj => Bio::EnsEMBL::Transcript->new(), accession => 'SO:0000234', name => 'mRNA' },
   { obj => Bio::EnsEMBL::Transcript->new(-biotype => 'processed_transcript'), accession => 'SO:0001877', name => 'lnc_RNA' },
   { obj => Bio::EnsEMBL::Transcript->new(-biotype => 'tRNA'), accession => 'SO:0000253', name => 'tRNA' },
   # exons
   { obj => Bio::EnsEMBL::Exon->new, accession => 'SO:0000147', name => 'exon' },
   # slices
   { obj => $db->get_SliceAdaptor->fetch_by_region('chromosome', '20', 30_270_000, 31_200_000), accession => 'SO:0000340', name => 'chromosome' },
   # simple features
   { obj => Bio::EnsEMBL::SimpleFeature->new(), accession => 'SO:0001411', name => 'biological_region' },
   # misc features
   { obj => Bio::EnsEMBL::MiscFeature->new(), accession => 'SO:0001411', name => 'biological_region' },
   # repeat feature
   { obj => Bio::EnsEMBL::RepeatFeature->new(), accession => 'SO:0000657', name => 'repeat_region' },
  ];

throws_ok { $mapper->to_accession( bless({}, 'Someclass')) } 
  qr /not found/, "to_accession raises exception";

map { 
  is($mapper->to_accession($_->{obj}), $_->{accession}, sprintf "%s maps to correct accession", ref($_->{obj})) and
    is($mapper->to_name($_->{obj}), $_->{name}, sprintf "%s maps to correct name", ref($_->{obj})) }
  @{$mappings};

done_testing();

