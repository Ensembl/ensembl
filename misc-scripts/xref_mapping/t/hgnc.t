=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

use strict;
use warnings;
use Test::More;
use Test::Exception;
use FindBin '$Bin';

use Xref::Test::TestDB;
use Bio::EnsEMBL::Test::MultiTestDB;
use XrefParser::HGNCParser;

use_ok 'XrefParser::HGNCParser';

my $db = Xref::Test::TestDB->new();

my $config = $db->config;

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $core_dba = $multi_db->get_DBAdaptor('core');

my $species = $db->schema->resultset('Species')->create({
  species_id  => 9606,
  taxonomy_id => 9606,
  name        => 'homo_sapiens'
});

# populate used sources
$db->schema->resultset('Source')->populate([
  {
    source_id            => 45,
    name                 => 'HGNC',
    priority_description => 'ccds',
  },{
    source_id            => 46,
    name                 => 'HGNC',
    priority_description => 'entrezgene_manual',
  },{
    source_id            => 47,
    name                 => 'HGNC',
    priority_description => 'refseq_manual',
  },{
    source_id            => 48,
    name                 => 'HGNC',
    priority_description => 'ensembl_manual',
  },{
    source_id            => 49,
    name                 => 'HGNC',
    priority_description => 'desc_only',
  },{
    source_id            => 54,
    name                 => 'LRG_HGNC_notransfer',
  },{
    source_id            => 23,
    name                 => 'EntrezGene',
  },{
    source_id            => 146,
    name                 => 'WikiGene',
  }
]);

use_ok 'XrefParser::HGNCParser';

my $parser = XrefParser::HGNCParser->new($db->dbh);

isa_ok( $parser, 'XrefParser::HGNCParser' );

$parser->run_script( {
 source_id  => 46,
 species_id => 9606,
 file       => "$Bin/test-data/hgnc.tsv",
 dba        => $core_dba,
 verbose    => 1
} );

# Test if all the rows were inserted
is($db->schema->resultset('Xref')->count, 22, "All 22 HGNC xrefs were inserted");
is($db->schema->resultset('Synonym')->count, 42, "All 42 HGNC synonyms were inserted");

ok(
  $db->schema->resultset('Xref')->check_direct_xref({
    accession   => 'HGNC:49894',
    label       => 'AARSP1',
    description => 'alanyl-tRNA synthetase pseudogene 1',
    source_id   => 48,
    species_id  => 9606
  }),
 'Sample HGNC direct Xref has been inserted'
);

my $synomym_rs = $db->schema->resultset('Synonym')->search({
  xref_id => 2
});

my @synonyms;
my @expected_synonyms = qw( A1BG-AS A1BGAS FLJ23569 NCRNA00181 );
while ( my $synonym_result = $synomym_rs->next ) {
  push @synonyms, $synonym_result->synonym;
}

foreach my $test_synonym (sort @synonyms) {
  my $expected_synonym = shift @expected_synonyms;
  is( $test_synonym, $expected_synonym, "Synonym $expected_synonym of xref_id 2 is correct");
}

done_testing();

