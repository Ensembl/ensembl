=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

my $db = Xref::Test::TestDB->new();
my %config = %{ $db->config };

use_ok 'XrefParser::EntrezGeneParser';

# add EntrezGene/WikiGene source to the db
$db->schema->populate(
  'Source',
  [
   [ qw{ name } ],
   [ 'EntrezGene' ],
   [ 'WikiGene' ],
  ]
);

my ($source_id, $species_id) = ( 1, 9606 );

throws_ok {
  XrefParser::EntrezGeneParser->new($db->dbh)->run({
    species_id => $species_id,
    files      => [ "$Bin/test-data/entrez-gene_info" ],
  })
} qr/source_id/, 'Run without source';

throws_ok {
  XrefParser::EntrezGeneParser->new($db->dbh)->run({
    source_id  => $source_id,
    files      => [ "$Bin/test-data/entrez-gene_info" ],
  })
} qr/species_id/, 'Run without species_id';

throws_ok {
  XrefParser::EntrezGeneParser->new($db->dbh)->run({
    source_id  => $source_id,
    species_id => $species_id,
  })
} qr/files/, 'Run without files';


my $parser = XrefParser::EntrezGeneParser->new($db->dbh);

isa_ok( $parser, 'XrefParser::EntrezGeneParser' );

$parser->run({
 source_id  => $source_id,
 species_id => $species_id,
 files      => [ "$Bin/test-data/entrez-gene_info" ],
  });

my @xrefs = $db->schema->resultset('Xref')->search(
  {
   accession   => "10",
   label       => 'NAT2',
   description => 'N-acetyltransferase 2',
   species_id  => $species_id,
  }
);

# parser insert both EntrezGene and WikiGene xrefs
is( scalar @xrefs, 2, 'Sample EntrezGene/WikiGene xrefs inserted' );

# check synonyms have been inserted as well
foreach my $syn ( qw/ AAC2 NAT-2 PNAT / ) {
  ok( $db->schema->resultset('Synonym')->search( { synonym => $syn } ), 'Synonym inserted' );
}

@xrefs = $db->schema->resultset('Xref')->search(
  {
   accession   => "13",
   label       => 'AADAC',
   description => 'arylacetamide deacetylase',
   species_id  => $species_id,
  }
);

is( scalar @xrefs, 2, 'Sample EntrezGene/WikiGene xrefs inserted' );

foreach my $syn ( qw/ CES5A1 DAC / ) {
  ok( $db->schema->resultset('Synonym')->search( { synonym => $syn } ), 'Synonym inserted' );
}

# Test if all 10 entries were inserted (twice)
is($db->schema->resultset('Xref')->count, 20, "All entries inserted");

done_testing();

