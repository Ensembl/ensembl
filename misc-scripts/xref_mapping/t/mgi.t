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

use strict;
use warnings;

use Test::More;
use Test::Exception;
use Test::Warnings;

use FindBin '$Bin';

use Xref::Test::TestDB;

my $db = Xref::Test::TestDB->new();

# populate the synonyms to test
$db->schema->populate( 'Synonym',
  [ [qw/xref_id synonym/], [ 1, 'Test-Synonym1' ], [ 1, 'Test-Synonym2' ] ] );

use_ok 'XrefParser::MGIParser';

my $parser = XrefParser::MGIParser->new($db->dbh);
isa_ok( $parser, 'XrefParser::MGIParser' );

$parser->run({
  source_id  => 55,
  species_id => 10090,
  files      => ["$Bin/test-data/MGI_mini.rpt"],
});

ok(
  $db->schema->resultset('Xref')->check_direct_xref(
    {
      accession   => 'MGI:1915733',
      label       => '1110002O04Rik',
      description => 'RIKEN cDNA 1110002O04 gene',
      source_id   => 55,
      species_id  => 10090,
      info_type   => 'DIRECT'
    }
  ),
  'Sample mouse direct Xref has been inserted'
);

ok(
  $db->schema->resultset('GeneDirectXref')->find(
    {
      ensembl_stable_id => "ENSMUSG00000102531"
    }
  ),
  'Sample mouse gene direct Xref has been inserted'
);


is($db->schema->resultset('Xref')->count, 10, "All 10 rows were inserted in to Xref");
is($db->schema->resultset('GeneDirectXref')->count, 10, "All 10 rows were inserted in to GeneDirectXref");
is($db->schema->resultset('Synonym')->count, 2, "Synonym count remained the same");


done_testing();
