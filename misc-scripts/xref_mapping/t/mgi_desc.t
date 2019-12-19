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

use strict;
use warnings;

use Test::More;
use Test::Exception;
use Test::Warnings;

use FindBin '$Bin';

use Xref::Test::TestDB;

my $db     = Xref::Test::TestDB->new();

use_ok 'XrefParser::MGI_Desc_Parser';

my $parser = XrefParser::MGI_Desc_Parser->new($db->dbh);
isa_ok( $parser, 'XrefParser::MGI_Desc_Parser' );

$parser->run({
  source_id  => 58,
  species_id => 10090,
  files      => ["$Bin/test-data/MGIdesc_mini.rpt"],
});

# check xrefs
ok(
  $db->schema->resultset('Xref')->check_direct_xref(
    {
      accession   => 'MGI:1341858',
      label       => '03B03F',
      description => 'DNA segment, 03B03F (Research Genetics)',
      source_id   => 58,
      species_id  => 10090,
      info_type   => 'MISC'
    }
  ),
  'Sample mouse misc description Xref has been inserted'
);

is( $db->schema->resultset('Xref')->count, 10, "All 10 rows were inserted in to Xref" );

# check synonym
ok(
  $db->schema->resultset('Synonym')->check_synonym(
    {
      xref_id => 10,
      synonym => "Ecrg4"
    }
  ),
  'Synonym Ecrg4 has been inserted'
);

is( $db->schema->resultset('Synonym')->count, 2, "Inserted two synonyms" );

done_testing();
