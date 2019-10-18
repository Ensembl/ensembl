=head1 LICENSE

See the NOTICE file distributed with this work for additional information
regarding copyright ownership.
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
use Test::Warnings;

use FindBin '$Bin';

use Xref::Test::TestDB;

my $db = Xref::Test::TestDB->new();

use_ok 'XrefParser::XenopusJamboreeParser';

my $parser = XrefParser::XenopusJamboreeParser->new($db->dbh);

isa_ok( $parser, 'XrefParser::XenopusJamboreeParser' );

$parser->run({
  source_id  => 150,
  species_id => 8364,
  files      => ["$Bin/test-data/xenopusjamboree.txt"],
});

ok(
 $db->schema->resultset('Xref')->check_direct_xref(
  {
   accession   => "XB-GENE-478054",
   label       => 'trnt1',
   description => 'tRNA nucleotidyl transferase, CCA-adding, 1',
   source_id   => 150,
   species_id  => 8364
  }
 ),
 'Sample frog direct Xref has been inserted'
);

# Test if all the rows were inserted
is($db->schema->resultset('Xref')->count, 10, "All 10 rows were inserted");

# Test the description parsing
my $input_desc = "Putative ortholog of g2/mitotic-specific cyclin B3, 3 of 14";
my $expected_desc = "Putative ortholog of g2/mitotic-specific cyclin B3";
is(XrefParser::XenopusJamboreeParser::parse_description($input_desc), $expected_desc, "Desc parsing ok");

done_testing();

