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
use FindBin '$Bin';

use Xref::Test::TestDB;
use XrefParser::VGNCParser;

use_ok 'Xref::Test::TestDB';

my $db = Xref::Test::TestDB->new();

my $config = $db->config;


my $species = $db->schema->resultset('Species')->create({
  species_id  => 9598,
  taxonomy_id => 9598,
  name        => 'Pan troglodytes'
});

my $vgnc_source = $db->create_db_row('Source',{
  name => 'VGNC',
  priority_description => 'Test VGNC source',
  priority => 10
});


use_ok 'XrefParser::VGNCParser';

my $parser = XrefParser::VGNCParser->new($db->dbh);

isa_ok( $parser, 'XrefParser::VGNCParser' );

$parser->run( {
 source_id  => $vgnc_source->source_id,
 species_id => 9598,
 files      => ["$Bin/test-data/vgnc.txt"],
 verbose    => 1
} );


# Test if all the rows were inserted
is($db->schema->resultset('Xref')->count, 6, "All 6 rows were inserted");

ok(
  $db->schema->resultset('Xref')->check_direct_xref({
    accession   => 'VGNC:14659',
    label       => 'CYYR1',
    description => 'cysteine and tyrosine rich 1',
    source_id   => $vgnc_source->source_id,
    species_id  => 9598
  }),
 'Sample chimpanzee direct Xref has been inserted'
);

is($db->schema->resultset('Synonym')->count, 1, "Sample synonym was inserted");

my $parser_no_file = XrefParser::VGNCParser->new($db->dbh);

throws_ok{
  $parser_no_file->run( {
   source_id  => $vgnc_source->source_id,
   species_id => 9598,
   files      => [],
  } );
} qr/No file name/, 'No file provided throws error' ;

done_testing();

