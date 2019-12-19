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
use lib "$Bin/";

use Xref::Test::TestDB;

my $db     = Xref::Test::TestDB->new();
my %config = %{ $db->config };

use_ok 'XrefParser::HPAParser';

# initialize the parser
my $parser = XrefParser::HPAParser->new($db->dbh);

isa_ok( $parser, 'XrefParser::HPAParser' );

$parser->run({
  source_id  => 11,
  species_id => 9601,
  files      => ["$Bin/test-data/hpa.txt"],
});

ok(
    $db->schema->resultset('Xref')->check_direct_xref(
        {
            accession  => "1",
            label      => 'CAB000001',
            info_type  => 'DIRECT',
            source_id  => 11,
            species_id => 9601
        }
    ),
    'Sample hpa direct Xref has been inserted'
);

# Test if all the distinct rows were inserted in to xref
is( $db->schema->resultset('Xref')->count,
    2, "2 rows with distinct accessions were inserted" );

# Test if all the rows were inserted in to translation_direct_xref
is( $db->schema->resultset('TranslationDirectXref')->count,
    10, "10 rows were inserted in to translation_direct_xref" );

done_testing();

