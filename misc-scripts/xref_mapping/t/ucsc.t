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
use Test::Warnings;

use FindBin '$Bin';

use Xref::Test::TestDB;

my $db = Xref::Test::TestDB->new();

use_ok 'XrefParser::UCSC_human_parser';

# add EntrezGene/WikiGene source to the db
$db->schema->populate(
  'Source',
  [
   [ qw/name/ ],
   [ 'UCSC_human' ],
  ]
);

my ($source_id, $species_id) = ( 1, 9606 );

my $parser = XrefParser::UCSC_human_parser->new($db->dbh);

isa_ok( $parser, 'XrefParser::UCSCParser' );

$parser->run({
  source_id  => $source_id,
  species_id => $species_id,
  files      => [ "$Bin/test-data/ucsc.txt" ],
});

# Test if all 10 entries were inserted
is($db->schema->resultset('CoordinateXref')->count, 10, 'All entries inserted');

my @xrefs = $db->schema->resultset('CoordinateXref')->search(
  {
   accession   => 'uc031tla.1',
   chromosome  => '1',
   txStart     => '17369',
   txEnd       => '17436',
   strand      => '-1',
   exonStarts  => '17369',
   exonEnds    => '17436',
   species_id  => $species_id,
  }
);

is(scalar @xrefs, 1, 'Xref uc031tla.1 exists and is all as expected');

done_testing();
