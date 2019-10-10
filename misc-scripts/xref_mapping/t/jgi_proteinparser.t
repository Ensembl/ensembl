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

## no critic 'RequireFilenameMatchesPackage'
package Bio::EnsEMBL::Xref::Test::Parser::JGI_ProteinParser;


use strict;
use warnings;

use Test::More;
use Test::Exception;
use Test::Warnings;

use FindBin '$Bin';

use Xref::Test::TestDB;
use_ok 'XrefParser::JGI_ProteinParser';


my $db = Xref::Test::TestDB->new();

# add ciona species related sources to the db
$db->schema->populate(
  'Source',
  [
   [ qw/name/ ],
   [ 'cint_jgi_v1' ],
   [ 'cint_aniseed_v1' ],
  ]
);

my ($source_id, $species_id) = ( 1, 7719 );
my $parser;

throws_ok {
  $parser = XrefParser::JGI_ProteinParser->new($db->dbh);
  $parser->run({
    species_id => $species_id,
    files      => [ "$Bin/test-data/jgi_prot.fasta" ],
  });
} qr/source_id/, 'No source specified';

throws_ok {
  $parser = XrefParser::JGI_ProteinParser->new($db->dbh);
  $parser->run({
    source_id  => $source_id,
    files      => [ "$Bin/test-data/jgi_prot.fasta" ],
  });
} qr/species_id/, 'No species specified';

throws_ok {
  $parser = XrefParser::JGI_ProteinParser->new($db->dbh);
  $parser->run({
    source_id  => $source_id,
    species_id => $species_id,
  });
} qr/files/, 'No files specified';

$parser = XrefParser::JGI_ProteinParser->new($db->dbh);
isa_ok( $parser, 'XrefParser::JGI_ProteinParser' );

$parser->run({
  source_id  => $source_id,
  species_id => $species_id,
  files      => [ "$Bin/test-data/jgi_prot.fasta" ],
});

my @xrefs = $db->schema->resultset('Xref')->search(
  {
   accession   => "130003",
   label       => '130003',
   source_id   => $source_id,
   species_id  => $species_id
  }
);

is( scalar @xrefs, 1, 'Xref inserted' );

my $xref_id = $xrefs[0]->xref_id;

# sequence should have been inserted as well
my $primary_xref = $db->schema->resultset('PrimaryXref')->find(
  {
   xref_id => $xref_id
  }
);

my $sequence = 'MPPKKKKEVEKPPLILGRLGTSLKIGIVGLPNVGKSTFFNVLTKSEASAENFPFCTIDPNESRVPVPDER' .
  'WEFLCKYHKPASKVPAFLSVVDIAGLVKGANEGQGLGNAFLSHISGCDAIFHMTRAFDDAEVVHVEGDVN' .
  'PVRDLEIIQEELRLKDVEHLTKRLAELEKVYSRGGEKKYKLEFETLSKIKTLLVDEKKPVRDGEWGGKEI' .
  'EVLNEHLFLTSKPQIYLVNLSEKDYIRKKNKWLMKIKTWVTENDSSAILIPFSGAFELKLAEMADDAERK' .
  'AYLEEQYKDSVGSALSKIVVTGFKCLGLQYFFTAGADEVKAWTIKTGFLAPQAAGRIHTDFEKGFIMAEV' .
  'MKFSDFKELGSESAVKSAGKYRQQGRNYIVEDGDIIFFKFNTPSQPKKK*MSQLAEMADDAERKAYLEEQ' .
  'YKDSVGSALSKIVVTGFKCLGLQYFFTAGADEVKAWTIKTGFLAPQAAGRIHTDFEKGFIMAEVMKFSDF' .
  'KELGSESAVKSAGKYRQQGRNYIVEDGDIIFFKFNTPSQPKKK*';

ok($primary_xref, 'Xref sequence inserted');
is($primary_xref->sequence, $sequence, 'Correct xref sequence');
is($primary_xref->sequence_type, 'peptide', 'Correct xref sequence type');

# Test if all the rows were inserted
is($db->schema->resultset('Xref')->count, 10, "All 10 rows were inserted");

done_testing();

