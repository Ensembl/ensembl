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
package Bio::EnsEMBL::Xref::Test::Parser::MIMParser;

use strict;
use warnings;

use Test::More;
use Test::Exception;
#FIXME
use Test::Warnings 'allow_warnings';

use English '-no_match_vars';
use FindBin '$Bin';
use Readonly;

use Xref::Test::TestDB;

use XrefParser::MIMParser;


Readonly my $SOURCE_ID_OMIM        => 60;
Readonly my $SOURCE_ID_OMIM_GENE   => 61;
Readonly my $SOURCE_ID_OMIM_MORBID => 62;
Readonly my $SPECIES_ID_HUMAN      => 9606;

Readonly my $TI_CONTENT => << 'EOF';
#100100 PRUNE BELLY SYNDROME; PBS
;;ABDOMINAL MUSCLES, ABSENCE OF, WITH URINARY TRACT ABNORMALITY AND
CRYPTORCHIDISM;;
EAGLE-BARRETT SYNDROME; EGBRS;
EOF

# The test file contains 3 phenotype-only entries, 1 gene-only entry
# and 1 gene/phenotype entry
Readonly my $NUMBER_OF_GENE_XREFS => 2;
Readonly my $NUMBER_OF_MORBID_XREFS => 4;

# The test file contains 3 "MOVED TO" entries. Of those, two form a
# chain ultimately pointing to a gene/phenonype entry, meaning they
# should appear independently for both sources. Conversely, the last
# one points to a "REMOVED FROM DATABASE" entry and thus should not be
# inserted at all.
Readonly my $NUMBER_OF_SYNONYMS => 4;


my $db = Xref::Test::TestDB->new();

my $config = $db->config();


my $parser;

subtest 'extract_ti()' => sub {

  # First some invalid records...

  my $TI_TRUNCATED = << 'EOF'
*FIELD* NO
100100
*FIELD* TI
100100 PRU
EOF
    ;
  ok( ! defined XrefParser::MIMParser::extract_ti( $TI_TRUNCATED ),
      'Reports failure on record truncated mid-TI'
   );

  my $NO_TI = << 'EOF'
*FIELD* NO
100100
*FIELD* TX
fnord
*RECORD*
EOF
    ;
  ok( ! defined XrefParser::MIMParser::extract_ti( $NO_TI ),
      'Reports failure on record with no (recognisable) TI field'
   );

  # ...and then some which should work

  my $TI_MIDDLE_OF_RECORD = << "EOF"
*FIELD* NO
100100
*FIELD* TI
$TI_CONTENT
*FIELD* TX
fnord
EOF
    ;
  is( XrefParser::MIMParser::extract_ti( $TI_MIDDLE_OF_RECORD ),
      $TI_CONTENT,
      'Extracts TI from field in the middle of record'
   );

  my $TI_END_OF_RECORD = << "EOF"
*FIELD* NO
100100
*FIELD* TI
$TI_CONTENT
*RECORD*
EOF
    ;
  is( XrefParser::MIMParser::extract_ti( $TI_END_OF_RECORD ),
      $TI_CONTENT,
      'Extracts TI from field at the end of record'
   );

  my $TI_END_OF_FINAL_RECORD = << "EOF"
*FIELD* NO
100100
*FIELD* TI
$TI_CONTENT
*THEEND*
EOF
    ;
  is( XrefParser::MIMParser::extract_ti( $TI_END_OF_FINAL_RECORD ),
      $TI_CONTENT,
      'Extracts TI from field at the end of file'
   );

};


subtest 'parse_ti()' => sub {
  my ( $type_symbol, $omim_number, $description );

  ( $type_symbol, $omim_number, $description )
    = XrefParser::MIMParser::parse_ti( $TI_CONTENT );
  is( $type_symbol, q{#}, 'Extracted type symbol' );
  is( $omim_number, '100100', 'Extracted OMIM ID' );
  like( $description, qr{ \A PRUNE[ ]BELLY .+ EGBRS; }msx,
        'Description begins and ends where expected' );
  like( $description, qr{ PBS;;ABDO .+ DISM;;EAGLE }msx,
        'Newlines removed from mid-description' );
  like( $description, qr{ ABNORMALITY[ ]AND[ ]CRYPTORCHIDISM }msx,
        'Mid-description newline removal preserves word boundaries' );

  ( $type_symbol, $omim_number, $description )
    = XrefParser::MIMParser::parse_ti(
      '100050 AARSKOG SYNDROME'
    );
  is( $type_symbol, q{}, 'Correctly handled null type symbol' );

  ( $type_symbol )
    = XrefParser::MIMParser::parse_ti(
      'AARSKOG SYNDROME'
    );
  is( $type_symbol, undef, 'Reported error on missing OMIM ID');

  ( $type_symbol, $omim_number, $description )
    = XrefParser::MIMParser::parse_ti(
      '_100050 AARSKOG'
    );
  is( $type_symbol, undef, 'Reported error on invalid type symbol');

  ( $type_symbol )
    = XrefParser::MIMParser::parse_ti(
      '100050'
    );
  is( $type_symbol, undef, 'Reported error on missing description');

};


subtest 'Missing required source IDs' => sub {

  $parser = XrefParser::MIMParser->new($db->dbh);

  # As BaseParser stands in August 2019, attempting to retrieve source
  # ID for a source that is not present in the database merely
  # produces a warning; it is up to the parsers themselves to abort.
  # Seeing as the whole point of this test is to confirm that we *do*
  # abort, simply ignore the warnings.
  allow_warnings( 1 );
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_OMIM,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/omim-mini.txt" ],
  }); },
             qr{ \A Failed[ ]to[ ]retrieve[ ]MIM[ ]source[ ]IDs }msx,
             'Throws on source IDs missing from DB' );
  # IMPORTANT, this does not reset itself upon end of current scope
  allow_warnings( 0 );

};


# We will need these later
$db->schema->populate( 'Source', [
  [ qw{ source_id name } ],
  [ $SOURCE_ID_OMIM_GENE,   'MIM_GENE' ],
  [ $SOURCE_ID_OMIM_MORBID, 'MIM_MORBID' ],
] );


subtest 'Problems with input file' => sub {

  $parser = XrefParser::MIMParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_OMIM,
    species_id => $SPECIES_ID_HUMAN,
    files      => [],
  }); },
             qr{ \A No[ ]file[ ]name[ ] }msx,
             'Throws on no file name' );

  $parser = XrefParser::MIMParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_OMIM,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/omim-NONEXISTENT.txt" ],
  }); },
             qr{\A Could[ ]not[ ]find[ ]either[ ]' }msx,
             'Throws on missing input file' );

};

subtest 'Malformed input' => sub {

  # Most of this is already taken care of by extract_ti() and
  # parse_ti() tests, here we just make sure the parser actually
  # aborts on errors from those two.

  $parser = XrefParser::MIMParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_OMIM,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/omim-truncated.txt" ],
  }); },
             qr{ \A Failed[ ]to[ ]extract[ ]TI[ ]field[ ] }msx,
             'Throws on failure to extract TI field from record' );

  $parser = XrefParser::MIMParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_OMIM,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/omim-badTIfield.txt" ],
  }); },
             qr{ \A Failed[ ]to[ ]extract[ ]record[ ]type[ ] }msx,
             'Throws on failure to parse extracted TI field' );

};

subtest 'Successful inserts' => sub {

  $parser = XrefParser::MIMParser->new($db->dbh);
  isa_ok( $parser, 'XrefParser::MIMParser', 'Instantiated OMIM parser' );
  lives_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_OMIM,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/omim-mini.txt" ],
  }); }, 'Parsed OMIM data without errors' );
};

subtest 'Xrefs' => sub {

  is( $db->schema->resultset('Xref')->count({
        source_id => $SOURCE_ID_OMIM_GENE
      }),
      $NUMBER_OF_GENE_XREFS,
      'All OMIM Gene Map xrefs have been inserted' );
  is( $db->schema->resultset('Xref')->count({
        source_id => $SOURCE_ID_OMIM_MORBID
      }),
      $NUMBER_OF_MORBID_XREFS,
      'All OMIM Morbid Map xrefs have been inserted' );

  ok(
    $db->schema->resultset('Xref')->check_direct_xref({
      accession   => '100640',
      label       => 'ALDEHYDE DEHYDROGENASE 1 FAMILY, MEMBER A1; ALDH1A1 [*100640]',
      description => 'ALDEHYDE DEHYDROGENASE 1 FAMILY, MEMBER A1; ALDH1A1;;'
        . 'ALDEHYDE DEHYDROGENASE 1; ALDH1;;'
        . 'ACETALDEHYDE DEHYDROGENASE 1;;'
        . 'ALDH, LIVER CYTOSOLIC;;'
        . 'RETINAL DEHYDROGENASE 1; RALDH1',
      source_id   => $SOURCE_ID_OMIM_GENE,
      species_id  => $SPECIES_ID_HUMAN,
      info_type   => 'UNMAPPED',
    }),
    'Sample OMIM Gene Map entry inserted as expected'
  );

  ok(
    $db->schema->resultset('Xref')->check_direct_xref({
      accession   => '200150',
      label       => 'CHOREOACANTHOCYTOSIS; CHAC [+200150]',
      description => 'CHOREOACANTHOCYTOSIS; CHAC;;'
        . 'LEVINE-CRITCHLEY SYNDROME;;'
        . 'ACANTHOCYTOSIS WITH NEUROLOGIC DISORDER;;'
        . 'NEUROACANTHOCYTOSIS;;'
        . 'CHOREA-ACANTHOCYTOSIS',
      source_id   => $SOURCE_ID_OMIM_MORBID,
      species_id  => $SPECIES_ID_HUMAN,
      info_type   => 'UNMAPPED',
    }),
    'Sample OMIM Morbid Map entry inserted as expected'
  );

};

subtest 'Synonyms' => sub {
  my $rs;
  my $matching_xref;

  is( $db->schema->resultset('Synonym')->count, $NUMBER_OF_SYNONYMS,
      'Expected number of synonyms' );

  $rs = $db->schema->resultset('Xref')->search({
    accession => '200150',
  });
  $matching_xref = $rs->next;
  ok( $db->schema->resultset('Synonym')->check_synonym({
        xref_id => $matching_xref->xref_id,
        synonym => '100500',
      }),
      'Synonym from beginning of a chain of renames points at the real entry'
    );
  ok( $db->schema->resultset('Synonym')->check_synonym({
        xref_id => $matching_xref->xref_id,
        synonym => '100650',
      }),
      'Synonym from the middle of a chain of renames points at the real entry'
    );

  is ( $db->schema->resultset('Synonym')->count({
         synonym => '100680',
       }), 0,
       'Synonyms not created for ultimately removed entries'
     );

};

subtest 'Replay safety' => sub {

  $parser = XrefParser::MIMParser->new($db->dbh);
  lives_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_OMIM,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/omim-mini.txt" ],
  }); }, 'Re-parsed OMIM data without errors' );

  is( $db->schema->resultset('Xref')->count({
        source_id => $SOURCE_ID_OMIM_GENE
      }),
      $NUMBER_OF_GENE_XREFS,
      'No new OMIM Gene Map xrefs inserted by the replay' );

  is( $db->schema->resultset('Xref')->count({
        source_id => $SOURCE_ID_OMIM_MORBID
      }),
      $NUMBER_OF_MORBID_XREFS,
      'No new OMIM Morbid Map xrefs inserted by the replay' );

  is( $db->schema->resultset('Synonym')->count, $NUMBER_OF_SYNONYMS,
      'No new synonyms inserted by the replay' );

  # Ideally we would also make sure the replay has not modified
  # existing entries, no quick way of doing so though.

};


done_testing();


1;
