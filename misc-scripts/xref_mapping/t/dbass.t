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
package Bio::EnsEMBL::Xref::Test::Parser::DBASSParser;

use strict;
use warnings;

use Test::More;
use Test::Exception;
use Test::Warnings 'warnings';

use FindBin '$Bin';
use Readonly;

use Xref::Test::TestDB;

use XrefParser::DBASSParser;


Readonly my $SOURCE_ID_DBASS3       => 18;
Readonly my $SOURCE_ID_DBASS5       => 19;
Readonly my $SPECIES_ID_HUMAN       => 9606;

Readonly my $NUMBER_OF_MAPPED_XREFS => 18;
Readonly my $NUMBER_OF_SYNONYMS     => 7;


my $db = Xref::Test::TestDB->new();

my $config = $db->config();


my $parser;

subtest 'Problems with input file' => sub {

  $parser = XrefParser::DBASSParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_DBASS3,
    species_id => $SPECIES_ID_HUMAN,
    files      => [],
  }); },
             qr{ \A No[ ]file[ ]name[ ] }msx,
             'Throws on no file name' );

  $parser = XrefParser::DBASSParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_DBASS3,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/dbass-NONEXISTENT.txt" ],
  }); },
             qr{\A Could[ ]not[ ]find[ ]either[ ]' }msx,
             'Throws on no input file missing' );

};

subtest 'Malformed header' => sub {
  my $QR_MALFORMED_HEADER = qr{ \A Malformed[ ]or[ ]unexpected[ ]header[ ] }msx;

  $parser = XrefParser::DBASSParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_DBASS3,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/dbass-badHeader-tooFewCols.txt" ],
  }); },
             $QR_MALFORMED_HEADER,
             'Throws on too few header columns' );

  $parser = XrefParser::DBASSParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_DBASS3,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/dbass-badHeader-tooManyCols.txt" ],
  }); },
             $QR_MALFORMED_HEADER,
             'Throws on too many header columns' );

  $parser = XrefParser::DBASSParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_DBASS3,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/dbass-badHeader-wrongName1.txt" ],
  }); },
             $QR_MALFORMED_HEADER,
             'Throws on wrong name of first column' );

  $parser = XrefParser::DBASSParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_DBASS3,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/dbass-badHeader-wrongName2.txt" ],
  }); },
             $QR_MALFORMED_HEADER,
             'Throws on wrong name of second column' );

  $parser = XrefParser::DBASSParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_DBASS3,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/dbass-badHeader-wrongName3.txt" ],
  }); },
             $QR_MALFORMED_HEADER,
             'Throws on wrong name of third column' );

  $parser = XrefParser::DBASSParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_DBASS3,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/dbass-badHeader-mixedUpEnds.txt" ],
  }); },
             $QR_MALFORMED_HEADER,
             'DBASS IDs and names for two different ends' );

};

subtest 'Malformed data' => sub {
  my $QR_MALFORMED_DATA
    = qr{ '[ ]has[ ]an[ ]incorrect[ ]number[ ]of[ ]columns[ ] }msx;

  $parser = XrefParser::DBASSParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_DBASS3,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/dbass-badData-tooFewCols.txt" ],
  }); },
             $QR_MALFORMED_DATA,
             'Throws on too few data columns' );

  $parser = XrefParser::DBASSParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_DBASS3,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/dbass-badData-tooManyCols.txt" ],
  }); },
             $QR_MALFORMED_DATA,
             'Throws on too many data columns' );

};

subtest 'Successful inserts' => sub {

  $parser = XrefParser::DBASSParser->new($db->dbh);
  isa_ok( $parser, 'XrefParser::DBASSParser', 'Instantiated DBASS3 parser' );
  lives_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_DBASS3,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/dbass-mini3.txt" ],
  }); }, 'Parsed DBASS3 data without errors' );

  $parser = XrefParser::DBASSParser->new($db->dbh);
  isa_ok( $parser, 'XrefParser::DBASSParser', 'Instantiated DBASS5 parser' );
  lives_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_DBASS5,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/dbass-mini5.txt" ],
  }); }, 'Parsed DBASS5 data without errors' );

};

subtest 'Xrefs' => sub {

  is( $db->schema->resultset('Xref')->count,
      $NUMBER_OF_MAPPED_XREFS,
      "All mapped and no unmapped xrefs have been inserted" );

  # We have chosen xrefs with synonyms to make sure they have been
  # correctly removed from respective labels
  # FIXME: this method is misnamed, it checks xrefs and not *direct* xrefs
  # FIXME: this isn't informative enough - it doesn't show WHERE the mismatch is if there is one
  ok(
     $db->schema->resultset('Xref')->check_direct_xref({
       accession   => '45',
       label       => 'FALDH',
       description => undef,
       source_id   => $SOURCE_ID_DBASS3,
       species_id  => $SPECIES_ID_HUMAN,
       info_type   => 'DIRECT',
     }),
     'Sample DBASS3 entry inserted as expected'
   );
  ok(
     $db->schema->resultset('Xref')->check_direct_xref({
       accession   => '206',
       label       => 'DAX1',
       description => undef,
       source_id   => $SOURCE_ID_DBASS5,
       species_id  => $SPECIES_ID_HUMAN,
       info_type   => 'DIRECT',
     }),
     'Sample DBASS5 entry inserted as expected'
   );

};

subtest 'Direct-xref links' => sub {
  my $rs;
  my $matching_xref;

  is( $db->schema->resultset('GeneDirectXref')->count,
      $NUMBER_OF_MAPPED_XREFS,
      "Number of xrefs and direct-xref links is the same" );

  $rs = $db->schema->resultset('Xref')->search({
    accession => '2',
  });
  $matching_xref = $rs->next;
  $rs = $db->schema->resultset('GeneDirectXref')->search({
    general_xref_id => $matching_xref->xref_id,
  });
  is( $rs->next->ensembl_stable_id, 'ENSG00000198691',
      'A mapped xref has a matching direct-xref link' );
};

subtest 'Synonyms' => sub {
  my $rs;
  my $matching_xref;

  is( $db->schema->resultset('Synonym')->count,
      $NUMBER_OF_SYNONYMS, "All synonyms have been inserted" );

  $rs = $db->schema->resultset('Xref')->search({
    accession => '2',
  });
  $matching_xref = $rs->next;
  ok(
     $db->schema->resultset('Synonym')->check_synonym({
       xref_id => $matching_xref->xref_id,
       synonym => 'ABCA4',
     }),
     'Correct synonym assignment for "foo (bar)" syntax'
   );

  $rs = $db->schema->resultset('Xref')->search({
    accession => '118',
  });
  $matching_xref = $rs->next;
  ok(
     $db->schema->resultset('Synonym')->check_synonym({
       xref_id => $matching_xref->xref_id,
       synonym => 'HCF-1A',
     }),
     'Correct synonym assignment for "foo(bar)" syntax'
   );

  $rs = $db->schema->resultset('Xref')->search({
    accession => '206',
  });
  $matching_xref = $rs->next;
  ok(
     $db->schema->resultset('Synonym')->check_synonym({
       xref_id => $matching_xref->xref_id,
       synonym => 'NR0B1',
     }),
     'Correct synonym assignment for "foo/bar" syntax'
   );
};

subtest 'Replay safety' => sub {

  $parser = XrefParser::DBASSParser->new($db->dbh);
  lives_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_DBASS3,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/dbass-mini3.txt" ],
  }); }, 'Re-parsed DBASS3 data without errors' );

  $parser = XrefParser::DBASSParser->new($db->dbh);
  lives_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_DBASS5,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/dbass-mini5.txt" ],
  }); }, 'Re-parsed DBASS5 data without errors' );

  is( $db->schema->resultset('Xref')->count,
      $NUMBER_OF_MAPPED_XREFS,
      "No new xrefs inserted by the replay" );

  is( $db->schema->resultset('GeneDirectXref')->count,
      $NUMBER_OF_MAPPED_XREFS,
      "No new direct-xref links inserted by the replay" );

  is( $db->schema->resultset('Synonym')->count,
      $NUMBER_OF_SYNONYMS, "No new synonyms inserted by the replay" );

  # Ideally we would also make sure the replay has not modified
  # existing entries, no quick way of doing so though.

};

done_testing();


1;
