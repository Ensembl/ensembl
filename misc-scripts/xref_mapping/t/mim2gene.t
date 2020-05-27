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

## no critic 'RequireFilenameMatchesPackage'
package Bio::EnsEMBL::Xref::Test::Parser::Mim2GeneParser;

use strict;
use warnings;

use Test::More;
use Test::Exception;
use Test::Warnings;

use English '-no_match_vars';
use FindBin '$Bin';

use Xref::Test::TestDB;

use XrefParser::Mim2GeneParser;


my $SOURCE_ID_ENTREZGENE  = 23;
my $SOURCE_ID_MIM2GENE    = 59;
my $SOURCE_ID_OMIM_GENE   = 61;
my $SOURCE_ID_OMIM_MORBID = 62;
my $SPECIES_ID_HUMAN      = 9606;

# Increase this by 1 once the parser has been updated to consider synonyms
my $NUMBER_OF_DEPENDENT_LINKS      = 3;
# Decrease this by 1 if/when BaseAdaptor::get_valid_codes() has begun to
# consider synonyms
my $NUMBER_OF_STILL_UNMAPPED_XREFS = 5;

my $db = Xref::Test::TestDB->new();

# Add some dummy OMIM and EntrezGene xrefs with synonyms. Note that we
# should *not* add corresponding source entries at this point yet, as
# there is a test below which checks what happens if the parser cannot
# retrieve appropriate source IDs.

$db->schema->populate( 'Xref', [
  [ qw{ xref_id accession source_id species_id info_type } ],
  [  1, '100050', $SOURCE_ID_OMIM_MORBID, $SPECIES_ID_HUMAN, 'UNMAPPED' ], # unmapped
  [  2, '100640', $SOURCE_ID_OMIM_GENE,   $SPECIES_ID_HUMAN, 'UNMAPPED' ], # dependent
  [  3, '100100', $SOURCE_ID_OMIM_MORBID, $SPECIES_ID_HUMAN, 'UNMAPPED' ], # dependent
  [  4, '142830', $SOURCE_ID_OMIM_MORBID, $SPECIES_ID_HUMAN, 'UNMAPPED' ], # unmapped
  [  5, '142830', $SOURCE_ID_OMIM_GENE,   $SPECIES_ID_HUMAN, 'UNMAPPED' ], # unmapped
  [  6, '100660', $SOURCE_ID_OMIM_GENE,   $SPECIES_ID_HUMAN, 'UNMAPPED' ], # dependent
  [  7, '100300', $SOURCE_ID_OMIM_MORBID, $SPECIES_ID_HUMAN, 'UNMAPPED' ], # via synonym
  [  8, '999999', $SOURCE_ID_OMIM_GENE,   $SPECIES_ID_HUMAN, 'UNMAPPED' ], # not referenced

  [  9, '216',    $SOURCE_ID_ENTREZGENE,  $SPECIES_ID_HUMAN, 'DIRECT'   ], # <- 100640
  [ 10, '1131', $SOURCE_ID_ENTREZGENE,  $SPECIES_ID_HUMAN, 'DIRECT'   ], # <- 100100 
  [ 11, '218', $SOURCE_ID_ENTREZGENE,  $SPECIES_ID_HUMAN, 'DIRECT'   ], # 100660
  [ 12, '222222', $SOURCE_ID_ENTREZGENE,  $SPECIES_ID_HUMAN, 'DIRECT'   ], # not referenced <- via synonym
] );
$db->schema->populate( 'Synonym', [
  [ qw{ xref_id synonym } ],
  [  3, '121212' ],  # not referenced
  [  7, '57514' ],  # 100300
  [ 11, '3101'   ],  # not referenced
  [ 12, '3106'   ],  # not references
  [ 99, '100680' ],  # moved/removed (should be ignored; use nonexistent xref_id to provoke an error if it is)
] );

my $config = $db->config();


my $parser;

subtest 'Missing required source IDs' => sub {

  $parser = XrefParser::Mim2GeneParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_MIM2GENE,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/mim2gene-mini.txt" ],
  }); },
             qr{ \A No[ ]source_id[ ]for[ ]source_name=' }msx,
             'Throws on source IDs missing from DB' );

};


# We will need these later
$db->schema->populate( 'Source', [
  [ qw{ source_id name } ],
  [ $SOURCE_ID_ENTREZGENE,  'EntrezGene' ],
  [ $SOURCE_ID_OMIM_GENE,   'MIM_GENE' ],
  [ $SOURCE_ID_OMIM_MORBID, 'MIM_MORBID' ],
] );


subtest 'Problems with input file' => sub {

  $parser = XrefParser::Mim2GeneParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_MIM2GENE,
    species_id => $SPECIES_ID_HUMAN,
    files      => [],
  }); },
             qr{ \A No[ ]file[ ]name[ ] }msx,
             'Throws on no file name' );

  $parser = XrefParser::Mim2GeneParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_MIM2GENE,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/mim2gene-NONEXISTENT.txt" ],
  }); },
             qr{\A Could[ ]not[ ]find[ ]either[ ]' }msx,
             'Throws on missing input file' );

};

subtest 'Malformed header' => sub {
  my $QR_MALFORMED_HEADER = qr{ \A Malformed[ ]or[ ]unexpected[ ]header[ ] }msx;

  $parser = XrefParser::Mim2GeneParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_MIM2GENE,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/mim2gene-badHeader-wrongName1.txt" ],
  }); },
             $QR_MALFORMED_HEADER,
             'Throws on wrong name of first column' );

  $parser = XrefParser::Mim2GeneParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_MIM2GENE,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/mim2gene-badHeader-wrongName2.txt" ],
  }); },
             $QR_MALFORMED_HEADER,
             'Throws on wrong name of second column' );

  $parser = XrefParser::Mim2GeneParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_MIM2GENE,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/mim2gene-badHeader-wrongName3.txt" ],
  }); },
             $QR_MALFORMED_HEADER,
             'Throws on wrong name of third column' );

  $parser = XrefParser::Mim2GeneParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_MIM2GENE,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/mim2gene-badHeader-wrongName4.txt" ],
  }); },
             $QR_MALFORMED_HEADER,
             'Throws on wrong name of fourth column' );

  $parser = XrefParser::Mim2GeneParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_MIM2GENE,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/mim2gene-badHeader-wrongName5.txt" ],
  }); },
             $QR_MALFORMED_HEADER,
             'Throws on wrong name of fifth column' );

  $parser = XrefParser::Mim2GeneParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_MIM2GENE,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/mim2gene-badHeader-wrongName6.txt" ],
  }); },
             $QR_MALFORMED_HEADER,
             'Throws on wrong name of sixth column' );

};

subtest 'Malformed data' => sub {
  my $QR_WRONG_NUMBER_OF_COLUMNS
    = qr{ [ ]has[ ]an[ ]incorrect[ ]number[ ]of[ ]columns[ ] }msx;

  $parser = XrefParser::Mim2GeneParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_MIM2GENE,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/mim2gene-badData-tooFewCols.txt" ],
  }); },
             $QR_WRONG_NUMBER_OF_COLUMNS,
             'Throws on too few data columns' );

  $parser = XrefParser::Mim2GeneParser->new($db->dbh);
  throws_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_MIM2GENE,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/mim2gene-badData-tooManyCols.txt" ],
  }); },
             $QR_WRONG_NUMBER_OF_COLUMNS,
             'Throws on too many data columns' );

# add a test for valid value in type column if a xref is dependent on it
};

subtest 'Successful inserts' => sub {

  $parser = XrefParser::Mim2GeneParser->new($db->dbh);
  isa_ok( $parser, 'XrefParser::Mim2GeneParser',
          'Instantiated Mim2Gene parser' );
  lives_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_MIM2GENE,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/mim2gene-mini.txt" ],
  }); }, 'Parsed Mim2Gene map without errors' );

};

subtest 'Dependent-xref links' => sub {
  my $rs;
  my $matching_xref;
  my $matching_link;

  is( $db->schema->resultset('DependentXref')->count,
      $NUMBER_OF_DEPENDENT_LINKS,
      'Expected number of dependent-xref links' );

  $rs = $db->schema->resultset('Xref')->search({
    accession => '100100',
    source_id => $SOURCE_ID_OMIM_MORBID,
  });
  $matching_xref = $rs->next;
  is( $matching_xref->info_type, 'DEPENDENT',
      'info_type of indirectly linked xref updated correctly' );

  # Now let's have a closer look at the link
  $rs = $db->schema->resultset('DependentXref')->search({
    dependent_xref_id => $matching_xref->xref_id,
  });
  $matching_link = $rs->next;
  isnt( $matching_link, undef,
        'Indirectly linked xref has a dependent-xref link' );
  $rs = $db->schema->resultset('Xref')->search({
    xref_id => $matching_link->master_xref_id,
  });
  $matching_xref = $rs->next;
  is( $matching_xref->accession, '1131',
      'Dependent-xref link points to correct master accession' );
  is( $matching_xref->source_id, $SOURCE_ID_ENTREZGENE,
      'Master xref is an EntrezGene entry' );
  is( $matching_link->linkage_source_id, $SOURCE_ID_OMIM_MORBID,
      "Link's linkage_source_id is the dependent source_id" );
  is( $matching_link->linkage_annotation, $matching_xref->source_id,
      "Link's linkage_annotation is the master source_id" );

  # This may or may not change in the future
  $rs = $db->schema->resultset('Synonym')->search({
    synonym => '3106',
  });
  $matching_xref = $rs->next;
  is( $rs = $db->schema->resultset('DependentXref')->count({
        master_xref_id => $matching_xref->xref_id,
      }), 0, 'No dependent-xref link assigned to an EntrezGene synonym' );

};

subtest 'Non-mapping entries' => sub {
  my $rs;
  my $matching_xref;

  is( $db->schema->resultset('Xref')->count({
        info_type => 'UNMAPPED',
      }), $NUMBER_OF_STILL_UNMAPPED_XREFS,
      'Expected number of still-unmapped xrefs' );

  $rs = $db->schema->resultset('Xref')->search({
    accession => '100050',
    source_id => $SOURCE_ID_OMIM_MORBID,
  });
  $matching_xref = $rs->next;
  is( $matching_xref->info_type, 'UNMAPPED',
      'info_type of xref with Mim2Gene entry without links left alone' );

};

subtest 'Replay safety' => sub {

  $parser = XrefParser::Mim2GeneParser->new($db->dbh);
  lives_ok( sub { $parser->run({
    source_id  => $SOURCE_ID_MIM2GENE,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/mim2gene-mini.txt" ],
  }); }, 'Re-parsed Mim2Gene map without errors' );

  is( $db->schema->resultset('DependentXref')->count,
      $NUMBER_OF_DEPENDENT_LINKS,
      'No new dependent-xref links inserted by the replay' );

  is( $db->schema->resultset('Xref')->count({
        info_type => 'UNMAPPED',
      }), $NUMBER_OF_STILL_UNMAPPED_XREFS,
      'Number of still-unmapped xrefs unchanged by the replay' );

  # Ideally we would also make sure the replay has not modified
  # existing entries, no quick way of doing so though.

};


done_testing();


1;
