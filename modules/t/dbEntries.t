use lib 't';
use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 31;
}

use MultiTestDB;
use TestUtils qw ( debug test_getter_setter );

use Bio::EnsEMBL::DBEntry;

# switch on the debug prints

our $verbose = 0;

debug( "Startup test" );
#
# 1 Test started
#
ok(1);

my $multi = MultiTestDB->new();

my $db = $multi->get_DBAdaptor( "core" );

debug( "Test database instatiated" );

#
# 2 Database instatiated
#
ok( $db );


# some retrievals
my $dbEntryAdaptor = $db->get_DBEntryAdaptor();


my $sth = $db->prepare( 'select count(*) from object_xref where ensembl_object_type = "Translation"' );
$sth->execute();

my ( $xref_count )  = $sth->fetchrow_array();
my $db_entry_count = 0;
my $goxref_count = 0;
my $ident_count = 0;

$sth->finish();

my $ga = $db->get_GeneAdaptor();

my $all_gene_ids = $ga->list_dbIDs();
for my $gene_id ( @$all_gene_ids ) {
  my $gene = $ga->fetch_by_dbID( $gene_id );

  for my $tr ( @{$gene->get_all_Transcripts()} ) {
    my $tl = $tr->translation();
    my $dbentries = $dbEntryAdaptor->fetch_all_by_Translation( $tl );
    $db_entry_count += scalar( @{$dbentries});
    $goxref_count += grep { $_->isa( "Bio::EnsEMBL::GoXref" )} @$dbentries;
    $ident_count += grep {$_->isa( "Bio::EnsEMBL::IdentityXref" )} @$dbentries;
  }
}

debug( "Found $xref_count xrefs and $db_entry_count dblinks." );
debug( " $goxref_count GoXrefs, $ident_count identityXrefs." );

#
# 3 as many dblinks as entries in object_xref
#
ok( $db_entry_count == $xref_count );

#
# 4,5 correct number of GoXrefs and IdentityXrefs
#
ok( $goxref_count == 48 );
ok( $ident_count == 32 );


# try storing and retrieval

my $xref = Bio::EnsEMBL::DBEntry->new
  (
   -primary_id => "1",
   -dbname => "SWISSPROT",
   -release => "1",
   -display_id => "Ens related thing"
   );


my %goxref = %$xref;
my %identxref = %$xref;

my $goref = Bio::EnsEMBL::GoXref->new
  (
   -primary_id => "1",
   -dbname => "GO",
   -release => "1",
   -display_id => "Ens related GO"
   );
$goref->add_linkage_type( "IC" );

my $ident_xref = Bio::EnsEMBL::IdentityXref->new
  (
   -primary_id => "1",
   -dbname => "SPTREMBL",
   -release => "1",
   -display_id => "Ens related Ident"
   );

$ident_xref->query_identity( 100 );
$ident_xref->target_identity( 95 );


$multi->hide( "core", "object_xref", "xref", "identity_xref", "go_xref" );


my $gene = $ga->fetch_by_dbID( $all_gene_ids->[0] );
my $tr = $gene->get_all_Transcripts()->[0];
my $tl = $tr->translation();



$dbEntryAdaptor->store( $xref, $gene, "Gene" );
$dbEntryAdaptor->store( $xref, $tr, "Transcript" );
$dbEntryAdaptor->store( $goref, $tl, "Translation" );
$dbEntryAdaptor->store( $ident_xref, $tl, "Translation" );
$dbEntryAdaptor->store( $ident_xref, $tr, "Transcript" );

my ( $oxr_count, $go_count );

$sth = $db->prepare( "select count(*) from object_xref" );
$sth->execute();

( $oxr_count ) = $sth->fetchrow_array();
$sth->finish();

#
# 6 right number of object xrefs in db
#
debug( "object_xref_count = $oxr_count" );
ok( $oxr_count == 5 );


$sth = $db->prepare( "select count(*) from xref" );
$sth->execute();

( $xref_count ) = $sth->fetchrow_array();
$sth->finish();

#
# 7 number of xrefs right
#
debug( "Number of xrefs = $xref_count" );
ok( $xref_count == 3 );

$sth = $db->prepare( "select count(*) from go_xref" );
$sth->execute();

( $go_count ) = $sth->fetchrow_array();
$sth->finish();

#
# 8 number of go entries right
#
debug( "Number of go_xrefs = $go_count" );
ok( $go_count == 1 );

$sth = $db->prepare( "select count(*) from identity_xref" );
$sth->execute();

( $ident_count ) = $sth->fetchrow_array();
$sth->finish();


#
# 9 identity xrefs right
#

# the identity (query/target)values are not normalized ...
debug( "Number of identity_xrefs = $ident_count" );
ok( $ident_count == 2 );


$multi->restore();


#
# 10-12 Test that external synonyms and go evidence tags are retrieved
#
my $ta = $db->get_TranscriptAdaptor();
my $translation = $ta->fetch_by_dbID(21737)->translation;

#pull out an xref that should 2 external synonyms
my $xrefs = $dbEntryAdaptor->fetch_all_by_Translation($translation);
($xref) = grep {$_->dbID == 315} @$xrefs;
my @syns = grep {$_ eq 'syn1' || $_ eq 'syn2'} @{$xref->get_all_synonyms};
ok(@syns == 2);

#and also 2 evidence tags
if($xref && $xref->isa('Bio::EnsEMBL::GoXref')) {
  my @evtags = 
    grep {$_ eq 'IEA' || $_ eq 'IC'} @{$xref->get_all_linkage_types()};
  ok(@evtags == 2);
} else {
  ok(0);
}


$translation = $ta->fetch_by_dbID(21723)->translation;

$xrefs = $dbEntryAdaptor->fetch_all_by_Translation($translation);
($xref) = grep {$_->dbID == 257} @$xrefs;

if($xref && $xref->isa('Bio::EnsEMBL::GoXref')) {
  my ($evtag) = @{$xref->get_all_linkage_types()};
  ok($evtag eq 'IC');
} else {
  ok(0);
}


#test the idxref mapping code a bit
my $id_xref = Bio::EnsEMBL::IdentityXref->new
  (-QUERY_IDENTITY    => 80.4,
   -TARGET_IDENTITY   => 90.1,
   -SCORE             => 100,
   -EVALUE            => 1,
   -CIGAR_LINE        => '10M5D10M5I10M',
   -QUERY_START       => 1,
   -QUERY_END         => 35,
   -TRANSLATION_START => 10,
   -TRANSLATION_END   => 45,
   -ADAPTOR           => $dbEntryAdaptor);


my $mapper = $id_xref->get_mapper();

ok($mapper && ref($mapper) && $mapper->isa('Bio::EnsEMBL::Mapper'));

my @coords = $mapper->map_coordinates('external_id', 10, 100, 1, 'external');

ok(@coords == 5);

ok($coords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate'));
ok($coords[0]->start() == 19);
ok($coords[0]->end()   == 19);

ok($coords[1]->isa('Bio::EnsEMBL::Mapper::Coordinate'));
ok($coords[1]->start == 25);
ok($coords[1]->end   == 34);

ok($coords[2]->isa('Bio::EnsEMBL::Mapper::Gap'));
ok($coords[2]->length == 5);

ok($coords[3]->isa('Bio::EnsEMBL::Mapper::Coordinate'));
ok($coords[3]->start == 35);
ok($coords[3]->end == 44);

ok($coords[4]->isa('Bio::EnsEMBL::Mapper::Gap'));
ok($coords[4]->length == 65);



my $feature = Bio::EnsEMBL::Feature->new
  (-START => 10,
   -END   => 100,
   -STRAND => 1,
   -SEQNAME => 'SwisProtSeq');

@coords = $id_xref->map_feature($feature);

ok(@coords == 5);
ok($coords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate') &&
   $coords[0]->start() == 19 &&
   $coords[0]->end()   == 19);


my $features = $id_xref->transform_feature($feature);

ok(@$features == 3);

ok($features->[0]->start == 19 &&
   $features->[0]->end   == 19 &&
   $features->[1]->start == 25 &&
   $features->[1]->end   == 34 &&
   $features->[2]->start == 35 &&
   $features->[2]->end   == 44);


#
# test DBEntryAdaptor::fetch_by_db_accession and
#      DBEntryAdaptor::fetch_by_dbID
#

$xref = $dbEntryAdaptor->fetch_by_dbID(152202);

ok($xref->dbID == 152202);
ok($xref->display_id() eq 'C20orf125');
ok($xref->dbname() eq 'HUGO');
ok($xref->primary_id() eq '16118');


$xref = $dbEntryAdaptor->fetch_by_db_accession('Interpro', 'IPR000010');
ok($xref->dbID == 999999);
ok($xref->description eq 'Test interpro desc2');

$multi->hide('core', 'xref');

#make sure this still works when interpro entries missing in xref table
$xref = $dbEntryAdaptor->fetch_by_db_accession('Interpro', 'IPR000010');
ok($xref);
ok($xref->primary_id() eq 'IPR000010');

$multi->restore('core', 'xref');
