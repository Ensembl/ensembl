use lib 't';
use strict;
use warnings;

BEGIN { $| = 1;  
	use Test;
	plan tests => 14;
}

use MultiTestDB;
use TestUtils qw ( debug test_getter_setter );

use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::ProteinFeature;

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

my $analysis = Bio::EnsEMBL::Analysis->new
  ( -logic_name => 'protmap',
    -program => 'pmatch',
    -database => 'pmatch' );


my $ident_xref = Bio::EnsEMBL::IdentityXref->new
  (
   -primary_id => "1",
   -dbname => "SPTREMBL",
   -release => "1",
   -display_id => "Ens related Ident",
   -cigar_line => "10MDI2M3D2M3I5M",
   -translation_start => 1,
   -translation_end => 23,
   -hit_start => 1,
   -hit_end => 23,
   -score => 20.0,
   -evalue => 123,
   -analysis => $analysis
   );

$ident_xref->query_identity( 100 );
$ident_xref->target_identity( 95 );


$multi->hide( "core", "object_xref", "xref", "identity_xref", "go_xref", "analysis" );


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


my $protein_feature = Bio::EnsEMBL::ProteinFeature->new
  (
   -start => 1,
   -end => 23,
   -seqname => 'SP0001'
  );

my $mapper = $ident_xref->get_mapper();

debug( "Mapper ".ref( $mapper ) );

my $copy_features = $ident_xref->transform_feature( $protein_feature );

for my $feat ( @$copy_features ) {
  debug(join("\n", map({$_ ."->" . $feat->{$_} } keys %$feat) ));
  debug( "------------");
}

# 4 M-sections in the xref make 4 mapped sections

#
# 10 mapping features via alignment
#
ok( scalar( @$copy_features ) == 4 );

# 10M is the first matching section

#
# 11 mapping first section to 1-10
#
ok( $copy_features->[0]->start() == 1 &&
    $copy_features->[0]->end() == 10 );

$multi->restore();


#
# 12-14 Test that external synonyms and go evidence tags are retrieved
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


#
# Need more tests here ...
#



 
