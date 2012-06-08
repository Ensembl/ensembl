use strict;
use warnings;

use Test::More;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::DBSQL::DBEntryAdaptor;
use Bio::EnsEMBL::DBEntry;

note( "Startup test" );
#
# 1 Test started
#
ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $db = $multi->get_DBAdaptor( "core" );

note( "Test database instantiated" );

#
# 2 Database instatiated
#
ok( $db );

# some retrievals
my $dbEntryAdaptor = $db->get_DBEntryAdaptor();


my $sth = $db->dbc->prepare( 'select count(*) from object_xref where ensembl_object_type = "Translation"' );
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
    $goxref_count += grep { $_->isa( "Bio::EnsEMBL::OntologyXref" )} @$dbentries;
    $ident_count += grep {$_->isa( "Bio::EnsEMBL::IdentityXref" )} @$dbentries;
  }
}

note( "Found $xref_count xrefs and $db_entry_count dblinks." );
note( " $goxref_count GoXrefs, $ident_count identityXrefs." );

#
# 3 as many dblinks as entries in object_xref
#
note $db_entry_count."\t".$xref_count."\n";
ok( $db_entry_count == $xref_count );

#
# 4,5 correct number of GoXrefs and IdentityXrefs
#
note( "GoXrefs and IdentityXrefs: ".$goxref_count . " " . $ident_count);
ok( $goxref_count == 48 );
ok( $ident_count == 32 );

my $analysis_adap = $db->get_AnalysisAdaptor();
my $analysis = $analysis_adap->fetch_by_logic_name("RepeatMask");

# try storing and retrieval

my $xref = Bio::EnsEMBL::DBEntry->new
  (
   -primary_id => "1",
   -dbname => "Uniprot/SWISSPROT",
   -release => "1",
   -display_id => "Ens related thing",
   -primary_id_linkable => "0",
   -display_id_linkable => "1",
   -priority => "5",
   -db_display_name => "Nice friendly name",
   -info_type => "PROJECTION",
   -info_text => "from human gene ENSG0000011111",
   -type => "ARRAY",
    -analysis => $analysis
   );


my %goxref = %$xref;
my %identxref = %$xref;

my $ident_xref = Bio::EnsEMBL::IdentityXref->new
  (
   -primary_id => "1",
   -dbname => "Uniprot/SPTREMBL",
   -release => "1",
   -display_id => "Ens related Ident",
   -analysis => $analysis   );

$ident_xref->xref_identity( 100 );
$ident_xref->ensembl_identity( 95 );

my $goref = Bio::EnsEMBL::OntologyXref->new
  (
   -primary_id => "1",
   -dbname => "GO",
   -release => "1",
   -display_id => "Ens related GO",
   -analysis => $analysis
   );
$goref->add_linkage_type( "IC" ); # Linkage type on own
$goref->add_linkage_type( "ISS", $goref ); # Linkage type with source xref


$multi->hide( "core", "object_xref", "xref", "identity_xref", "ontology_xref" );


my $gene = $ga->fetch_by_dbID( $all_gene_ids->[0] );
my $tr = $gene->get_all_Transcripts()->[0];
my $tl = $tr->translation();

my $oxr_count = count_rows($db, 'object_xref');
$dbEntryAdaptor->store( $xref, $tr->dbID, "Transcript" );
$oxr_count = count_rows($db, 'object_xref');
$dbEntryAdaptor->store( $ident_xref, $tl->dbID, "Translation" );
#intentional duplicates should be filtered out and not increase row count
$dbEntryAdaptor->store( $ident_xref, $tl->dbID, "Translation" ); 
$dbEntryAdaptor->store( $ident_xref, $tl->dbID, "Translation" );
$dbEntryAdaptor->store( $ident_xref, $tl->dbID, "Translation" );
$dbEntryAdaptor->store( $ident_xref, $tl->dbID, "Translation" );


$oxr_count = count_rows($db, 'object_xref');
$dbEntryAdaptor->store( $goref, $tl->dbID, "Translation" );
$oxr_count = count_rows($db, 'object_xref');
$dbEntryAdaptor->store( $ident_xref, $tr->dbID, "Transcript" );
my ($go_count );

$oxr_count = count_rows($db, 'object_xref');
#
# 6 right number of object xrefs in db
#
note( "object_xref_count = $oxr_count" );
ok( $oxr_count == 4 );


$xref_count = count_rows($db, 'xref');
$sth->finish();

#
# 7 number of xrefs right
#
note( "Number of xrefs = $xref_count" );
ok( $xref_count == 3 );

#
# 8 number of go entries right
#
$go_count = count_rows($db, 'ontology_xref');
note( "Number of go_xrefs = $go_count" );
ok( $go_count == 2 );

#
# 9 identity xrefs right
#

$ident_count = count_rows($db, 'identity_xref');
# the identity (query/target)values are not normalized ...
note( "Number of identity_xrefs = $ident_count" );
ok( $ident_count == 2 );

# Check type storing and retrieval
ok($xref->type() eq 'ARRAY');

$multi->restore();

# test parallel insertions of identical xrefs

$xref = Bio::EnsEMBL::DBEntry->new
  (
   -primary_id => "1",
   -dbname => "Vega_gene",
   -release => "1",
   -display_id => "Ens fake thing",
   -primary_id_linkable => "0",
   -display_id_linkable => "1",
   -priority => "5",
   -db_display_name => "Nice unfriendly name",
   -info_type => "MISC",
   -info_text => "Concurrent insert",
   -type => "ARRAY",
    -analysis => $analysis
   );
# db connection must be severed for threads to access DB    
$dbEntryAdaptor->dbc->disconnect_if_idle();
use threads;

sub parallel_store {
    my $xref_id = $dbEntryAdaptor->store( $xref, $tr->dbID, "Transcript" );
    return $xref_id
}
   
my $thread1 = threads->create(\&parallel_store);
my $thread2 = threads->create(\&parallel_store);
my $thread3 = threads->create(\&parallel_store);

    
my @xref_ids;
@xref_ids = ($thread1->join,$thread2->join,$thread3->join);

note("Threaded xrefs: ".$xref_ids[0]." ".$xref_ids[1]." ".$xref_ids[2]);

# Test 10 - Verify that only one xref has been inserted under parallel inserts
ok($xref_ids[0] == 1000009 && $xref_ids[1] == $xref_ids[0] && $xref_ids[2] == $xref_ids[0]);

# Test 11 - Exception testing on ->store()

$xref = Bio::EnsEMBL::DBEntry->new
  (
   -primary_id => "1",
   -dbname => "Vega_gene",
   -release => "1",
   -display_id => "Ens fakiest thing",
   -primary_id_linkable => "0",
   -display_id_linkable => "1",
   -priority => "5",
   -db_display_name => "Nice unfriendly name",
   -info_type => "MISC",
   -info_text => "Full exception checking",
   -type => "ARRAY",
   -analysis => undef,
   );
   
my $xref_id = $dbEntryAdaptor->store($xref, undef, "Transcript");
note("Xref_id from insert: ".$xref_id);
ok($xref_id == 1000010);

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

#and also 2 evidence tags, and one source_xref
if($xref && $xref->isa('Bio::EnsEMBL::OntologyXref')) {
  my @evtags = 
    grep {$_ eq 'IEA' || $_ eq 'IC'} @{$xref->get_all_linkage_types()};
  ok(@evtags == 2);
  my @source_xrefs = 
    grep {UNIVERSAL::isa($_->[1],'Bio::EnsEMBL::DBEntry')}
	@{$xref->get_all_linkage_info};
  ok(@source_xrefs == 1);
} else {
  ok(0);
}


$translation = $ta->fetch_by_dbID(21723)->translation;

$xrefs = $dbEntryAdaptor->fetch_all_by_Translation($translation);
($xref) = grep {$_->dbID == 257} @$xrefs;

if($xref && $xref->isa('Bio::EnsEMBL::OntologyXref')) {
  my ($evtag) = @{$xref->get_all_linkage_types()};
  ok($evtag eq 'IC');
} else {
  ok(0);
}


#test the idxref mapping code a bit
my $id_xref = Bio::EnsEMBL::IdentityXref->new
  (-XREF_IDENTITY    => 80.4,
   -ENSEMBL_IDENTITY   => 90.1,
   -SCORE             => 100,
   -EVALUE            => 1,
   -CIGAR_LINE        => '10M5D10M5I10M',
   -XREF_START       => 1,
   -XREF_END         => 35,
   -ENSEMBL_START => 10,
   -ENSEMBL_END   => 45,
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

$multi->save('core', 'object_xref', 'identity_xref', 'ontology_xref');

#
# test the removal of dbentry associations
#

$translation = $ta->fetch_by_dbID(21723)->translation;

my $dbes = $translation->get_all_DBEntries();

my $all_count = @$dbes;
$go_count  = grep {$_->isa('Bio::EnsEMBL::OntologyXref')} @$dbes;
my $id_count  = grep {$_->isa('Bio::EnsEMBL::IdentityXref')} @$dbes;

my $all_total = count_rows($db, 'object_xref');
my $go_total  = count_rows($db, 'ontology_xref');
my $id_total  = count_rows($db, 'identity_xref');

print_dbEntries($dbes);


foreach my $dbe (@$dbes) {
  $dbEntryAdaptor->remove_from_object($dbe, $translation, 'Translation');
}

# make sure the appropriate rows were deleted

ok($all_total - $all_count == count_rows($db, 'object_xref'));
ok($go_total - $go_count   == count_rows($db, 'ontology_xref'));
ok($id_total - $id_count   == count_rows($db, 'identity_xref'));

$multi->restore('core', 'object_xref', 'identity_xref', 'ontology_xref');


# new type checks

# fetch_all_by Gene Transcript Translation

# translation info type

$translation = $ta->fetch_by_dbID(21737)->translation;
$dbes = $translation->get_all_DBEntries(undef,"MISC");

ok(@$dbes == 19); # test 44

$dbes = $translation->get_all_DBLinks(undef,"MISC");
ok(@$dbes == 19);

$dbes = $translation->get_all_DBLinks(undef,"ALT_TRANS");
ok(@$dbes == 9);

$dbes = $translation->get_all_DBLinks(undef,"ARRAY");
ok(@$dbes == 1);

$dbes = $translation->get_all_DBLinks(undef,"RUBBISH");
ok(@$dbes == 0);



# transcript external type

my $transcript = $ta->fetch_by_dbID(21737);
$dbes = $transcript->get_all_DBEntries(undef,"MISC");

ok(@$dbes == 0); # test 49

$dbes = $transcript->get_all_DBLinks(undef,"MISC");
ok(@$dbes == 19);

$dbes = $transcript->get_all_DBLinks(undef,"ALT_TRANS");
ok(@$dbes == 9);

$dbes = $transcript->get_all_DBLinks(undef,"ARRAY");
ok(@$dbes == 1);

$dbes = $transcript->get_all_DBLinks(undef,"RUBBISH");
ok(@$dbes == 0);



# gene external type

$gene = $ga->fetch_by_dbID(18272);
$dbes = $gene->get_all_DBEntries(undef,"MISC");

ok(@$dbes == 0); # test 54

$dbes = $gene->get_all_DBLinks(undef,"MISC");
ok(@$dbes == 19);

$dbes = $gene->get_all_DBLinks(undef,"ALT_TRANS");
ok(@$dbes == 9);

$dbes = $gene->get_all_DBLinks(undef,"ARRAY");
ok(@$dbes == 1);

$dbes = $gene->get_all_DBLinks(undef,"RUBBISH");
ok(@$dbes == 0);

# test fetch_all_by_source

$xrefs = $dbEntryAdaptor->fetch_all_by_source("%Uniprot%");
ok(@{$xrefs} == 23);  #test 60

my $db_name = $dbEntryAdaptor->get_db_name_from_external_db_id(4100);
ok($db_name eq 'UniGene');

# Test multiple inserts for empty descriptions
{
  $multi->hide('core', 'xref', 'object_xref');
  my @basic_args = (-PRIMARY_ID => 'AAAA', -DBNAME => 'Uniprot/SWISSPROT', -RELEASE => 1);
  my $entry_no_desc = Bio::EnsEMBL::DBEntry->new(@basic_args, -DESCRIPTION => q{});
  my $no_desc_id = $dbEntryAdaptor->store($entry_no_desc, $gene->dbID(), 'Gene');
  is_rows(1, $db, 'xref', 'where description = ?', [q{}]);
  is_rows(1, $db, 'object_xref');
  my $no_desc_id_again = $dbEntryAdaptor->store($entry_no_desc, $gene->dbID(), 'Gene');
  is($no_desc_id_again, $no_desc_id, 'Checking the ID is consistent between store() invocations');
  is_rows(1, $db, 'xref', 'where description = ?', [q{}]);
  is_rows(1, $db, 'object_xref');
  is_rows(0, $db, 'object_xref', 'where xref_id =?', [0]);
  
  $multi->restore('core', 'xref', 'object_xref');
}

sub print_dbEntries {
  my $dbes = shift;

  foreach my $dbe (@$dbes) {
    if($dbe->isa('Bio::EnsEMBL::IdentityXref')) {
      note("IDXref");
    } elsif($dbe->isa('Bio::EnsEMBL::OntologyXref')) {
      note("GOXref");
    } elsif($dbe->isa('Bio::EnsEMBL::DBEntry')) {
      note("DBEntry");
    } else {
      note("UNKNOWN dbentry type");
    }

    note(" ".$dbe->dbname()."-".$dbe->display_id()."\n");
  }

  note(scalar(@$dbes). " total");

}

done_testing();