use lib 't';
use strict;
use warnings;
use vars qw( $verbose );

BEGIN { $| = 1;
	use Test;
	plan tests => 35;
}

use MultiTestDB;
use TestUtils qw( debug test_getter_setter );
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Slice;

my $multi = MultiTestDB->new();

$verbose = 0; #set to true to turn on debug print outs

ok( $multi );


my $db = $multi->get_DBAdaptor('core' );

my $sa = $db->get_SliceAdaptor();

my $slice = $sa->fetch_by_region('chromosome', "20", 30_249_935, 31_254_640 );

my $genes = $slice->get_all_Genes();

my $translates = 1;
my $utr_trans;

debug( "Checking if all Transcripts translate and transform" );

for my $gene ( @$genes ) {
  for my $trans ( @{$gene->get_all_Transcripts()} ) {
    if( $trans->translate()->seq() =~ /\*./ ) {
      $translates = 0;
      debug( $trans->stable_id()." does not translate." );
      last;
    }
    if( $trans->coding_region_start() != $trans->start() &&
	$trans->coding_region_end() != $trans->end() ) {
      $utr_trans = $trans->stable_id();
    }
  }

  $gene = $gene->transform( "contig" );
  next if( ! $gene );	
  for my $trans ( @{$gene->get_all_Transcripts()} ) {
    if( $trans->translate()->seq() =~ /\*./ ) {
      $translates = 0;
      debug( $trans->stable_id()." does not translate." );
      last;
    }
  }
}

debug( "utr Transcript is $utr_trans" );

ok( $translates );

#
# Try pulling off genes from an NTContig and making sure they still translate.
# This is a fairly good test of the chained coordinate mapping since the
# transcripts are stored in chromosomal coordinates and there is no direct
# mapping path to the supercontig coordinate system.
#
my $supercontig = $sa->fetch_by_region('supercontig', "NT_028392");

my $transcripts = $supercontig->get_all_Transcripts();

debug( "Checking if all transcripts on NTContig NT_028392 translate and transform" );

$translates = 1;
for my $trans ( @$transcripts ) {
  if( $trans->translate()->seq() =~ /\*./ ) {
    $translates = 0;
    debug( $trans->stable_id()." does not translate." );
    last;
  }
  debug($trans->stable_id() .  ":" . $trans->translate()->seq() . "\n");
}

ok($translates);


my $ta = $db->get_TranscriptAdaptor();

debug ("Transcript->list_dbIDs");
my $ids = $ta->list_dbIDs();
ok (@{$ids});

debug ("Transcript->list_stable_ids");
my $stable_ids = $ta->list_stable_ids();
ok (@{$stable_ids});

my $tr = $ta->fetch_by_stable_id( "ENST00000217347" );

$tr = $tr->transform('contig');

ok( $tr );

debug ( "External transcript name: " . $tr->external_name );
ok ( $tr->external_name eq "MAPRE1");

debug ( "External transcript dbname: " . $tr->external_db );
ok ( $tr->external_db eq 'HUGO' );

debug ( "Display_xref_id: " . $tr->display_xref->dbID() );
ok ( $tr->display_xref->dbID() == 97759 );
ok( test_getter_setter( $tr, "display_xref", 42 ));

ok( test_getter_setter( $tr, "dbID", 100000 ));
ok( test_getter_setter( $tr, "type", "NOVEL" ));
ok( $tr->translation->isa( "Bio::EnsEMBL::Translation" ));

debug( "start() == ".$tr->start() );
ok( $tr->start() == 79874 );

debug( "end() == ".$tr->end() );
ok( $tr->end() == 110306 );

debug( "spliced_seq->substr == \"".substr( $tr->spliced_seq(),0, 10 )."\"" );
ok( substr( $tr->spliced_seq(), 0, 10 ) eq "ACGAGACGAA" ); 

debug( "translateable_seq->substr == \"".substr( $tr->translateable_seq(),0,10 )."\"" );
ok( substr( $tr->translateable_seq(),0,10 ) eq "ATGGCAGTGA" );

debug( "coding_region_start() == ".$tr->coding_region_start() );
ok( $tr->coding_region_start() == 85834 );

debug( "coding_region_end() == ".$tr->coding_region_end() );
ok( $tr->coding_region_end() == 108631 );

debug( "pep2genomic: ".($tr->pep2genomic( 10,20 ))[0]->start());
my @pepcoords = $tr->pep2genomic( 10, 20 );
ok( $pepcoords[0]->start() == 85861 );

my $t_start = $tr->start;
my $t_end   = $tr->end;
my $t_strand = $tr->get_all_Exons->[0]->strand;
my $pep_len = ($tr->cdna_coding_end - $tr->cdna_coding_start + 1) / 3;

my $coord_num = 0;
my $coding_region_start = $tr->coding_region_start;
my $coding_region_end   = $tr->coding_region_end;
#expect coordinate for each exon with coding sequence 
foreach my $e (@{$tr->get_all_Exons}) {
  if($e->end > $coding_region_start && $e->start < $coding_region_end) {
    $coord_num++;
  }
}

my @coords = ();
my @gaps = ();
my $coord_len = 0;
foreach my $c ($tr->genomic2pep($t_start, $t_end, $t_strand)) {
  if($c->isa('Bio::EnsEMBL::Mapper::Gap')) {
    push @gaps, $c;
  } else {
    push @coords, $c;
  }
}

debug("expecting $coord_num coords (in pep coords), got " . scalar(@coords));
ok(scalar(@coords) == $coord_num);
my ($last_coord) = reverse(@coords); 
debug("expecting peptide length: $pep_len, got".$last_coord->end);
ok($pep_len == $last_coord->end);

debug( "start Exon: ".$tr->start_Exon->stable_id() );
debug( "end Exon: ".$tr->end_Exon->stable_id() );

debug( "cdna_coding_start: ". $tr->cdna_coding_start );
ok($tr->cdna_coding_start == 65);
ok(test_getter_setter($tr, 'cdna_coding_start', 99));

debug( "five_prime_utr: ".substr( $tr->five_prime_utr()->seq(), -5 , 5 ));
ok( substr( $tr->five_prime_utr()->seq(), -5, 5) eq "CGAAG" ); 

debug( "cdna_coding_end: ". $tr->cdna_coding_end );
ok($tr->cdna_coding_end == 868);
ok(test_getter_setter($tr, 'cdna_coding_end', 102));

debug( "three_prime_utr: ".substr( $tr->three_prime_utr()->seq(), -5, 5  ));
ok( substr( $tr->three_prime_utr()->seq(), -5, 5  ) eq "TTCAA");

debug( "Transcript has: ". scalar( @{$tr->get_all_Exons()} ). " Exons" );
ok( scalar( @{$tr->get_all_Exons()} ) == 7 );

debug( "Flushing Exons" );
$tr->flush_Exons();

ok( scalar( @{$tr->get_all_Exons()} ) == 0 );


# get a fresh tr to check the update method
$tr = $ta->fetch_by_stable_id( "ENST00000217347" );

$multi->save('core', 'transcript');

# the first update should have no effect
$ta->update($tr);

my $up_tr = $ta->fetch_by_stable_id( "ENST00000217347" );
ok ( $up_tr->display_xref->dbID() == 97759 );

my $dbentryAdaptor = $db->get_DBEntryAdaptor();

$tr->display_xref($dbentryAdaptor->fetch_by_dbID( 614 ));
$ta->update($tr);

$up_tr = $ta->fetch_by_stable_id( "ENST00000217347" );
ok ( $up_tr->display_xref->dbID() == 614 );

$multi->restore('core', 'transcript');



my $interpro = $ta->get_Interpro_by_transid("ENST00000252021");
foreach my $i (@$interpro) {
  debug($i);
}
###currently no interpro info in the test db
ok(@$interpro == 1);

#
# test fetch_all_by_external_name
#

($tr) = @{$ta->fetch_all_by_external_name('BAB15482')};
ok($tr && $tr->stable_id eq 'ENST00000262651');

#
# test fetch_by_translation_id
#

$tr = $ta->fetch_by_translation_id(21734);

ok($tr && $tr->stable_id eq 'ENST00000201961');
