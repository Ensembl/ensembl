use lib 't';

BEGIN { $| = 1;  
	use Test;
	plan tests => 18;
}

my $loaded = 0;
END {print "not ok 1\n" unless $loaded;}

use MultiTestDB;

my $verbose = 0;

$loaded = 1;

ok(1);

my $multi = MultiTestDB->new();

ok($multi);


my $db = $multi->get_DBAdaptor( "core" );

$cadp = $db->get_RawContigAdaptor();
$contig = $cadp->fetch_by_dbID(317101);

ok($contig);


$contig = $cadp->fetch_by_name('AL353092.6.1.25010');

ok($contig->dbID == 339816);

ok($contig->name eq 'AL353092.6.1.25010');

$clone = $contig->clone();

ok($clone);

ok($clone->embl_id eq 'AL353092');

ok($contig->seq);

debug( "Contig length: ".$contig->length() );
ok( $contig->length() == 25010 );

my $subseq = $contig->subseq( 1, 5, -1 );
debug( "subseq: $subseq" );

ok( $subseq eq "GAACT" );


my $ptrans = $contig->get_all_PredictionTranscripts();
ok( $ptrans );

my $rFeatures = $contig->get_all_RepeatFeatures();
debug( "repeat features: ".@{$rFeatures} );

ok( @$rFeatures == 55 );

my $repseq  = $contig->get_repeatmasked_seq( );
debug( "Repeatmasked: ".substr( $repseq->seq(), 0, 50 ) );
debug( "  isa ".ref( $repseq ));

ok( $repseq->isa( "Bio::PrimarySeq" ));


my $sFeatures = $contig->get_all_SimilarityFeatures( "swall" );
debug( "SimilarityFeatures swall: ".@{$sFeatures} );
ok( @$sFeatures == 123 );

my $daFeatures = $contig->get_all_DnaAlignFeatures( );
debug( "DnaAlignFeatures ".@{$daFeatures} );
ok( @$daFeatures == 392 );


my $pFeatures = $contig->get_all_ProteinAlignFeatures( );
debug( "Protein Align Features: ".@{$pFeatures} );
ok( @$pFeatures == 136 );

my $simFeatures = $contig->get_all_SimpleFeatures( );
debug( "Simple Features: ".@{$simFeatures} );

debug( "Embloffset: ". $contig->embl_offset());
ok( $contig->embl_offset() == 1 );

my $extFeatures = $contig->get_all_ExternalFeatures();
debug( "External Features: ".@{$extFeatures} );
ok( @$extFeatures == 0 );

sub debug {
  my $txt = shift;
  if( $verbose ) {
    print STDERR $txt,"\n";
  }
}



