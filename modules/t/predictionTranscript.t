use lib 't';

BEGIN { $| = 1;  
	use Test;
	plan tests => 4;
}

my $loaded = 0;
END {print "not ok 1\n" unless $loaded;}

use MultiTestDB;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;

$loaded = 1;
my $verbose = 0;

ok(1);

my $multi = MultiTestDB->new();

my $db = $multi->get_DBAdaptor( "core" );


ok( $db );

my $sa = $db->get_SliceAdaptor();

my $slice = $sa->fetch_by_chr_start_end("20", 30_252_000, 31_252_001 );

my $p_transs = $slice->get_all_PredictionTranscripts();

debug( "Have ".scalar( @$p_transs )." prediction transcripts." );

ok( scalar( @$p_transs ) ==  27 );

for my $t ( @{$p_transs} ) {
  debug( "Transcript start ".$t->start()." end ". $t->end );
  debug( "Translation ".$t->translate()->seq() );

  my $exons = $t->get_all_Exons();
  for my $e ( @{$exons} ) {
    if( ! defined $e ) {
      debug(  "Undefined exon " );
      next;
    } 
    debug( "     Exons start ".$e->start()." end ".
	   $e->end()." strand ".$e->strand());
  }
}

debug( "All ".scalar( @$p_transs )." prediction transcripts exercised" );

ok(1);




sub debug {
  my $txt = shift;
  if( $verbose ) {
    print STDERR $txt,"\n";
  }
}
