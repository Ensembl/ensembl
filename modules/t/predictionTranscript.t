use lib 't';

BEGIN { $| = 1;  
	use Test;
	plan tests => 16;
}

my $loaded = 0;
END {print "not ok 1\n" unless $loaded;}

use MultiTestDB;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;

$loaded = 1;

ok(1);

my $multi = MultiTestDB->new();

my $db = $multi->get_DBAdaptor( "core" );


ok( $db );

my $sa = $db->get_SliceAdaptor();

my $slice = $sa->fetch_by_chr_start_end("20", 30_252_000, 31_252_001 );

my $p_transs = $slice->get_all_PredictionTranscripts();

# print STDERR "Number of Transcripts: ",

ok( scalar( @$p_transs ) ==  25 );



print STDERR "Exon count ", scalar( @{$exons} ),"\n";

for my $t ( @{$p_transs} ) {
  my $exons = $t->get_all_Exons();
    for my $e ( @{$exons} ) {
      if( ! defined $e ) {
	print "Undefined exon\n";
        next;
      } 
      print "     Exons start ", $e->start(), " end ", $e->end(), " strand ", $e->strand(),"\n";
  }
    print "start ",$t->start(), " end ", $t->end,"\n";
  print STDERR "Translation ",$t->translate()->seq(),"\n";
}




ok(1);
