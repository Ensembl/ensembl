use strict;

BEGIN { $| = 1;
	use Test ;
	plan tests => 14;
}

my $loaded = 0;
END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0; #set to 1 to turn on debug printouts

$loaded = 1;
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

ok(1);

my $db = $multi->get_DBAdaptor( 'core' );

ok($db);

my $stable_id = 'ENST00000217347';
my $transcript_adaptor = $db->get_TranscriptAdaptor();
my $transcript = 
  $transcript_adaptor->fetch_by_stable_id($stable_id);


my @exons = (@{$transcript->get_all_Exons()});  
my @introns = (@{$transcript->get_all_Introns()});  

my $i=0;
foreach my $intron (@introns){
  ok($intron->prev_Exon->end == $intron->start-1);
  ok($intron->next_Exon->start == $intron->end+1);
}


$multi->restore();

