use lib 't';
use strict;

BEGIN { $| = 1;  
	use Test ;
	plan tests => 9;
}

my $loaded = 0;
END {print "not ok 1\n" unless $loaded;}

use MultiTestDB;

my $verbose = 1;

$loaded = 1;
my $multi = MultiTestDB->new();

ok(1);

my $db = $multi->get_DBAdaptor( 'core' );

ok($db);


# Exon specific tests

my $exonad = $db->get_ExonAdaptor();
my $rca = $db->get_RawContigAdaptor();

my $contig = $rca->fetch_by_dbID( 469270 );
ok($exonad);

my $exon = Bio::EnsEMBL::Exon->new();



$exon->start(10);
ok($exon->start == 10);

$exon->end(20);
ok($exon->end == 20);

$exon->strand(1);
ok($exon->strand == 1);

$exon->phase(0);
ok($exon->phase == 0);

$exon->contig( $contig );
# should try to store (!)
$exon->end_phase( -1 );

$multi->hide( "core", "exon" );

$exonad->store($exon);

ok($exon->dbID() == 1);

# now test fetch_by_dbID

my $newexon = $exonad->fetch_by_dbID($exon->dbID);

ok($newexon);

$multi->restore();


sub debug {
  my $txt = shift;
  if( $verbose ) {
    print STDERR $txt,"\n";
  }
}


