use lib 't';

BEGIN { $| = 1;  
	use Test;
	plan tests => 19;
}

use MultiTestDB;
use TestUtils qw(debug);

use Bio::EnsEMBL::Utils::Exception qw(verbose);

######################################################################
# 
# Clone is a deprecated class but needed for backwards 
# compatibility.  These tests ensure that it actually works,
# but verbosity is turned off to avoid all of the deprecated warnings
#
#######################################################################

verbose(-1);

our $verbose = 0; #set to 1 for debug printing

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

my $base_count = $contig->get_base_count();
debug( "BaseCount: ". 
       join(' ', map({"$_ => $base_count->{$_}\n"} keys(%$base_count))));
my $a = $base_count->{'a'};
my $c = $base_count->{'c'};
my $t = $base_count->{'t'};
my $g = $base_count->{'g'};
my $n = $base_count->{'n'};
my $gc_content = $base_count->{'%gc'};
ok( $a == 6395 &&
    $c == 6070 &&
    $t == 6539 &&
    $g == 6006 &&
    $n == 0 &&
    $gc_content == 48.28 && 
    $a + $t + $c + $g + $n == $contig->length());



verbose(0);