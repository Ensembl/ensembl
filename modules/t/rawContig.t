use lib 't';

BEGIN { $| = 1;  
	use Test;
	plan tests => 8;
}

my $loaded = 0;
END {print "not ok 1\n" unless $loaded;}

use MultiTestDB;

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


# $contig->subseq( start end strand)
# $contig->get_repeatmasked_seq( logic_name soft_mask_enable_flag )
# $contig->get_all_PredictionTranscripts()
# $contig->get_all_RepeatFeatures()
# $contig->get_all_SimilarityFeatures( logicname score )
# $contig->get_all_DnaAlignFeatures( logic_name score )
# $contig->get_all_ProteinAlignFeatures( logic_name score )
# $contig->get_all_SimpleFeatures( logic_name score )
# $contig->embl_offset()
# $contig->clone()
# $contig->get_all_ExternalFeatures()

# we should be able to store RawContigs, shouldnt we?



