use strict;

BEGIN { $| = 1;  
	use Test;
	plan tests => 17;
}

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Map::DitagFeature;
use Bio::EnsEMBL::Analysis;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor( 'core' );

my $dfa   = $db->get_DitagFeatureAdaptor;

my $slice         = $db->get_SliceAdaptor->fetch_by_region('chromosome','20');
my $analysis      = $db->get_AnalysisAdaptor->fetch_by_logic_name('DitagAlign' );

my $dbID          = 4828567;
my $ditag_id      = 3278337;
my $qstart        = 120635196;
my $qend          = 120635214;
my $qstrand       = 1;
my $tstart        = 1;
my $tend          = 19;
my $tstrand       = 1;
my $ditag_side    = 'L';
my $ditag_pair_id = 1;
my $cigar_line    = '19M';


######
# 1  #
######

#test new

my $feature = Bio::EnsEMBL::Map::DitagFeature->new(
					-slice         => $slice,
					-start         => $qstart,
					-end           => $qend,
					-strand        => $qstrand,
					-hit_start     => $tstart,
					-hit_end       => $tend,
					-hit_strand    => $tstrand,
					-ditag_id      => $ditag_id,
					-ditag_side    => $ditag_side,
					-ditag_pair_id => $ditag_pair_id,
					-cigar_line    => $cigar_line,
					-analysis      => $analysis,
 					);
ok($feature && $feature->isa('Bio::EnsEMBL::Map::DitagFeature'));

########
# 2-14 #
########

#test ditag_id, ditag_side, hit_start, hit_end,
# hit_strand, cigar_line, start, end, strand,
# dbID, sequence, slice, ditag_pair_id

my $ditagFeatures = $dfa->fetch_by_ditagID($ditag_id);
my $ditagFeature  = $ditagFeatures->[0];

ok(defined $ditagFeature && $ditagFeature->isa('Bio::EnsEMBL::Map::DitagFeature'));

ok($ditagFeature->ditag_id      == $ditag_id);
ok($ditagFeature->slice->isa('Bio::EnsEMBL::Slice'));
ok($ditagFeature->ditag_pair_id == $ditag_pair_id);
ok($ditagFeature->ditag_side    eq $ditag_side);
ok($ditagFeature->hit_start     == $tstart);
ok($ditagFeature->hit_end       == $tend);
ok($ditagFeature->hit_strand    eq $tstrand);
ok($ditagFeature->cigar_line    eq $cigar_line);
ok($ditagFeature->start         == $qstart);
ok($ditagFeature->end           == $qend);
ok($ditagFeature->strand        eq $qstrand);
ok($ditagFeature->dbID          == $dbID);
ok(length($ditagFeature->sequence) > 10);

######
# 13 #
######

#test fetch_ditag

my $ditag = $ditagFeature->fetch_ditag();

ok(defined $ditag && $ditag->isa('Bio::EnsEMBL::Map::Ditag'));
ok($ditag->dbID == $ditag_id);

1;
