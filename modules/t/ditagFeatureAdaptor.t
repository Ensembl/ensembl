use strict;

BEGIN { $| = 1;  
	use Test ;
	plan tests => 12;
}

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Map::Ditag;
use Bio::EnsEMBL::Analysis;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor( 'core' );

my $region        = 11;
my $ditag_id      = 1;
my $qstart        = 120635196;
my $qend          = 120635214;
my $qstrand       = 1;
my $tstart        = 3;
my $tend          = 19;
my $tstrand       = 1;
my $ditag_side    = 'L';
my $ditag_pair_id = 1;
my $cigar_line    = '17M',
my $type          = "ZZ13";

my $dbID          = 4828567;
my $other_ditag   = 3278337;
#469273

my $slice         = $db->get_SliceAdaptor->fetch_by_region('chromosome', $region);
my $analysis      = $db->get_AnalysisAdaptor->fetch_by_logic_name('DitagAlign' );

######
# 1  #
######

#test constructor

my $dfa = $db->get_DitagFeatureAdaptor;
ok($dfa && ref $dfa);

######
# 2  #
######

#test construction

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


#######
# 3-4 #
#######

#test store

#hide the contents of ditag_feature table
$multi->hide('core', 'ditag_feature');

$dfa->store($feature);
ok($feature->dbID && $feature->adaptor == $dfa);

my $testfeature = $dfa->fetch_by_dbID($feature->dbID);
ok($testfeature && $testfeature->isa('Bio::EnsEMBL::Map::DitagFeature'));

#unhide table
$multi->restore('core', 'ditag_feature');

########
# 5-11 #
########

#test fetch methods

#test fetch all
my $dfs = $dfa->fetch_all();
ok(scalar @$dfs && $dfs->[0]->isa('Bio::EnsEMBL::Map::DitagFeature'));

#test fetch by dbID
my $df = $dfa->fetch_by_dbID($dbID);
ok($df && $df->isa('Bio::EnsEMBL::Map::DitagFeature') && $df->dbID == $dbID);

#test fetch by ditagID
$dfs = $dfa->fetch_by_ditagID($other_ditag);
ok((scalar @$dfs == 2) && $dfs->[0]->isa('Bio::EnsEMBL::Map::DitagFeature')
	&& $dfs->[0]->ditag_id == $other_ditag);

#test fetch by type
$dfs = $dfa->fetch_all_by_type($type);
ok((scalar @$dfs) && $dfs->[0]->isa('Bio::EnsEMBL::Map::DitagFeature')
	&& $dfs->[0]->fetch_ditag->type eq $type);

# test fetch all by slice
$slice = $db->get_SliceAdaptor->fetch_by_region('chromosome', $region,
                                                $qstart, $qend);
$dfs = $dfa->fetch_all_by_Slice($slice);
ok(scalar(@$dfs) && $dfs->[0]->isa('Bio::EnsEMBL::Map::DitagFeature'));

#test fetch_grouped
$dfs = $dfa->fetch_grouped('', $type);
ok(scalar @$dfs);
ok($dfs->[0]->{'ditag_id'} && $dfs->[0]->{'start'} && $dfs->[0]->{'end'});

######
# 12 #
######

#test list_dbIDs

my $dbIDs = $dfa->list_dbIDs();
ok(scalar @$dbIDs);

1;
