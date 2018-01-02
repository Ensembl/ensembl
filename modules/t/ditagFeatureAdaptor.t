# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use Test::More;
use Test::Warnings;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Map::Ditag;
use Bio::EnsEMBL::Analysis;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor( 'core' );

my $region        = '11';
my $ditag_id      = 3278356;
my $qstart        = 83225874;
my $qend          = 83236347;
my $qstrand       = -1;
my $tstart        = 3;
my $tend          = 19;
my $tstrand       = 1;
my $ditag_side    = 'L';
my $ditag_pair_id = 1;
my $cigar_line    = '17M',
my $tag_library   = 'ZZ13';
my $logic_name    = 'DitagAlign';
my $dbID          = 4828567;
my $other_ditag   = 3278337;

my $slice         = $db->get_SliceAdaptor->fetch_by_region('chromosome', $region);
my $analysis      = $db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);

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
					-ditag_id      => undef,
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
$multi->hide('core', 'ditag');

#create a ditag object to attach to the ditagFeature object fist
my $name      = "101A01-2";
my $type      = "ZZ13";
my $tag_count = 2;
my $sequence  = "GAGAACTTGGACCGCAGAGAATACACACAAATCAAACC";
my $new_ditag = Bio::EnsEMBL::Map::Ditag->new (
                                               -name     => $name, 
                                               -type     => $type,
					       -count    => $tag_count,
                                               -sequence => $sequence, 
                                              );
my @ditags = ( $new_ditag );
$db->get_DitagAdaptor->store(\@ditags);

#attach to ditagFeature
$feature->ditag($new_ditag);
$feature->ditag_id($new_ditag->dbID);

$dfa->store($feature);
ok(defined $feature->dbID && $feature->adaptor == $dfa);

my $testfeature = $dfa->fetch_by_dbID($feature->dbID);
ok($testfeature && $testfeature->isa('Bio::EnsEMBL::Map::DitagFeature'));

#unhide table
$multi->restore('core', 'ditag_feature');
$multi->restore('core', 'ditag');

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
$dfs = $dfa->fetch_all_by_ditagID($other_ditag);
ok((scalar @$dfs == 2) && $dfs->[0]->isa('Bio::EnsEMBL::Map::DitagFeature')
	&& $dfs->[0]->ditag_id == $other_ditag);

#test fetch by type
$dfs = $dfa->fetch_all_by_type($tag_library);
ok((scalar @$dfs) && $dfs->[0]->isa('Bio::EnsEMBL::Map::DitagFeature')
	&& $dfs->[0]->ditag->type eq $tag_library);

# test fetch all by slice
$slice = $db->get_SliceAdaptor->fetch_by_region('chromosome', $region,
                                                $qstart, $qend);
#use slice only
$dfs = $dfa->fetch_all_by_Slice($slice);
ok(scalar(@$dfs) && $dfs->[0]->isa('Bio::EnsEMBL::Map::DitagFeature'));

#use tag-library
$dfs = $dfa->fetch_all_by_Slice($slice, $tag_library);
ok(scalar(@$dfs) && $dfs->[0]->isa('Bio::EnsEMBL::Map::DitagFeature'));

#use logic-name
$dfs = $dfa->fetch_all_by_Slice($slice, '', $logic_name);
ok(scalar(@$dfs) && $dfs->[0]->isa('Bio::EnsEMBL::Map::DitagFeature'));

######
# 12 #
######

#test list_dbIDs

my $dbIDs = $dfa->list_dbIDs();
ok(scalar @$dbIDs);

done_testing();
