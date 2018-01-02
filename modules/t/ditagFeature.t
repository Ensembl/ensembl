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
# 2-15 #
########

#test ditag_id, ditag_side, hit_start, hit_end,
# hit_strand, cigar_line, start, end, strand,
# dbID, sequence, slice, ditag_pair_id

my $ditagFeatures = $dfa->fetch_all_by_ditagID($ditag_id);
my ($ditagFeature)  = grep { $_->cigar_line() eq $cigar_line } @{$ditagFeatures};

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

#########
# 16-17 #
#########

#test fetch_ditag

my $ditag = $ditagFeature->ditag();

ok(defined $ditag && $ditag->isa('Bio::EnsEMBL::Map::Ditag'));
ok($ditag->dbID == $ditag_id);

######
# 18 #
######

#test get_ditag_location

my ($rstart, $rend, $rstrand) = $ditagFeature->get_ditag_location();
ok(defined($rstart) && defined($rend) && defined($rstrand));

done_testing();
