# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
use warnings;

use Test::More;
use Test::Warnings;
use Test::Exception;
use Data::Dump::Color qw/dump/;

use Bio::EnsEMBL::Test::MultiTestDB;

# Test to ensure that get_nearest_Feature calls return the correct choice all the time
# Test data is based off empty test DB and SimpleFeatures to simplify the problem of edge cases

my $db = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $db->get_DBAdaptor('empty');
my $dbc = $dba->dbc();

my $sfa = $dba->get_SimpleFeatureAdaptor();
my $sa = $dba->get_SliceAdaptor();

my $ref_slice = $sa->fetch_by_seq_region_id(469294);
my $par_slice = $sa->fetch_by_seq_region_id(469283);

my $analysis = $dba->get_AnalysisAdaptor()->fetch_by_dbID(4);

my $a = Bio::EnsEMBL::SimpleFeature->new(
    -start => 30300100,
    -end => 30300110,
    -strand => 1,
    -slice => $par_slice,
    -analysis => $analysis,
    -display_label => 'Test me!',
    #midpoint = 30300105
);

my $b = Bio::EnsEMBL::SimpleFeature->new(
    -start => 30300107,
    -end => 30300117,
    -strand => 1,
    -slice => $par_slice,
    -analysis => $analysis,
    -display_label => 'Downstream overlap',
    #midpoint = 30300112
);

my $c = Bio::EnsEMBL::SimpleFeature->new(
    -start => 30300260,
    -end => 30300270,
    -strand => 1,
    -slice => $par_slice,
    -analysis => $analysis,
    -display_label => 'Downstream far',
    #midpoint = 30300265
);

my $d = Bio::EnsEMBL::SimpleFeature->new(
    -start => 9999950,
    -end => 9999951,
    -strand => 1,
    -slice => $ref_slice,
    -analysis => $analysis,
    -display_label => 'Upstream far, assembly exception',
    #midpoint = 9999950
);

my $e = Bio::EnsEMBL::SimpleFeature->new(
    -start => 9999850,
    -end => 9999851,
    -strand => -1,
    -slice => $ref_slice,
    -analysis => $analysis,
    -display_label => 'Reverse strand',
    #midpoint = ??
);

my $f = Bio::EnsEMBL::SimpleFeature->new(
    -start => 30300050,
    -end => 30300160,
    -strand => 1,
    -slice => $par_slice,
    -analysis => $analysis,
    -display_label => 'Enveloping',
    #midpoint = 30300105
);

my @simple_features = ($a,$b,$c,$d,$e,$f);
$sfa->store(@simple_features);

cmp_ok(scalar(@{$sfa->fetch_all}),'==',6,'verify successful storage of test features');
#($self, $feat, $prime, $stranded, $stream, $num_feats, $max_dist, $stream_from_midpoint)
my ($results,$distances) = @{ $sfa->fetch_all_nearest_by_Feature($a,undef,undef,-1,2,1000,1) };
foreach (@$results) {
    note($_->display_label);
}
foreach (@$distances) {
    note($_);
}

# Test primeless (= test from midpoint), both strands, downstream, range 1000
($results,$distances) = @{ $sfa->fetch_all_nearest_by_Feature($a,undef,undef,-1,3,1000,undef) };
print $results->[0]->display_label . " overlaps ";
print $distances->[0] . "\n";
print $results->[1]->display_label . " overlaps ";
print $distances->[1] . "\n";
print $results->[2]->display_label . " overlaps ";
print $distances->[2] . "\n";
is($results->[0]->display_label, 'Test me!', 'First feature is itself');
is($results->[1]->display_label, 'Downstream overlap', 'Second feature is closest downstream');
is($results->[2]->display_label, 'Enveloping', 'Third feature is enveloping');
is($distances->[0], '0', 'Closest feature is on the same point');
is($distances->[1], '7', 'Downstream is 7 away');
is($distances->[2], '50', 'Enveloping end is 50 away');
($results, $distances) = @{ $sfa->fetch_all_nearest_by_Feature($a, undef, undef, 1, 2, 1000, undef) };
is($results->[0]->display_label, 'Test me!', 'First feature is itself');
is($results->[1]->display_label, 'Enveloping', 'Second feature upstream is enveloping');
is($distances->[1], '50', 'Enveloping is 50 away');

note(dump($distances));
done_testing;
