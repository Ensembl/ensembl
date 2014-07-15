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

# Test distance calculations
cmp_ok($sfa->_compute_midpoint($a), '==', 30300105, 'Midpoints calculated correctly');
my $distance = $sfa->_compute_nearest_end($a->start,$sfa->_compute_midpoint($a),$a->end,$c->start,$sfa->_compute_midpoint($c),$c->end);
cmp_ok($distance, '==', 155, 'A->C distance correct, i.e. nearest edge to reference middle');
my @results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $a, -STREAM => -1, -LIMIT => 2) };

cmp_ok($results[0]->[2], '==', 0,'Nearest feature is very very close indeed');
cmp_ok($results[1]->[2], '==', 0,'Next nearest feature is also exceedingly close');

# print_what_you_got(\@results);


@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $a,-STREAM => -1, -LIMIT => 5) };
cmp_ok(scalar(@results), '==', 4,'Found all features downstream');
is($results[0]->[0]->display_id, 'Enveloping','Found enveloping feature first');
#print_what_you_got(\@results);

@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $a, -STREAM => 1, -LIMIT => 5) };
is($results[0]->[0]->display_id, 'Enveloping','Found enveloping feature first upstream');
# print_what_you_got(\@results);
cmp_ok(scalar(@results), '==', 3,'Found only overlapping features upstream');

@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $a, -STREAM => 1, -LIMIT => 5, -NOT_OVERLAPPING => 1) };
cmp_ok(scalar(@results), '==', 0, 'Found nothing upstream after excluding overlaps');

@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $a, -STREAM => -1, -LIMIT => 5, -NOT_OVERLAPPING => 1) };
is($results[0]->[0]->display_id, 'Downstream far','Found downstream feature first downstream without overlaps');
# print_what_you_got(\@results);
cmp_ok(scalar(@results),'==', 1, 'Found distant feature only downstream, after excluding overlaps');

@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $a, -LIMIT => 5, -NOT_OVERLAPPING => 1) };
# print_what_you_got(\@results);
cmp_ok(scalar(@results),'==', 1, 'Found distant feature only upstream and downstream, after excluding overlaps');

@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $a, -LIMIT => 5) };
# print_what_you_got(\@results);
cmp_ok(scalar(@results), '==', 4,'Found all features in area');

@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $a, -LIMIT => 5, -THREE_PRIME => 1) };

@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $a, -LIMIT => 5, -FIVE_PRIME => 1) };

@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $d, -LIMIT => 6) };
cmp_ok(scalar(@results), '==', 6, 'Assembly exception brings in HAP data');


@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $d, -LIMIT => 6, -STRAND => 1) };
cmp_ok(scalar(@results), '==', 5, 'Strand restriction removes reverse stranded feature');


@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $e, -LIMIT => 6, -STRAND => 1) };
# print_what_you_got(\@results);
cmp_ok(scalar(@results), '==', 1, 'Reverse strand restriction removes all candidates');

@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $e, -LIMIT => 6, -STREAM => 1, -THREE_PRIME => 1) };
# print_what_you_got(\@results);
cmp_ok(scalar(@results), '==', 5, 'Upstream of reverse stranded feature gives all five features');

sub print_what_you_got {
    my $results = shift;
    note ("Results:");
    if (scalar(@$results) == 0) { note("No hits"); return;}
    for (my $i =0 ; $i<scalar(@$results);$i++) {
        my ($feature,$overlap,$distance) = @{$results->[$i]};
        note("Feature: ".$feature->display_id." at ".$distance.". Overlap? ".$overlap);
    }
}

done_testing;
