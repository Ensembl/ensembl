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
use warnings;

use Test::More;
use Test::Warnings;
use Test::Exception;

use Bio::EnsEMBL::Test::MultiTestDB;

# Test to ensure that fetch_all_nearest_by_Feature calls return the correct choice all the time
# Test data is based off empty test DB and SimpleFeatures to simplify the problem of edge cases

my $db  = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $db->get_DBAdaptor('empty');
my $dbc = $dba->dbc(); 

my $analysis  = $dba->get_AnalysisAdaptor()->fetch_by_dbID(4);
my $sfa       = $dba->get_SimpleFeatureAdaptor();
my $sa        = $dba->get_SliceAdaptor();
my $ref_slice = $sa->fetch_by_seq_region_id(469294);#chromosome:NCBI33:Y:1:58368225:1
my $par_slice = $sa->fetch_by_seq_region_id(469283);#chromosome:NCBI33:20:1:62842997:1

note("Main test slice(par):\t".$par_slice->name);
note("Alt  test slice(ref):\t".$ref_slice->name);

$db->hide('empty', 'simple_feature'); #Create an empty back up.

my $a = Bio::EnsEMBL::SimpleFeature->new(
  -start => 30300100,
  -end => 30300110,
  -strand => 1,
  -slice => $par_slice,
  -analysis => $analysis,
  -display_label => 'a',
  #midpoint = 30300105
  );

# Downstream of a, overlapping
my $b = Bio::EnsEMBL::SimpleFeature->new(
  -start => 30300107,
  -end => 30300117,
  -strand => 1,
  -slice => $par_slice,
  -analysis => $analysis,
  -display_label => 'b',
  #midpoint = 30300112
  );

# Downstream of a, by a long way
my $c = Bio::EnsEMBL::SimpleFeature->new(
  -start => 30300260,
  -end => 30300270,
  -strand => 1,
  -slice => $par_slice,
  -analysis => $analysis,
  -display_label => 'c',
  #midpoint = 30300265
  );

# Feature on an assembly exception
my $d = Bio::EnsEMBL::SimpleFeature->new(
  -start => 9999950,
  -end => 9999951,
  -strand => 1,
  -slice => $ref_slice,
  -analysis => $analysis,
  -display_label => 'd',
  #midpoint = 9999950
  );

# Feature on an assembly exception, reverse stranded
my $e = Bio::EnsEMBL::SimpleFeature->new(
  -start => 9999850,
  -end => 9999851,
  -strand => -1,
  -slice => $ref_slice,
  -analysis => $analysis,
  -display_label => 'e',
  #midpoint = ??
  );

my $f = Bio::EnsEMBL::SimpleFeature->new(
  -start => 30300050,
  -end => 30300160,
  -strand => 1,
  -slice => $par_slice,
  -analysis => $analysis,
  -display_label => 'f',
  #midpoint = 30300105
  );



# Test feature store
my @simple_features = ($a,$b,$c,$d,$e,$f);
$sfa->store(@simple_features);
cmp_ok(scalar(@{$sfa->fetch_all}),'==',6,'verify successful storage of test features');

# Test distance calculations
cmp_ok($sfa->_compute_midpoint($a), '==', 30300105, 'Midpoints calculated correctly');

my ($distance,$effective) = $sfa->_compute_nearest_end($a->start, 
                                          $sfa->_compute_midpoint($a), 
                                          $a->end, 
                                          $c->start,
                                          $sfa->_compute_midpoint($c),
                                          $c->end                     );

cmp_ok($distance, '==', 150, 'A->C distance correct, i.e. nearest edge to reference middle');
my @results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $a, -UPSTREAM => 1, -LIMIT => 2) };

print_what_you_got(\@results);
cmp_ok($results[0]->[1], '==', 0,'Nearest feature is very very close indeed'); # Found $f
# cmp_ok($results[1]->[2], '==', 0,'Next nearest feature is also exceedingly close'); # Found $f

@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $a,-DOWNSTREAM => 1, -LIMIT => 5) };
cmp_ok(scalar(@results), '==', 3,'Found all features downstream');
is($results[0]->[0]->display_id, 'f','Found enveloping feature first');

#arguable b is closer here on midpoint.

# print_what_you_got(\@results);

@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $a, -UPSTREAM => 1, -LIMIT => 5) };
is($results[0]->[0]->display_id, 'f','Found enveloping feature first upstream');
# print_what_you_got(\@results);
cmp_ok(scalar(@results), '==', 1,'Found only one overlapping feature upstream');

@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $a, -UPSTREAM => 1, -LIMIT => 5, -NOT_OVERLAPPING => 1) };
cmp_ok(scalar(@results), '==', 0, 'Found nothing upstream after excluding overlaps');

@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $a, -DOWNSTREAM => -1, -LIMIT => 5, -NOT_OVERLAPPING => 1) };
# print_what_you_got(\@results);
is($results[0]->[0]->display_id, 'c','Found downstream feature first downstream without overlaps');
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

@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $d, -LIMIT => 6, -SAME_STRAND => 1) };
cmp_ok(scalar(@results), '==', 5, 'Strand restriction removes reverse stranded feature');

@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $e, -LIMIT => 6, -OPPOSITE_STRAND => 1) };
# print_what_you_got(\@results);
cmp_ok(scalar(@results), '==', 1, 'Reverse strand restriction removes all candidates');

@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $e, -LIMIT => 6, -UPSTREAM => 1, -THREE_PRIME => 1) };
print_what_you_got(\@results);
cmp_ok(scalar(@results), '==', 5, 'Upstream of reverse stranded feature gives all five features');

# Test iterative search

@results = @{ $sfa->fetch_all_by_outward_search(-FEATURE => $a, -LIMIT => 5, -RANGE => 10, -MAX_RANGE => 160)};
cmp_ok($results[3]->[0]->display_id, 'eq', 'c', 'Check for a feature found outside its initial search range');


################################ NJs Tests ################################

# DO NOT DELETE THE FOLLOWING NOTES !!!!

# Distance metric requirements
#
# A comparable metric is used for all distances reported. i.e. ordering of features sharing 
# the same distance (0bp or otherwise) may use a different metric, so long as it is internally 
# consistent and logical.
#

# Midpoint issues
#
# 2 cannot distinguish between features sharing the same midpoint but with different 
# start/end loci. b is clearly closer to q than a and c, but midpoints are the same:
#      >a>
#    >>>b>>>
#            >q>
#     <<c<<    
# KT - sorted with additional criteria.

# 3 can report erroenous results. it's easy to spot that the distances of a, b and 
#   c should be 5, 2 and 3 repectively. midpoint will give same distance to each.
#   and in the following example will even select b which is obviously futher away:
#   >>>>>>a>>>>>>
#           >b>
#                  <<<q<<<

# 4 Distance can currently seem to be 1 out due to rounding of the midpoint for features with 
#   an even length.
# KT - true, but not fatal for ordering purposes. It is consistent.

# Primes in parenthesis would be specified by -PRIME, hence over-riding the default:
#   Upstream:    Query 5' to Target 3'(5')
#   Downstream:  Query 3' to Target 5'(3')

# Features overlapping query
# 
# Order overlapping features sensibly is non-trivial, even if we measure from midpoint. 
# Consider closest of a..d upstream to Q:
# >>>a>>>
#   >>>>>>b>>>>>>
#  >>>>>>>>>>>>c>>>>>>>>>>>>
#      >>dd>>
#          >>>>>>>>>>f>>>>>>>>>>>  
#                >>>>g>>>>>
#                >>>>>>>>>>>>>h>>>>>>>>>>  
#      >>>>>>>>>>Q>>>>>>>>
#Bearing in mind the final distance reported should be 0 for overlaps, or d 5' to target 3' where target 3' 
#does not overlap.
#suggested target 5' order, so: b c a.

#Now consider streamless/primeless query Q:
#      >>a>>
#       >>>>>>>>>b>>
#       >c>
#     >>>>>>>Q>>>>>>>>
#  >>>>>>d>
#    >>>e>>
#   >>>>>>>f>>>>>>>>>>>>>>
# >>>>>>>>>g>>>>>>>>>>>

### Issues/Todos ###
#
#6 For streamed queries we probably do not want to return features where the target measurement point is enclosed
#  by the source feature. This suggests we should change the functionality of -NOT_OVERLAPPING to -INCLUDE_ENCLOSED
#  Including overlaps should be the default for non-streamed queries

# 10 STRAND could be renamed, as this is a little ambiguous to. Could be -TARGET_STRAND ?
#    Also take +/- aswell as 1/-1?


note("Nathan's tests");
my $g = Bio::EnsEMBL::SimpleFeature->new(
  -start => 30300121,
  -end => 30300128,
  -strand => 1,
  -slice => $par_slice,
  -analysis => $analysis,
  -display_label => 'g',
  );
  

my $h= Bio::EnsEMBL::SimpleFeature->new(
  -start => 30300116,
  -end => 30300141,
  -strand => 1,
  -slice => $par_slice,
  -analysis => $analysis,
  -display_label => 'h',
  );

@simple_features = ($g, $h);
$sfa->store(@simple_features);
cmp_ok(scalar(@{$sfa->fetch_all}),'==',8,'verify successful storage of more test features');

### Approximate feature landscape
# Nbp represents omitted bps, so >>>>10bp>>>, represents 17bps.
# >n> represents midpoint, or >nn> if it stradles two bases
#  chromosome:NCBI33:20:1:62842997:1
#                
#a               >>>>aa>>>>        
#b                     >>>>bb>>>> 
#c                                                                 95bp       >>>>cc>>>>
#f >40bp>>>>>>>>>>>>>>f>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#g                                  >>>g>>> 
#h                             >>>>>>>>>>>hh>>>>>>>>>>>>   


@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $h,
                                                 -UPSTREAM  => 1, 
                                                 -LIMIT   => 5  ) };
print_what_you_got(\@results);

# Results:
# Feature: (b)Downstream overlap at -11. Overlap? 1
# Feature: (f)Enveloping at -23. Overlap? 1
# Feature: Test me! at -18. Overlap? 0

#g is due to it being right next to the midpoint. But it is not upstream
# KT agreed.

cmp_ok(scalar(@results), '==', 3, 'Found 3 features which have some region upstream');
is($results[0]->[0]->display_id, 'b', 'Closest feature(b) is overlapping source 5\'');
#b is only closest if we implement 5' ordering of overlapping upstream features
is($results[1]->[0]->display_id, 'f', 'Next closest feature(f) is overlapping source 5\'');


ok((($results[0]->[1] == 0) &&
    ($results[1]->[1] == 0)   ), 'Overlapping upstream features have distance of 0');

#Do we need to add it another overlapping feature with identical start to f , but a closer 3', to test 2ndary sort

is($results[2]->[0]->display_id, 'a', 'Next closest upstream feature(a) does not overlap source 5\'');

cmp_ok(abs($results[2]->[1]), '==', 6, 'Closest non-overlapping upstream feature(a) is 6bp? away (h 5\' to a 3\')'); # 5bp is the separation, between boundaries


@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE     => $h,
                                                 -UPSTREAM    => 1, 
                                                 -THREE_PRIME => 1, 
                                                 -LIMIT       => 5  ) };
print_what_you_got(\@results);

#f should not be reported at all
#a should be first with 5 i.e. measuring from h 5' to a 3'? 

cmp_ok(scalar(@results), '==', 1, 'Found 1 upstream feature(a) 3\'');

#Refresh the table
$db->restore('empty', 'simple_feature');
$db->hide('empty', 'simple_feature');


my $i = Bio::EnsEMBL::SimpleFeature->new(
  -start => 1,
  -end => 13,
  -strand => 1,
  -slice => $par_slice,
  -analysis => $analysis,
  -display_label => 'i',
  );

my $j = Bio::EnsEMBL::SimpleFeature->new(
  -start => 1,
  -end => 13,
  -strand => 1,
  -slice => $par_slice,
  -analysis => $analysis,
  -display_label => 'j',
  );

my $query = Bio::EnsEMBL::SimpleFeature->new(
  -start => 1,
  -end => 13,
  -strand => 1,
  -slice => $par_slice,
  -analysis => $analysis,
  -display_label => 'Q',
  );

$sfa->store($i, $j, $query);


#These are now not valid as we have changed the way overlapping features are measured
#see below
#@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE     => $query ) };
#print_what_you_got(\@results);
#is($results[0]->[0]->display_id, 'j', 'Closest feature(j) with shared mid-point, based on closer 5\'');
#cmp_ok($results[0]->[2], '==', 0, 'Distance is based on k 5\' to  j 5\'');

#@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE     => $query,
#                                                 -THREE_PRIME => 1 ) };
#print_what_you_got(\@results);
#is($results[0]->[0]->display_id, 'j', '3\' query gives closest feature(j) with overlapping 3\'');
#cmp_ok($results[0]->[2], '==', 0, '3\' query with overlapping 3\' distance is 0');
#cmp_ok($results[1]->[2], '==', 1, '3\' query with non-overlapping 3\' distance is 1');
#should be j=0 then i=1



my $k = Bio::EnsEMBL::SimpleFeature->new(
  -start => 15,
  -end => 17,
  -strand => -1,
  -slice => $par_slice,
  -analysis => $analysis,
  -display_label => 'k',
  );

my $l = Bio::EnsEMBL::SimpleFeature->new(
  -start => 15,
  -end => 17,
  -strand => 1,
  -slice => $par_slice,
  -analysis => $analysis,
  -display_label => 'l',
  );

my $m = Bio::EnsEMBL::SimpleFeature->new(
  -start => 15,
  -end => 19,
  -strand => 1,
  -slice => $par_slice,
  -analysis => $analysis,
  -display_label => 'm',
  );

$sfa->store($k, $l ,$m);

# >>>>>>i>>>>>>
#     >>j>>
#               <k<
#               >l>
#               >>m>>
#     >>Q>>


@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $query,
                                                 -DOWNSTREAM  => 1, 
                                                 -LIMIT => 5) };
print_what_you_got(\@results);

cmp_ok(scalar(@results), '==', 3, 'Found 3 downstream features');


@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $query,
                                                 -DOWNSTREAM  => 1,
                                                 -SAME_STRAND  => 1 , 
                                                 -LIMIT => 5) };

cmp_ok(scalar(@results), '==', 2, 'Found 2 downstream +ve strand features');

@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $query,
                                                 -DOWNSTREAM  => 1,
                                                 -OPPOSITE_STRAND  => 1, 
                                                 -LIMIT => 5) };
print_what_you_got(\@results);
cmp_ok(scalar(@results), '==', 1, 'Found 1 downstream +ve strand features');


is($results[0]->[0]->display_id, 'k', 'Closest downstream -ve strand feature(m),');
cmp_ok($results[0]->[1], '==', 2, 'Closest downstream -ve strand distance is based on 3\' to 3\'');


sub print_what_you_got {
  my $results = shift;
  note ("Results: ".scalar @$results." features");
  if (scalar(@$results) == 0) { note("No hits"); return;}
  for (my $i =0; $i<scalar(@$results);$i++) {
    my ($feature,$distance) = @{$results->[$i]};
    no warnings 'uninitialized';
    note("Feature: ".$feature->display_id." at ".$distance);
  }
}

# test utility functions get_nearest_Gene and 
$dba = $db->get_DBAdaptor('core');
my $tra = $dba->get_TranscriptAdaptor();
my $ga = $dba->get_GeneAdaptor();
my $enst = $tra->fetch_by_stable_id('ENST00000300415');
my $gene;
($gene,$distance) = @{ $ga->fetch_nearest_by_Feature($enst) };
is($gene->stable_id, 'ENSG00000101331','Simple usecase for fetch_nearest_by_Feature');
my $ensg = $enst->get_nearest_Gene;
# diag(Dumper $ensg);
is($ensg->stable_id,'ENSG00000101331','Simple usecase for get_nearest_Gene');

done_testing;

