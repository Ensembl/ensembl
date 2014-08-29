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

my $db  = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $db->get_DBAdaptor('empty');
my $dbc = $dba->dbc(); 
#NJ Testing connection? If so need to get db_handle too? Should also be eval'd as a test
#Or does MutlTestDB already do this, if so drop.

my $analysis  = $dba->get_AnalysisAdaptor()->fetch_by_dbID(4);
my $sfa       = $dba->get_SimpleFeatureAdaptor();
my $sa        = $dba->get_SliceAdaptor();
my $ref_slice = $sa->fetch_by_seq_region_id(469294);#chromosome:NCBI33:Y:1:58368225:1
my $par_slice = $sa->fetch_by_seq_region_id(469283);#chromosome:NCBI33:20:1:62842997:1

#Using warn here causes not ok 19 - no (unexpected) warnings (via done_testing)
note("Main test slice(par):\t".$par_slice->name);
note("Alt  test slice(ref):\t".$ref_slice->name);


#NJ Probably need to use hive/restore to manage table contents between tests
#So let's put in a hide here for the empty table, just for convenience?
$db->hide('empty', 'simple_feature'); #Create an empty back up.

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



# Test feature store
my @simple_features = ($a,$b,$c,$d,$e,$f);
$sfa->store(@simple_features);
cmp_ok(scalar(@{$sfa->fetch_all}),'==',6,'verify successful storage of test features');

# Test distance calculations
cmp_ok($sfa->_compute_midpoint($a), '==', 30300105, 'Midpoints calculated correctly');

my $distance = $sfa->_compute_nearest_end($a->start, 
                                          $sfa->_compute_midpoint($a), 
                                          $a->end, 
                                          $c->start,
                                          $sfa->_compute_midpoint($c),
                                          $c->end                     );

cmp_ok($distance, '==', 155, 'A->C distance correct, i.e. nearest edge to reference middle');
my @results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $a, -STREAM => -1, -LIMIT => 2) };

cmp_ok($results[0]->[2], '==', 0,'Nearest feature is very very close indeed');
cmp_ok($results[1]->[2], '==', 0,'Next nearest feature is also exceedingly close');

# print_what_you_got(\@results);


@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $a,-STREAM => -1, -LIMIT => 5) };
cmp_ok(scalar(@results), '==', 4,'Found all features downstream');
is($results[0]->[0]->display_id, 'Enveloping','Found enveloping feature first');
#NJ what is distance here? If we are measuring down stream, what are we measuring to and from? Should be midpoint?
#arguable b is closer here on midpoint.

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

#NJ What are the distances here?

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






#NJs Tests

# DO NOT DELETE THE FOLLOWING NOTES !!!!

# Distance metric requirements
#
# A comparable metric is used for all distances reported. i.e. ordering of features sharing 
# the same distance (0bp or otherwise) may use a different metric, so long as it is internally 
# consistent and logical.
#

# Midpoint issues
#
# 1 Is not the expected distance for the user, which would be to the nearest end, unless
#   a prime has been specified
#
# 2 Cannot distinguish between features sharing the same midpoint but with different 
# start/end loci. b is clearly closer to Q than a and c, but midpoints are the same:
#      >a>
#    >>>b>>>
#            >Q>
#     <<c<<
#
# 3 Can report erroenous results. It's easy to spot that the distances of a, b and 
#   c should be 5, 2 and 3 repectively. Midpoint will give same distance to each.
#   And in the following example will even select b which is obviously futher away:
#   >>>>>>a>>>>>>
#           >b>
#                  <<<Q<<<
#
# 4 Distance can currently seem to be 1 out due to rounding of the midpoint for features with 
#   an even length.
#
# Given these issues, it is best to measure to the specified /nearest end or 'prime2prime'.
# Primes in parenthesis would be specified by -PRIME, hence over-riding the default:
#   Upstream:    Query 5' to Target 3'(5')
#   Downstream:  Query 3' to Target 5'(3')
# 

# Features overlapping query
# 
# Given the obvious flaws of the mid-point metric and the suggested new prime2prime 
# measurement, features overlapping the query pose somewhat of a problem. It's extremely
# hard to define a point to measurement to for PRIMEless queries and a point to measure
# from for STREAMless queries. The problem gets even more gnarly when considerng query 
# enclosed and query enclosing features. I would not suggest working through the case below,
# rather just accept, that all overlap distance should be set to 0, and it's just not worth
# working out any scheme to order in, as there is no sensible way to do it which will be 
# generic enough to apply in all contexts and be internally consistent. Which ever scheme you 
# pick at least some of the ordering comes out whacky.
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
#What measure should be used?
#Bearing in mind the final distance reported should be 0 for overlaps, or d 5' to target 3' where target 3' 
#does not overlap.
#The only possible measure seems to be target 5' order, so: b c a.
#This is at least standardised/preditable/testable, even if the logic is questionable.
#What if feature share a start loci? FFS! Then 2ndary order key of 3'?

#Now consider streamless/primeless query Q:
#      >>a>>
#       >>>>>>>>>b>>
#       >c>
#     >>>>>>>Q>>>>>>>>
#  >>>>>>d>
#    >>>e>>
#   >>>>>>>f>>>>>>>>>>>>>>
# >>>>>>>>>g>>>>>>>>>>>


#Measure closest primes between query and target? 
#Then a would be closer, but b is arguably closer here due to coverage.
#mid point would resolve this but cannot resolve a vs c
#what about summing the differences and sorting based on that?
# in that case the order would be b a c. 
#
#but how are these resolved when we have partial overlaps, such as d
#currently partial overlaps and enclosures are not considered together
#is c closer than d or e? which has a 5'-5' distance of just 1?
# is f closer than d? is g closer than d?
#so there is definitely an intersect between apparent closeness between enclosing features and partial overlaps
#
#This is way too gnarly, just set to 0 and return in arbitrary order for all overlaps?
#or just use mid-point of ordering? despite it's inability to resolve features whith shared midpoint.

### Issues/Todos ###
#
#1 It seems we're getting the source feature back. This needs filtering out.
#2 Ordering is wrong and seemingly arbitrary when dealing with negative numbers. Drop
#  negative numbers, such that the distances are standardised.
#  This will be useful for filtering, which will be a common task, and leave up/downsteam
#  detection to users, as this will be much less common?
#4 Distances always seem to be from midpoint. Convert these to prime2prime.
#5 Distances should be set to 0 for overlapping features, as it is non-trivial to define what is nearest. 
#  I suspect this is already done, but I haven't dug into the code to see whether this is always the case.
#  (Could optionally order features sharing same distance (most likely 0bp) by midpoint, although gives ambiguous results)
#6 For streamed queries we probably do not want to return features where the target measurement point is enclosed
#  by the source feature. This suggests we should change the functionality of -NOT_OVERLAPPING to -INCLUDE_ENCLOSED
#  Including overlaps should be the default for non-streamed queries
# 8 Strand tests seem lacking, so we need to add a target feature on the negative strand for each context?
# 9 Change -STREAM values to up|down
# 10 STRAND could be renamed, as this is a little ambiguous to. Could be -TARGET_STRAND ?
#    Also take +/- aswell as 1/-1?
# 11 Change -THREE/FIVEPRIME to -TARGET_PRIME and take values of 5|3. This should simplify the interface
#    and make it more obviosu which prime we're dealing with here.
# 11 Simplify return type drop overlap [ [$feature, $overlap, $distance] ... ]
#    Distance of 0 indicates overlap
#    If we features count exceeds limit due to some sharing the same distance, still return all, an dlet the user handle it.
# 12 Test search range
# 13 Add iterator
# 13 Do we need to be able pass a list of features, rather than querying for the, such that we can pass a heterogenous
#    set of features?
#    This will require moving part of the code to a calculate_distances(probably already exists), which could be in a separate
#    utils module?



#Adding features here so we don't change the results of the previous tests
#Using the feature var names as the display_ids the previous display_ids only make sense 
#in the context of a given test and it also makes it easier for reference to the feature image mockup

my $g = Bio::EnsEMBL::SimpleFeature->new(
  -start => 30300121,
  -end => 30300128,
  -strand => 1,
  -slice => $par_slice,
  -analysis => $analysis,
  -display_label => 'g',
  );
  

my $h= Bio::EnsEMBL::SimpleFeature->new(-start => 30300116,
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
#NJs new features, should we do this on -1?
#g                                  >>>g>>> 
#h                             >>>>>>>>>>>hh>>>>>>>>>>>>   
#It would probably be useful to copy this before each test and delete the unused features

@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $h,
                                                 -STREAM  => 1, 
                                                 -LIMIT   => 5  ) };
print_what_you_got(\@results);

# Results:
# Feature: h at 0. Overlap? 1
# Feature: g at 0. Overlap? 1 
# Feature: (b)Downstream overlap at -11. Overlap? 1
# Feature: (f)Enveloping at -23. Overlap? 1
# Feature: Test me! at -18. Overlap? 0

#g is due to it being right next to the midpoint. But it is not upstream
#so the code isn't doign the right filtering when a stream is defined.
#The distances for b and f make no sense to me
#also should we even have a - value prefix? This might be useful for differentiatiing up/down stream but makes it harder to filter.
#probably better to let the use figure that out?


cmp_ok(scalar(@results), '==', 3, 'Found 3 features which have some region upstream');
is($results[0]->[0]->display_id, 'b', 'Closest feature(b) is overlapping source 5\'');
#b is only closest if we implement 5' ordering of overlapping upstream features
is($results[1]->[0]->display_id, 'f', 'Next closest feature(f) is overlapping source 5\'');


ok((($results[0]->[2] == 0) &&
    ($results[1]->[2] == 0)   ), 'Overlapping upstream features have distance of 0');

#Do we need to add it another overlapping feature with identical start to f , but a closer 3', to test 2ndary sort

is($results[2]->[0]->display_id, 'a', 'Next closest upstream feature(a) does not overlap source 5\'');

#probably need skip statement here based on success of previous test?
cmp_ok($results[2]->[2], '==', 5, 'Closest non-overlapping upstream feature(a) is 5bp away (h 5\' to a 3\')');


@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE     => $h,
                                                 -STREAM      => 1, 
                                                 -THREE_PRIME => 1, 
                                                 -LIMIT       => 5  ) };
print_what_you_got(\@results);

# Results:
# Feature: (f)Enveloping at 32. Overlap? 1
# Feature: (a)Test me! at -18. Overlap? 0

#f should not be reported at all
#a should be first with 5 i.e. measuring from h 5' to a 3'? 
#a -18 is correct when measuring from midpoint.

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
#         >>Q>>


@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $query,
                                                 -STREAM  => -1, ) };
print_what_you_got(\@results);

cmp_ok(scalar(@results), '==', 3, 'Found 3 downstream features');

#test distances here are all 1?

@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $query,
                                                 -STREAM  => -1,
                                                 -STRAND  => 1 ) };

cmp_ok(scalar(@results), '==', 2, 'Found 2 downstream +ve strand features');

#test display labels are correct and distances here are all 1?

@results = @{ $sfa->fetch_all_nearest_by_Feature(-FEATURE => $query,
                                                 -STREAM  => -1,
                                                 -STRAND  => -1 ) };

cmp_ok(scalar(@results), '==', 1, 'Found 1 downstream +ve strand features');


is($results[0]->[0]->display_id, 'm', 'Closest downstream -ve strand feature(m),');
cmp_ok($results[0]->[2], '==', 1, 'Closest downstream -ve strand distance is based on 3\' to 3\'');


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

