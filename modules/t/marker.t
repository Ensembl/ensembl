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
use Bio::EnsEMBL::Map::Marker;
use Bio::EnsEMBL::Map::MarkerSynonym;
use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0; #set to 1 to turn on debug printouts

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $multi->get_DBAdaptor( 'core' );

my $dbID = 1;
my $marker_adaptor = $db->get_MarkerAdaptor;

my $left_primer = 'ACCTGTTG';
my $right_primer = 'CCCCTTTGGG';
my $min_primer_dist = 75;
my $max_primer_dist = 100;
my $priority = 80;



######
# 1  #
######

#test the constructor
my $marker = Bio::EnsEMBL::Map::Marker->new
  ($dbID, $marker_adaptor, $left_primer, $right_primer, $min_primer_dist,
   $max_primer_dist, $priority);

ok($marker && ref $marker && $marker->isa('Bio::EnsEMBL::Map::Marker'));


#######
# 2-3 #
#######

#test dbID
ok($marker->dbID == $dbID);
ok(&test_getter_setter($marker, 'dbID', 123));

#######
# 4-5 #
#######

#test adaptor
ok($marker->adaptor == $marker_adaptor);
ok(&test_getter_setter($marker, 'adaptor', undef));

#######
# 6-7 #
#######

#test left primer
ok($marker->left_primer eq $left_primer);
ok(&test_getter_setter($marker, 'left_primer', 'ACCCTTTGGG'));

#######
# 8-9 #
#######
 
#test right primer
ok($marker->right_primer eq $right_primer);
ok(&test_getter_setter($marker, 'right_primer', 'CCTTGGGTTTT'));

#########
# 10-11 #
#########

#test min primer dist
ok($marker->min_primer_dist == $min_primer_dist);
ok(&test_getter_setter($marker, 'min_primer_dist', 130));

##########
# 12-13  #
##########

#test max primer dist
ok($marker->max_primer_dist == $max_primer_dist);
ok(&test_getter_setter($marker, 'max_primer_dist', 180));

#########
# 14-15 #
#########

#test priority
ok($marker->priority == $priority);
ok(&test_getter_setter($marker, 'priority', 50));


#######
# 16 #
######

#test add_MarkerSynonyms & get_all_MarkerSynonyms
my $ms = Bio::EnsEMBL::Map::MarkerSynonym->new(123, 'genbank', 'Z123456');

$marker->add_MarkerSynonyms($ms);
ok($marker->get_all_MarkerSynonyms->[0] == $ms);

######
# 17 #
######

# test flush markerSynonyms & get_all_MarkerSynonyms
$marker->flush_MarkerSynonyms;

#expect lazy loaded synonyms
debug("got " . scalar(@{$marker->get_all_MarkerSynonyms}) . 
      "synonyms, expected 9");
ok(scalar(@{$marker->get_all_MarkerSynonyms}) == 9);

######
# 18 #
######

#test display marker synonym
ok(&test_getter_setter($marker, 'display_MarkerSynonym', $ms));

######
# 19 #
######

#test add_MapLocations
my $ml = Bio::EnsEMBL::Map::MapLocation->new('genethon', 1, 12.3, 0.05);

$marker->flush_MapLocations;
$marker->add_MapLocations($ml);

ok($marker->get_all_MapLocations->[0] == $ml);

######
# 20 #
######

#test flush_MapLocations & get_all_MapLocations
$marker->flush_MapLocations;
debug("got " . scalar(@{$marker->get_all_MapLocations}) . 
      " locations, expected 2");
#expect lazy-loaded map locations
ok(scalar @{$marker->get_all_MapLocations} == 2 );


######
# 21 #
######

#test get_MapLocation
ok($marker->get_MapLocation('genethon') &&
   $marker->get_MapLocation('marshfield'));

#########
# 22-23 #
#########

#test get_MarkerFeatures
my @mfs = @{$marker->get_all_MarkerFeatures};

ok(scalar(@mfs) == 1);
my $mf = shift @mfs;
ok($mf->start == 5769 && $mf->end == 5959);

done_testing();
