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

our $verbose = 0; #set to 1 to turn on debug printouts

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $multi->get_DBAdaptor( 'core' );


######
# 1  #
######

#test get_MarkerAdaptor
my $marker_adaptor = $db->get_MarkerAdaptor;
ok($marker_adaptor && ref $marker_adaptor);

#####
# 2 #
#####

#test fetch_by_dbID
my $marker = $marker_adaptor->fetch_by_dbID(1);

ok($marker->dbID == 1 && 
   $marker->left_primer && 
   $marker->right_primer &&
   $marker->min_primer_dist &&
   $marker->max_primer_dist);


#######
# 3-5 #
#######

#test fetch_all_by_synonym

my $markers = $marker_adaptor->fetch_all_by_synonym('D1S243');
ok(scalar(@$markers) == 1 && 
   $markers->[0]->dbID == 1 &&
   $markers->[0]->left_primer eq $marker->left_primer);

$markers = $marker_adaptor->fetch_all_by_synonym('Z16979', 'genbank');
ok(scalar(@$markers == 1));

$markers = $marker_adaptor->fetch_all_by_synonym('Z16979', 'other');
ok(scalar(@$markers == 0));

#####
# 6 #
#####

$marker_adaptor->fetch_attributes($marker);
ok(scalar(@{$marker->get_all_MarkerSynonyms}) == 9 &&
   scalar(@{$marker->get_all_MapLocations}) == 2);

#####
# 7 #
#####

#test fetch_all

$markers = $marker_adaptor->fetch_all;

debug('expecting 100 markers, got ' . scalar(@$markers));
ok(scalar(@$markers) == 100);

done_testing();
