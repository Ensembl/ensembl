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
use Bio::EnsEMBL::Map::MarkerFeature;
use Bio::EnsEMBL::Map::MarkerSynonym;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0; #set to 1 to turn on debug printouts

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $multi->get_DBAdaptor( 'core' );

######
# 1  #
######

#test constructor

my $adaptor = $db->get_MarkerFeatureAdaptor;
my $dbID = 111;
my $start = 100;
my $end   = 10;
my $slice = 
  $db->get_SliceAdaptor->fetch_by_region('contig','AL359765.6.1.13780');

my $analysis = Bio::EnsEMBL::Analysis->new;
my $marker_id = 1;
my $mapweight = 1;


#####
# 1 #
#####

#test construction
my $mf = Bio::EnsEMBL::Map::MarkerFeature->new
  ($dbID, $adaptor, $start, $end, $slice, $analysis, $marker_id, 
   $mapweight);

ok($mf && ref $mf && $mf->isa('Bio::EnsEMBL::Map::MarkerFeature'));


#######
# 2-3 #
#######

#test dbID
ok($dbID == $mf->dbID);
ok(&test_getter_setter($mf, 'dbID', undef));

#######
# 4-5 #
#######

#test adaptor
ok($adaptor == $mf->adaptor);
ok(&test_getter_setter($mf, 'adaptor', undef));

#######
# 6   #
#######

#test marker lazy-loading

my $marker = $mf->marker;
ok($marker->dbID == $marker_id);

#######
# 7-10#
#######

#test contig, start, end, strand (inherited)

ok($slice == $mf->slice);
ok($start == $mf->start);
ok($end == $mf->end);
ok($mf->strand == 0);

#######
#  11 #
#######

ok($mapweight == $mf->map_weight);



my $ms = Bio::EnsEMBL::Map::MarkerSynonym->new(1234, 'unists', 'a marker');

$mf->marker()->display_MarkerSynonym($ms);
ok($mf->display_id() eq $mf->marker()->display_MarkerSynonym()->name());

done_testing();
