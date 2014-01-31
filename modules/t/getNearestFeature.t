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
);

my $b = Bio::EnsEMBL::SimpleFeature->new(
    -start => 30300107,
    -end => 30300117,
    -strand => 1,
    -slice => $par_slice,
    -analysis => $analysis,
    -display_label => 'Downstream overlap',
);

my $c = Bio::EnsEMBL::SimpleFeature->new(
    -start => 30300260,
    -end => 30300270,
    -strand => 1,
    -slice => $par_slice,
    -analysis => $analysis,
    -display_label => 'Downstream far',
);

my $d = Bio::EnsEMBL::SimpleFeature->new(
    -start => 9999950,
    -end => 9999951,
    -strand => 1,
    -slice => $ref_slice,
    -analysis => $analysis,
    -display_label => 'Upstream far, assembly exception',
);

my $e = Bio::EnsEMBL::SimpleFeature->new(
    -start => 9999850,
    -end => 9999851,
    -strand => -1,
    -slice => $ref_slice,
    -analysis => $analysis,
    -display_label => 'Reverse strand',
);

my $f = Bio::EnsEMBL::SimpleFeature->new(
    -start => 30300050,
    -end => 30300160,
    -strand => 1,
    -slice => $par_slice,
    -analysis => $analysis,
    -display_label => 'Enveloping',
);

my @simple_features = ($a,$b,$c,$d,$e,$f);
$sfa->store(@simple_features);

cmp_ok(scalar(@{$sfa->fetch_all}),'==',6,'verify successful storage of test features');
#($self, $feat, $prime, $stranded, $stream, $num_feats, $max_dist, $stream_from_midpoint)
my ($results,$distances) = @{ $sfa->fetch_nearest_by_Feature($a,undef,undef,undef,undef,1000,undef) };

#note(dump($results) );
note(dump($distances));
done_testing;
