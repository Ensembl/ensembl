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
use Bio::EnsEMBL::Test::MultiTestDB;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $multi->get_DBAdaptor( "core" );

my $sa = $db->get_SliceAdaptor();
my $ga = $db->get_GeneAdaptor();

ok(!$ga->db()->no_cache(), 'Checking cache is on');

my $cache_assert = sub {
  my ($expected) = @_;
  is(scalar(keys %{$ga->{_slice_feature_cache}}), $expected, sprintf('Asserting cache has %d element(s)', $expected));
};

my $run = sub {
  my $start = 30_249_935;
  my $end = 31_254_640;
  my $offset = 0; 
  my @regions = (
    [$start, $end + $offset++],
    [$start, $end + $offset++],
    [$start, $end + $offset++],
    [$start, $end + $offset++],
    [$start, $end + $offset++]
  );
  $ga->fetch_all_by_Slice($sa->fetch_by_region( "chromosome", "20", @{$regions[0]} ));
  
  $cache_assert->(1);  
  foreach my $region (@regions) {
    my $slice = $sa->fetch_by_region( "chromosome", "20", @{$region} );
    my $features = $ga->fetch_all_by_Slice($slice);
  }
  $cache_assert->(4);
};

$run->();
$ga->clear_cache();
$run->();

done_testing();
