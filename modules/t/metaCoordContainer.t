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

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose = 1;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $multi->get_DBAdaptor('core');

$multi->save('core', 'meta_coord');


#
# 1 - Can construct meta coord container
#

my $mcc = $db->get_MetaCoordContainer();
ok($mcc, 'We have a MetaCoordContainer');

{
  my @coord_systems = @{$mcc->fetch_all_CoordSystems_by_feature_type('exon')};
  is(scalar(@coord_systems), 1, 'We have only one coordinate system');
  is($coord_systems[0]->name, 'chromosome', 'Only coordinate system is chromosome');
  my $cs = $coord_systems[0];
  my $current_max = $mcc->fetch_max_length_by_CoordSystem_feature_type($cs, 'exon');
  my $less_length = $current_max - 1;
  my $greater_length = $current_max + 1;

  #Test if we try to set max to something smaller
  $mcc->add_feature_type($cs, 'exon', $less_length);
  is($mcc->fetch_max_length_by_CoordSystem_feature_type($cs, 'exon'), $current_max, 'Submission of a smaller length does not replace the current max');

  #Test if we try to set max to something bigger
  $mcc->add_feature_type($cs, 'exon', $greater_length);
  is($mcc->fetch_max_length_by_CoordSystem_feature_type($cs, 'exon'), $greater_length, 'Submission of a greater length does replace the current max');
}

{
  #Adding new maximum to a new coord system
  my $cs = $db->get_CoordSystemAdaptor->fetch_by_name('contig');
  my $count = count_rows($db, 'meta_coord');
  my $expected_max = -1;
  # Generate 10 random lengths
  for( my $i=0; $i<10; $i++ ) {
    my $length = int(rand( 1000) + 100);
    $mcc->add_feature_type($cs, 'exon', $length );
    $expected_max = $length if ( $length > $expected_max );
  }
  my $actual_max = $mcc->fetch_max_length_by_CoordSystem_feature_type( $cs, 'exon' );
  is($actual_max, $expected_max, 'The expected and actual maximums agree');
  is(count_rows($db, 'meta_coord'), ($count + 1), 'meta_coord has grown by one row');
}

{
  my @coord_systems = @{$mcc->fetch_all_CoordSystems_by_feature_type('exon')};
  is(scalar(@coord_systems), 2, 'We have two entries returned for exon');
  is_deeply([qw/chromosome contig/], [sort map { $_->name() } @coord_systems], 'We have the two expected coordinate systems returned');
}

$multi->restore('core', 'meta_coord');

done_testing();
