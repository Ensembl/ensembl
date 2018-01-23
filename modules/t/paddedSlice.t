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
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::PaddedSlice;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor('core');
my $csa   = $db->get_CoordSystemAdaptor();

my $coord_system = $csa->fetch_by_name('chromosome');

my $test_seq = 'ATGC'x5;
my $padded_seq = (q{N}x10).$test_seq.(q{N}x20);
my $test_slice = Bio::EnsEMBL::Slice->new(
  -seq_region_name  => 'misc',
  -seq_region_length => 50,
  -start            => 11,
  -end              => 30,
  -strand           => 1,
  -coord_system     => $coord_system,
  -seq              => $test_seq,
);
my $padded = Bio::EnsEMBL::PaddedSlice->new(-SLICE => $test_slice);

#Test full length seq
is($test_slice->seq(),  $test_seq,  'Checking seq comes back as expected');
is($padded->seq(),      $padded_seq,'Checking padded seq comes back as expected');

#Test padded subseqs
is($padded->subseq(1, 5), 'N'x5, 'Checking padded at start');
is($padded->subseq(46, 50), 'N'x5, 'Checking padded at end');
is($padded->subseq(10, 14), 'NATGC', 'Checking overlapping padded at start');
is($padded->subseq(15, 18), 'ATGC', 'Checking mid subseq as expected');
is($padded->subseq(28, 34), 'TGCNNNN', 'Checking overlapping padded at end');

#Checking attributes
is($test_slice->name(), $padded->name(), 'Checking names are the same');

done_testing();
