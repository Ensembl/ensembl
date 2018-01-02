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

use Bio::EnsEMBL::RepeatMaskedSlice;

use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor('core');


my $slice_adaptor = $db->get_SliceAdaptor();

my $slice = $slice_adaptor->fetch_by_region('chromosome', '20');

my $rm_slice = Bio::EnsEMBL::RepeatMaskedSlice->new
  (-START             => 30_980_001,
   -END               => 30_980_100,
   -STRAND            => 1,
   -ADAPTOR           => $slice_adaptor,
   -SEQ_REGION_NAME   => '20',
   -SEQ_REGION_LENGTH => $slice->seq_region_length(),
   -COORD_SYSTEM      => $slice->coord_system(),
   -SOFT_MASK         => 0);


my $seq = $rm_slice->seq();
ok(length($seq) == 100);

debug("rm_slice->seq           = $seq");

my $subseq = $rm_slice->subseq(2,100);
ok(length($subseq) == 99);
debug("rm_slice->subseq(2,100) =  $subseq");
ok($subseq eq substr($seq, 1));


$rm_slice = Bio::EnsEMBL::RepeatMaskedSlice->new
  (-START             => 30_980_001,
   -END               => 30_980_100,
   -STRAND            => 1,
   -ADAPTOR           => $slice_adaptor,
   -SEQ_REGION_NAME   => '20',
   -SEQ_REGION_LENGTH => $slice->seq_region_length(),
   -COORD_SYSTEM      => $slice->coord_system(),
   -SOFT_MASK         => 1);

my $sm_seq = $rm_slice->seq();

my $normal_seq = $slice->subseq(30_980_001,30_980_100);

#check that uppercasing softmasked makes the same sequence as non masked
ok(uc($sm_seq) eq $normal_seq && $sm_seq ne $normal_seq);

#check that replacing lowercase with Ns gives hardmasked
$sm_seq =~ tr/[a-z]/N/;
ok($sm_seq eq $seq);

#
# make sure everything works on negative strand too
#
$rm_slice = Bio::EnsEMBL::RepeatMaskedSlice->new
  (-START             => 30_980_001,
   -END               => 30_980_100,
   -STRAND            => -1,
   -ADAPTOR           => $slice_adaptor,
   -SEQ_REGION_NAME   => '20',
   -SEQ_REGION_LENGTH => $slice->seq_region_length(),
   -COORD_SYSTEM      => $slice->coord_system(),
   -SOFT_MASK         => 0);


$seq = $rm_slice->seq();
debug("(-ve) rm_slice->seq           = $seq");

$subseq = $rm_slice->subseq(2,100);
debug("(-ve) rm_slice->subseq(2,100) =  $subseq");
ok($subseq eq substr($seq, 1));

done_testing();
