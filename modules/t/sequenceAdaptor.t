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

use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Attribute;

our $verbose= 0;


my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');

my $CHR           = '20';
my $START         = 30_220_000;
my $END           = 31_200_000;
my $STRAND        = 1;

#
# Test fetch_by_Slice_start_end_strand
#

my $slice_adaptor = $db->get_SliceAdaptor;
my $seq_adaptor = $db->get_SequenceAdaptor();

{
  my $slice = $slice_adaptor->fetch_by_region('chromosome', $CHR, $START, $END);
  compare_compliments($slice, $seq_adaptor);

  #Bigger than 1Mb
  $slice = $slice_adaptor->fetch_by_region('chromosome', $CHR, $START, $START+2_000_000);
  compare_compliments($slice, $seq_adaptor);

  $slice = $slice_adaptor->fetch_by_region('clone','AL031658.11');
  compare_compliments($slice, $seq_adaptor);

  $slice = $slice_adaptor->fetch_by_region('supercontig','NT_028392');
  compare_compliments($slice, $seq_adaptor);

  $slice = $slice_adaptor->fetch_by_region('contig', 'AL031658.11.1.162976');
  compare_compliments($slice, $seq_adaptor);
}

$multi_db->save('core','seq_region_attrib');
{
  #Adding an insertion _rna_edit into the seq_region_attrib table. This should
  #not occur and we should fail to build the adaptor whilst it is still there.
  #We will attach it to the clone
  my $attribute_adaptor = $db->get_AttributeAdaptor();
  my $bad_attribute = Bio::EnsEMBL::Attribute->new(-CODE => '_rna_edit', -VALUE => '10 10 AT'); # 1bp edit but string is 2bp
  my $slice = $slice_adaptor->fetch_by_region('clone','AL031658.11');
  $attribute_adaptor->store_on_Slice($slice, [$bad_attribute]);

  throws_ok { $seq_adaptor->new($db); } qr/Edit length .+ substitutions/, 'Exception will be thrown when attempting to build a SequenceAdaptor with a mis-match edit length _rna_edit (only support subs)';
}
$multi_db->restore('core','seq_region_attrib');

## Use patch data to test LRG _rna_edits (it has an LRG, test core doesn't)
{
  my $patch_db    = $multi_db->get_DBAdaptor('patch');
  $slice_adaptor = $patch_db->get_SliceAdaptor();

  # Unedited base at this location is G, ensure it was changed to a T
  my $slice = $slice_adaptor->fetch_by_region(undef, 'LRG_11', 28690, 28690);
  is($slice->seq, 'T', '_rna_edits of LRG produces correct base');

  # Go up stream, check for off-by-one regression
  $slice = $slice_adaptor->fetch_by_region(undef, 'LRG_11', 28689, 28689);
  is($slice->seq, 'C', 'One bp up from LRG edit, emit one base, one base only please');

  # Go up stream, across the _rna_edit, is it correct?
  $slice = $slice_adaptor->fetch_by_region(undef, 'LRG_11', 28689, 28690);
  is($slice->seq, 'CT', 'One bp up from LRG edit and in to the edit');

  # And check downstream from the _rna_edit
  $slice = $slice_adaptor->fetch_by_region(undef, 'LRG_11', 28691, 28691);
  is($slice->seq, 'G', 'Down bp up from LRG edit, emit one base, one base only please');

  # And a multi-bp LRG edit
  $slice = $slice_adaptor->fetch_by_region(undef, 'LRG_11', 28790, 28795);
  is($slice->seq, 'TATCGC', 'Check spanning across an LRG edit');

  # And a multi-bp LRG edit
  $slice = $slice_adaptor->fetch_by_region(undef, 'LRG_11', 28790, 28791);
  is($slice->seq, 'TA', 'Check spanning across an LRG edit');

  # And a multi-bp LRG edit
  $slice = $slice_adaptor->fetch_by_region(undef, 'LRG_11', 28792, 28793);
  is($slice->seq, 'TC', 'Check spanning across an LRG edit');

  # And a multi-bp LRG edit
  $slice = $slice_adaptor->fetch_by_region(undef, 'LRG_11', 28820, 28825);
  is($slice->seq, 'AGTATT', 'Check large edit fully within sequence');

  $slice = $slice_adaptor->fetch_by_region(undef, 'LRG_11', 28822, 28822);
  is($slice->seq, 'T', 'Check sequence fully within a larger edit');

}

sub compare_compliments {
  my $slice = shift;
  my $seq_adaptor = shift;

  my $seq = ${$seq_adaptor->fetch_by_Slice_start_end_strand($slice,1,undef,1)};

  note('FORWARD STRAND SLICE SEQ for ' . $slice->name());
  note($slice->length);

  my $invert_seq = 
    ${$seq_adaptor->fetch_by_Slice_start_end_strand($slice->invert,1,undef,1)};

  note('REVERSE STRAND SLICE SEQ for ' . $slice->name());
  

  is(length($seq), $slice->length, 'sequence is correct length');

  $seq = reverse $seq;  #reverse complement seq
  $seq =~ tr/ACTG/TGAC/;

  #Only use ok here; is would just be crazy
  my $ok = ok($seq eq $invert_seq, 'revcom same as seq on inverted slice');
  if(! $ok ) {
    diag 'SEQ:    '.$seq;
    diag 'INVERT: '.$invert_seq;
  }
}

done_testing();
