use strict;
use warnings;

use Test::More;

use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Slice;

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


my $slice = $slice_adaptor->fetch_by_region('chromosome', $CHR, $START, $END);
compare_compliments($slice, $seq_adaptor);

$slice = $slice_adaptor->fetch_by_region('clone','AL031658.11');
compare_compliments($slice, $seq_adaptor);

$slice = $slice_adaptor->fetch_by_region('supercontig','NT_028392');
compare_compliments($slice, $seq_adaptor);

$slice = $slice_adaptor->fetch_by_region('contig', 'AL031658.11.1.162976');
compare_compliments($slice, $seq_adaptor);


sub compare_compliments {
  my $slice = shift;
  my $seq_adaptor = shift;

  my $seq = ${$seq_adaptor->fetch_by_Slice_start_end_strand($slice,1,undef,1)};

  debug('FORWARD STRAND SLICE SEQ for ' . $slice->name());
  debug($seq);

  my $invert_seq = 
    ${$seq_adaptor->fetch_by_Slice_start_end_strand($slice->invert,1,undef,1)};

  debug('REVERSE STRAND SLICE SEQ for ' . $slice->name());
  debug($invert_seq);

  ok(length($seq) == $slice->length); #sequence is correct length

  $seq = reverse $seq;  #reverse complement seq
  $seq =~ tr/ACTG/TGAC/;

  ok($seq eq $invert_seq); #revcom same as seq on inverted slice
}

done_testing();
