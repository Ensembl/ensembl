use strict;
use warnings;

use lib 't';

BEGIN { $| = 1;  
	use Test;
	plan tests => 49;
}

use TestUtils qw( debug );

use MultiTestDB;
use Bio::EnsEMBL::Slice;

our $verbose= 0;


my $multi_db = MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');


my $CHR           = '20';
my $START         = 30_270_000;
my $END           = 31_200_000;
my $STRAND        = 1;


#
# Test fetch_by_Slice_start_end_strand
#

my $slice_adaptor = $db->get_SliceAdaptor;
my $slice = $slice_adaptor->fetch_by_chr_start_end($CHR, $START, $END);


my $seq_adaptor = $db->get_SequenceAdaptor();

my $seq = $seq_adaptor->fetch_by_Slice_start_end_strand($slice,1,-1,1);

debug('FORWARD STRAND SLICE SEQ');
debug($seq);

my $invert_seq = 
  $seq_adaptor->fetch_by_Slice_start_end_strand($slice->invert,1,-1,1);

debug('REVERSE STRAND SLICE SEQ');
debug($invert_seq);

ok(length($seq) == $slice->length); #sequence is correct length

$seq = reverse $seq;  #reverse complement seq
$seq =~ tr/ACTG/TGAC/;

ok($seq eq $invert_seq); #revcom same as seq on inverted slice


