use strict;
use warnings;

use lib 't';

BEGIN { $| = 1;  
	use Test;
	plan tests => 27;
}


use MultiTestDB;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Slice;


#
#1 TEST - Slice Compiles
#
ok(1); 


my $CHR           = '20';
my $START         = 31_000_000;
my $END           = 31_200_000;
my $STRAND        = 1;
my $ASSEMBLY_TYPE = 'NCBI_30';
my $DBID          = 123;

my $multi_db = MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');


#
#2-5 TEST - Slice creation from adaptor
#
my $slice_adaptor = $db->get_SliceAdaptor;
my $slice = $slice_adaptor->fetch_by_chr_start_end($CHR, $START, $END);
ok($slice->chr_name eq $CHR);
ok($slice->chr_start == $START); 
ok($slice->chr_end == $END);
ok($slice->adaptor);
  

#
#6 TEST - Slice::new (empty)
#
$slice = new Bio::EnsEMBL::Slice(-empty => 1);
ok($slice);


#
#7-12 TEST - Slice::new
#
$slice = new Bio::EnsEMBL::Slice(-chr_name  => $CHR,
		   -chr_start => $START,
		   -chr_end   => $END,
		   -strand    => $STRAND,
		   -assembly_type => $ASSEMBLY_TYPE,
		   -dbid     => $DBID);



ok($slice->chr_name eq $CHR);
ok($slice->chr_start == $START);
ok($slice->chr_end == $END);
ok($slice->strand == $STRAND);
ok($slice->assembly_type eq $ASSEMBLY_TYPE);
ok($slice->dbID == $DBID);

#
#13 Test - Slice::adaptor
#
$slice->adaptor($slice_adaptor);
ok($slice->adaptor == $slice_adaptor);

#
#14 Test - Slice::dbID
#
$slice->dbID(10);
ok($slice->dbID==10);

#
#15-17 Test Slice::name
#
#verify that chr_name start and end are contained in the name
my $name = $slice->name;
ok($name =~/$CHR/);
ok($name =~/$START/);
ok($name =~/$END/);


#
#18 Test Slice::id
#
ok($slice->id eq $slice->name);


#
#19 Test Slice::length
#
ok($slice->length == ($END-$START + 1));


#
#20-22 Test Slice::invert
#
my $inverted_slice = $slice->invert;
ok($slice != $inverted_slice); #slice is not same object as inverted slice
#inverted slice on opposite strand
ok($slice->strand == ($inverted_slice->strand * -1)); 
#slice still on same strand
ok($slice->strand == $STRAND);


#
# 23-24 Test Slice::seq
#
my $seq = $slice->seq;
my $invert_seq = $slice->invert->seq;

print STDERR "SEQ=[$seq]\n";

ok(length($seq) == $slice->length); #sequence is correct length
print STDERR "[".length($seq)."] != [".$slice->length."]\n"; 
$seq = uc reverse $seq;  #reverse complement seq
$seq =~ s/ACTG/TGAC/g; 
ok($seq eq $invert_seq); #revcom same as seq on inverted slice

#
# 25-26 Test Slice::subseq
#
my $SPAN = 10;
my $sub_seq = $slice->subseq(-$SPAN,$SPAN);
my $invert_sub_seq = $slice->invert->subseq($slice->length + $SPAN, 
					    $slice->length - $SPAN);
ok(length $sub_seq == (2*$SPAN) + 1 ); 
$seq = uc reverse $seq;
$seq =~ s/ACTG/TGAC/g;
ok($seq eq $invert_seq);








