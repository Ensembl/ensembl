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

#
#1 TEST - Slice Compiles
#
ok(1); 


my $CHR           = '20';
my $START         = 30_270_000;
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
my $seq = uc $slice->seq;
my $invert_seq = uc $slice->invert->seq;

ok(length($seq) == $slice->length); #sequence is correct length

$seq = reverse $seq;  #reverse complement seq
$seq =~ tr/ACTG/TGAC/; 

ok($seq eq $invert_seq); #revcom same as seq on inverted slice

#
# 25-26 Test Slice::subseq
#
my $SPAN = 10;
my $sub_seq = uc $slice->subseq(-$SPAN,$SPAN);
my $invert_sub_seq = uc $slice->invert->subseq( $slice->length - $SPAN + 1, 
						$slice->length + $SPAN + 1);

ok(length $sub_seq == (2*$SPAN) + 1 ); 
$sub_seq = reverse $sub_seq;
$sub_seq =~ tr/ACTG/TGAC/;

ok($sub_seq eq $invert_sub_seq);

#
# 27 Test Slice::get_all_PredictionTranscripts
#
my $pts = $slice->get_all_PredictionTranscripts;
ok(scalar @$pts);


#
# 28 Test Slice::get_all_DnaAlignFeatures
#
my $count = 0;
my $dafs = $slice->get_all_DnaAlignFeatures;
ok(scalar @$dafs);
$count += scalar @$dafs;

#
# 29 Test Slice::get_all_ProteinAlignFeatures
#
my $pafs = $slice->get_all_ProteinAlignFeatures;
ok(scalar @$pafs);
$count += scalar @$pafs;

#
# 30 Test Slice::get_all_SimilarityFeatures
#
ok($count == scalar @{$slice->get_all_SimilarityFeatures});

#
# 31 Test Slice::get_all_SimpleFeatures
#
ok(scalar @{$slice->get_all_SimpleFeatures});

#
# 32 Test Slice::get_all_RepeatFeatures
#
ok(scalar @{$slice->get_all_RepeatFeatures});

#
# 33 Test Slice::get_all_Genes
#
ok(scalar @{$slice->get_all_Genes});

#
# 34 Test Slice::get_all_Genes_by_type
#
ok(scalar @{$slice->get_all_Genes_by_type('ensembl')});

#
# 35 Test Slice::chr_name
#
my $old_val = $slice->chr_name;
my $new_val = 'Y';
$slice->chr_name($new_val);
ok($slice->chr_name eq $new_val);
$slice->chr_name($old_val);

#
# 36 Test Slice::chr_start
#
$old_val = $slice->chr_start;
$new_val = 123;
$slice->chr_start($new_val);
ok($slice->chr_start == $new_val);
$slice->chr_start($old_val);

#
# 37 Test Slice::chr_end
#
$old_val = $slice->chr_end;
$new_val = 1234567;
$slice->chr_end($new_val);
ok($slice->chr_end == $new_val);
$slice->chr_end($old_val);

#
# 38 Test Slice::strand
#
$old_val = $slice->strand;
$new_val = $old_val * -1;
$slice->strand($new_val);
ok($slice->strand == $new_val);
$slice->strand($old_val);

#
# 39 Test Slice::assembly_type
#
$old_val = $slice->assembly_type;
$new_val = 'TEST';
$slice->assembly_type($new_val);
ok($slice->assembly_type eq $new_val);
$slice->assembly_type($old_val);


#
# 40 Test Slice::get_all_KaryotypeBands
#
ok(scalar @{$slice->get_all_KaryotypeBands});


#
# 41-42 Test Slice::get_Chromosome
#
my $chromo;
ok($chromo = $slice->get_Chromosome);
ok($chromo->chr_name eq $slice->chr_name);

#
# 43-44 Test Slice::get_RepeatMaskedSeq
#
$seq = $slice->seq;
ok(length($slice->get_repeatmasked_seq->seq) == length($seq));

my $softmasked_seq = $slice->get_repeatmasked_seq(['RepeatMask'], 1)->seq;

ok($softmasked_seq ne $seq);
ok(uc($softmasked_seq) eq $seq);

$softmasked_seq = $seq = undef;  

#
# 45 Test Slice::get_all_MapFrags
#
# ok(scalar @{$slice->get_all_MapFrags('cloneset')});

#
# 46 Test Slice::get_tiling_path
#

ok(scalar @{$slice->get_tiling_path});


my $super_slices = $slice->get_all_supercontig_Slices();


#
# 47-48 get_all_supercontig_Slices()
#
debug( "Supercontig starts at ".$super_slices->[0]->chr_start() );

ok( $super_slices->[0]->chr_start() == 29591966 );

debug( "Supercontig name ".$super_slices->[0]->name() );

ok( $super_slices->[0]->name() eq "NT_028392" );

#
# 49 get_base_count
#
my $hash = $slice->get_base_count;
my $a = $hash->{'a'};
my $c = $hash->{'c'};
my $t = $hash->{'t'};
my $g = $hash->{'g'};
my $n = $hash->{'n'};
my $gc_content = $hash->{'%gc'};

debug( "Base count: a=$a c=$c t=$t g=$g n=$n \%gc=$gc_content");
ok($a == 234371 
   && $c == 224761 
   && $t == 243734 
   && $g == 227135 
   && $n == 0 
   && $gc_content == 48.59 
   && $a+$c+$t+$g+$n == $slice->length);

