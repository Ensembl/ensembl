use strict;
use warnings;

use Test::More;

use Bio::EnsEMBL::Test::TestUtils qw/capture_std_streams/;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::CircularSlice;
use Bio::EnsEMBL::ProjectionSegment;

my $circular_invert_broken = 1;

#
# TEST - Slice Compiles
#
ok(1);


my $CHR           = '20';
my $START         = 62_800_000;
my $END           = 43_000;
my $STRAND        = 1;
my $SEQ_REGION_LENGTH = 62_842_997;
my $COORD_SYSTEM  = 'chromosome';

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new('circ');
my $db = $multi_db->get_DBAdaptor('core');

#
# TEST - Slice creation from adaptor
#
my $slice_adaptor = $db->get_SliceAdaptor;
my $csa = $db->get_CoordSystemAdaptor();

my $slice = $slice_adaptor->fetch_by_region($COORD_SYSTEM, $CHR, $START, $END);
my $slend = $slice_adaptor->fetch_by_region($COORD_SYSTEM, $CHR, $START, $SEQ_REGION_LENGTH);
my $slstart = $slice_adaptor->fetch_by_region($COORD_SYSTEM, $CHR, 1, $END);
is($slice->is_circular,1,"slice is circular");
is($slice->seq_region_name,$CHR,"seq region name $CHR");
is($slice->start, $START, "start ==$START");
is($slice->end, $END, "end == $END");
is($slice->seq_region_length, $SEQ_REGION_LENGTH, "length == $SEQ_REGION_LENGTH");
is($slice->adaptor, $slice_adaptor, "adaptor is adaptor");

#
#TEST - CircularSlice::new
#
my $coord_system = $csa->fetch_by_name($COORD_SYSTEM);


my $test_seq = 'ATGCATGCATGCATGCATGCATGC';
my $test_slice = new Bio::EnsEMBL::CircularSlice
  (-seq_region_name  => 'misc',
   -seq_region_length => 24,
   -start            => 1,
   -end              => 24,
   -strand           => 1,
   -coord_system     => $coord_system,
   -seq              => $test_seq,
   );




is($test_slice->length, 24, 'slice length 24');

my $hash = $test_slice->get_base_count;
my $a = $hash->{'a'};
my $c = $hash->{'c'};
my $t = $hash->{'t'};
my $g = $hash->{'g'};
my $n = $hash->{'n'};
my $gc_content = $hash->{'%gc'};

is($a == 6
   && $c == 6 
   && $t == 6 
   && $g == 6 
   && $n == 0 
   && $gc_content == 50 
   && $a+$c+$t+$g+$n == $test_slice->length,1,'counts of bases are correct');


#
# test that subseq works correctly with attached sequence
#

my $subseq = $test_slice->subseq(2, 6);

note("subseq = $subseq");
is($subseq,'TGCAT','subseq works correctly with attached sequence');

$subseq = $test_slice->subseq(2,6,-1);
is($subseq, 'ATGCA','complement of the subseq');
note("subseq = $subseq");

# test that subslice works correctly with attached sequence

my $sub_slice = $test_slice->sub_Slice(2, 6);
is($sub_slice->seq,'TGCAT', 'subslice works correctly with attached sequence');

# test that invert works correctly with attached sequence
is($sub_slice->invert()->seq(),'ATGCA','invert works correctly with attached sequence');


# test that slice can be created without db, seq or coord system
capture_std_streams(sub {
  my ($std_out_ref, $std_err_ref) = @_;
  $test_slice = Bio::EnsEMBL::CircularSlice->new(-SEQ_REGION_NAME => 'test', -START => 1, -END => 3);
  my $check = qr/MSG: CircularSlice without coordinate system/;
  like(${$std_err_ref}, $check, 'Checking we are still warning about lack of coordinate system');
  return;
});

is(ref $test_slice,'Bio::EnsEMBL::CircularSlice', 'slice can be created without db, seq or coord system');
is($test_slice->seq(),'NNN','sequence of created slice is NNN');

note("\$test_slice->name = " . $test_slice->name());

is($test_slice->name(),'::test:1:3:1','name of the created slice');


$slice = new Bio::EnsEMBL::CircularSlice
  (-seq_region_name   => $CHR,
   -seq_region_length => $SEQ_REGION_LENGTH,
   -start             => $START,
   -end               => $END,
   -strand            => $STRAND,
   -coord_system      => $coord_system);


is($slice->seq_region_name, $CHR, "seq region name is $CHR");

is($slice->start,             $START             ,"slice start $START");
is($slice->end,               $END               ,"slice end $END");
is($slice->strand,            $STRAND            ,"slice strand $STRAND");
is($slice->seq_region_length, $SEQ_REGION_LENGTH ,"slice seq region length $SEQ_REGION_LENGTH");

#
#Test - CircularSlice::adaptor
#
$slice->adaptor($slice_adaptor);
is($slice->adaptor, $slice_adaptor, 'slice adaptor convenience function');


#
#1 Test CircularSlice::name
#
#verify that chr_name start and end are contained in the name
my $name = $slice->name;
is($name,"chromosome:NCBI33:$CHR:$START:$END:$STRAND", 'chr_name start and end are contained in the name');

#
# Test CircularSlice::length
#
is($slice->length,($SEQ_REGION_LENGTH - ($START > $END ? $START - $END : $END-$START) + 1),'Slice::length');

# Test exception name
is($slice->assembly_exception_type(), 'REF', 'Type of region is REF');

#
# Test get_attributes
#

my ($pstart, $pend) = (30_000,40_000); 
my $plasmid = $slice_adaptor->fetch_by_region($COORD_SYSTEM,$CHR,$pstart,$pend);

my @attrib = @{$plasmid->get_all_Attributes('circular_seq')};
is(@attrib == 1 && $attrib[0]->value() == 1,1,'get_all_Attributes works using circular_seq attribute');

#
# Test expand
#

$plasmid = $plasmid->expand(100,100);
is(($plasmid->start == $pstart - 100) && ($plasmid->end() == $pend + 100), 1, 'expand: expand 100,100');

$plasmid = $plasmid->expand(-100,-100);
is(($plasmid->start == $pstart) && ($plasmid->end() == $pend), 1, 'expand: contract -100,-100');

$plasmid = $plasmid->expand(0,1000);
is(($plasmid->start == $pstart) && ($plasmid->end() == $pend + 1_000), 1, 'expand: extend end 0,1000');

$plasmid = $plasmid->expand(-1000, 0);
is(($plasmid->start == $pstart + 1_000) && ($plasmid->end() == $pend + 1_000), 1, 'expand: contract start -1000,0');
#
# Test expand across origin
#
my $len = $slice_adaptor->fetch_by_region($COORD_SYSTEM,$CHR)->length;
($pstart, $pend) = ($len - 10_000, 10_000); 
$plasmid = $slice_adaptor->fetch_by_region($COORD_SYSTEM,$CHR,$pstart,$pend);

$plasmid = $plasmid->expand(100,100);
is(($plasmid->start == $pstart - 100) && ($plasmid->end() == $pend + 100), 1, 'expand across origin: 100,100 ');

$plasmid = $plasmid->expand(-100,-100);
is(($plasmid->start == $pstart) && ($plasmid->end() == $pend), 1, 'expand across origin: -100,-100');

$plasmid = $plasmid->expand(0,1000);
is(($plasmid->start == $pstart) && ($plasmid->end() == $pend + 1_000), 1, 'expand across origin: 0,1000');

$plasmid = $plasmid->expand(-1000, 0);
is(($plasmid->start == $pstart + 1_000) && ($plasmid->end() == $pend + 1_000), 1, 'expand across origin: -1000,0');
#
# Test expand from start and across
#
($pstart, $pend) = (10, 10_000); 
$plasmid = $slice_adaptor->fetch_by_region($COORD_SYSTEM,$CHR,$pstart,$pend);

$plasmid = $plasmid->expand(100,100);
is(($plasmid->start == $len - 90) && ($plasmid->end() == $pend + 100),1,'expand near start: 100,100 ');

$plasmid = $plasmid->expand(-100,-100);
is(($plasmid->start == $pstart) && ($plasmid->end() == $pend),1,'expand near start: -100,-100 ');
#
# Test expand with end at origin
#
($pstart, $pend) = ($len - 10_000, $len); 
$plasmid = $slice_adaptor->fetch_by_region($COORD_SYSTEM,$CHR,$pstart,$pend);

$plasmid = $plasmid->expand(0,1000);
is(($plasmid->start == $pstart) && ($plasmid->end() == 1000),1,'expand near end: 0,1000');

$plasmid = $plasmid->expand(0, -1000);
is(($plasmid->start == $pstart) && ($plasmid->end() == $pend),1,'expand near end: 0,-1000');

#
# Test CircularSlice::invert
#
my $inverted_slice = $slice->invert;
isnt($slice, $inverted_slice,'slice is not same object as inverted slice');
#inverted slice on opposite strand
is($slice->strand,($inverted_slice->strand * -1),'inverted slice strand is correct'); 
#slice still on same strand
is($slice->strand,$STRAND, 'slice still on same strand');


#
# Test CircularSlice::seq
#
SKIP: {
  skip 'Cannot test invert due to known bug', 2 if $circular_invert_broken;
  
  my $seq = uc $slice->seq;
  my $invert_seq = uc $slice->invert->seq;
  
  is(length($seq), $slice->length,'invert: sequence is correct length');
  
  $seq = reverse $seq;  #reverse complement seq
  $seq =~ tr/ACTG/TGAC/; 
  
  is($seq, $invert_seq, 'invert: revcom same as seq on inverted slice');
}

#
# Test CircularSlice::subseq
#
SKIP: {
  skip 'Cannot test invert due to known bug', 2 if $circular_invert_broken;
  my $SPAN = 10;
  my $sub_seq = uc $slice->subseq(-$SPAN,$SPAN);
  my $invert_sub_seq = uc $slice->invert->subseq( $slice->length - $SPAN + 1, 
            $slice->length + $SPAN + 1);
  
  is(length $sub_seq, $slice->length + (2 * $SPAN), 'invert: sub seq (-10,10) has correct length'); 
  $sub_seq = reverse $sub_seq;
  $sub_seq =~ tr/ACTG/TGAC/;
  
  is($sub_seq, $invert_sub_seq, 'invert: manually inverted sub seq and sub seq inverted by function are equivalent');
}

#
# Test CircularSlice::get_all_PredictionTranscripts
# jh15 20120704: unfortunately we have none at this time in e_coli_k12
my $pts = $slice->get_all_PredictionTranscripts;
is(scalar @$pts, 0, 'number of prediction transcripts');

#
# Test CircularSlice::get_seq_region_id
#
is($slice->get_seq_region_id(),469283,'seq_region_id');

#
# Test CircularSlice::get_all_DnaAlignFeatures
# jh15 20120704: unfortunately we have none at this time in e_coli_k12
my $count = 0;
my $dafs = $slice->get_all_DnaAlignFeatures;
is(scalar @$dafs, 107, 'number of dna align features');
$count += scalar @$dafs;

#
# Test CircularSlice::get_all_ProteinAlignFeatures
# jh15 20120704: unfortunately we have none at this time in e_coli_k12
my $pafs = $slice->get_all_ProteinAlignFeatures;
is(scalar @$pafs, 5, 'number of protein align features');
$count += scalar @$pafs;

#
# Test CircularSlice::get_all_SimilarityFeatures
#
is($count, scalar @{$slice->get_all_SimilarityFeatures}, 'all simililarity features');

#
#  Test CircularSlice::get_all_SimpleFeatures
#
is(scalar @{$slice->get_all_SimpleFeatures},0, 'all simple features');

#
#  Test CircularSlice::get_all_RepeatFeatures
#
is(scalar @{$slice->get_all_RepeatFeatures},0, 'all repeat features');

#
#  Test CircularSlice::get_all_Genes
#
my $total = scalar @{$slstart->get_all_Genes};
$total += scalar @{$slend->get_all_Genes};
is($total, scalar @{$slice->get_all_Genes}, 'spanning origin: all genes');

#
#  Test CircularSlice::get_all_Genes_by_type
#
$total = scalar @{$slstart->get_all_Genes_by_type('protein_coding')};
$total += scalar @{$slend->get_all_Genes_by_type('protein_coding')};
is($total, scalar @{$slice->get_all_Genes_by_type('protein_coding')},'spanning origin: all genes by type');

#
#  Test CircularSlice::get_all_Transcripts
#
$total = scalar @{$slstart->get_all_Transcripts};
$total += scalar @{$slend->get_all_Transcripts};
is($total, scalar @{$slice->get_all_Transcripts},'spanning origin: all transcripts');



#
# Test CircularSlice::get_all_KaryotypeBands
#
 $total = scalar @{$slstart->get_all_KaryotypeBands};
 $total += scalar @{$slend->get_all_KaryotypeBands};
is($total, scalar @{$slice->get_all_KaryotypeBands},'spanning origin: all karyotype bands');


#
# Test Slice::get_RepeatMaskedSeq
#
is(length($slice->get_repeatmasked_seq->seq), length($slice->seq),'repeatmask: length is correct');

my $softmasked_seq = $slice->get_repeatmasked_seq(['RepeatMask'], 1)->seq;

SKIP : {
  skip 'Cannot test due to known invert bug', 2 if $circular_invert_broken;
  isnt($softmasked_seq, $slice->seq,'repeatmask: masked sequence is different');
  is(uc($softmasked_seq), $slice->seq, 'repeatmask: manually reverted soft mask is same as original seq');
}

#
# Test Slice::get_all_MiscFeatures
#
ok(scalar @{$slice->get_all_MiscFeatures()});

#
# Test CircularSlice::get_all_MiscFeatures
#
$total = scalar @{$slstart->get_all_MiscFeatures};
$total += scalar @{$slend->get_all_MiscFeatures};
is($total, scalar @{$slice->get_all_MiscFeatures}, 'across origin: all misc features');

#
# Test CircularSlice::project
#

my @segments = @{$slice->project( 'seqlevel' )};
is(scalar @segments,1, 'project: number of segments' );

eval {
  my @sub_slices = map { $_->to_Slice() } @segments;
  my @starts = map { $_->from_start() } @segments;
  my @ends = map { $_->from_end() } @segments;
};

if( $@ ) {
  debug( "to_Slice call failed on segment of projection" );
  ok(0);
} else {
  ok(1)
}

#
# get_base_count
# ! not reimplemented in CircularSlice
$hash = $slice->get_base_count;
$a = $hash->{'a'};
$c = $hash->{'c'};
$t = $hash->{'t'};
$g = $hash->{'g'};
$n = $hash->{'n'};
$gc_content = $hash->{'%gc'};

note( "Base count: a=$a c=$c t=$t g=$g n=$n \%gc=$gc_content");
is($a, 10977, 'base count a');
is($c, 8792, 'base count c');
is($t, 13113, 'base count t');
is($g, 9473, 'base count g');
is($n, 43643, 'base count n');
is($gc_content, 43.12, 'gc content'); 
is($a+$c+$t+$g+$n, $slice->length, 'total bases == length');

$slice = $slice_adaptor->fetch_by_region($COORD_SYSTEM, $CHR, 10, 30);

my $sr_slice = $slice->seq_region_Slice();

is($sr_slice->start(), 1, 'start == 1');
is($sr_slice->end(), $slice->seq_region_length, 'end == length');
is($sr_slice->strand(), 1, 'strand == 1');

# synonym tests
note("START syn test");
my $multi = $multi_db;
$multi->save("core", "seq_region_synonym");

note("get slice");
$slice = $slice_adaptor->fetch_by_region($COORD_SYSTEM, $CHR, 1, 10);

my @alt_names = @{$slice->get_all_synonyms()};
foreach my $syn (@alt_names){
  note("syn\t".$syn->name."\n");
}
note("altnames ".scalar(@alt_names)."\n");
is(scalar @alt_names, 2, 'number of alt names');


$slice->add_synonym("PrimaryChromosome");
@alt_names = @{$slice->get_all_synonyms()};

is(scalar @alt_names, 3, 'number of alt names after adding one');


#slcie aleady stored so need to store syns
my $syn_adap =  $db->get_SeqRegionSynonymAdaptor; 
foreach my $syn (@alt_names){
 $syn_adap->store($syn);	
}

$slice = $slice_adaptor->fetch_by_region($COORD_SYSTEM, $CHR, 1, 10);

@alt_names = @{$slice->get_all_synonyms()};

is(scalar @alt_names, 3, 'number of alt names after saving and reloading');

$multi->restore();



$multi->save("core", 'seq_region_synonym');

$slice = $slice_adaptor->fetch_by_region($COORD_SYSTEM, $CHR, 1, 10);

@alt_names = @{$slice->get_all_synonyms()};

is(scalar @alt_names, 2, 'number of seq region synonyms');

$slice->add_synonym("1ish");

@alt_names = @{$slice->get_all_synonyms()};

is(scalar @alt_names, 3, 'number of seq region synonyms after adding one');

foreach my $syn (@alt_names){
 $syn_adap->store($syn);	
}


$slice = $slice_adaptor->fetch_by_region($COORD_SYSTEM, $CHR, 1, 10);

@alt_names = @{$slice->get_all_synonyms()};

is(scalar @alt_names, 3, 'number of seq region synonyms after saving and reloading');

$multi->restore();


#Test assembly exception type on HAP
my $pla_slice = $slice_adaptor->fetch_by_region($COORD_SYSTEM, $CHR, 1, 10);
is($pla_slice->assembly_exception_type(), 'REF', 'Ensuring reference regions are REF');

done_testing();
