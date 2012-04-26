use strict;
use warnings;

use Test::More tests => 64;

use Bio::EnsEMBL::Test::TestUtils;
use IO::String;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::ProjectionSegment;

our $verbose = 0;

#
# TEST - Slice Compiles
#
ok(1);


my $CHR           = '20';
my $START         = 30_270_000;
my $END           = 31_200_000;
my $STRAND        = 1;
my $SEQ_REGION_LENGTH = 50e6;

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');

#
# TEST - Slice creation from adaptor
#
my $slice_adaptor = $db->get_SliceAdaptor;
my $csa = $db->get_CoordSystemAdaptor();

my $slice = $slice_adaptor->fetch_by_region('chromosome', $CHR, $START, $END);
ok($slice->seq_region_name eq $CHR);
ok($slice->start == $START);
ok($slice->end == $END);
ok($slice->seq_region_length == 62842997);
ok($slice->adaptor == $slice_adaptor);


#
#TEST - Slice::new
#
my $coord_system = $csa->fetch_by_name('chromosome');


my $test_seq = 'ATGCATGCATGCATGCATGCATGC';
my $test_slice = new Bio::EnsEMBL::Slice
  (-seq_region_name  => 'misc',
   -seq_region_length => 24,
   -start            => 1,
   -end              => 24,
   -strand           => 1,
   -coord_system     => $coord_system,
   -seq              => $test_seq,
   );




ok($test_slice->length == 24);

my $hash = $test_slice->get_base_count;
my $a = $hash->{'a'};
my $c = $hash->{'c'};
my $t = $hash->{'t'};
my $g = $hash->{'g'};
my $n = $hash->{'n'};
my $gc_content = $hash->{'%gc'};

ok($a == 6
   && $c == 6 
   && $t == 6 
   && $g == 6 
   && $n == 0 
   && $gc_content == 50 
   && $a+$c+$t+$g+$n == $test_slice->length);


#
# test that subseq works correctly with attached sequence
#

my $subseq = $test_slice->subseq(2, 6);

debug("subseq = $subseq");
ok($subseq eq 'TGCAT');

$subseq = $test_slice->subseq(2,6,-1);
ok($subseq eq 'ATGCA');
debug("subseq = $subseq");

# test that subslice works correctly with attached sequence

my $sub_slice = $test_slice->sub_Slice(2, 6);
ok($sub_slice->seq eq 'TGCAT');

# test that invert works correctly with attached sequence
ok($sub_slice->invert()->seq() eq 'ATGCA');


# test that slice can be created without db, seq or coord system
{
  my $warnings = q{};
  my $new_stderr = IO::String->new(\$warnings);
  my $oldfh = select(STDERR);
  local *STDERR = $new_stderr;
  $test_slice = Bio::EnsEMBL::Slice->new('-seq_region_name' => 'test',
                                         '-start'           => 1,
                                         '-end'             => 3);
  my $check = qr/MSG: Slice without coordinate system/;
  like($warnings, $check, 'Checking we are still warning about lack of coordinate system');
}

ok($test_slice);
ok($test_slice->seq() eq 'NNN');

debug("\$test_slice->name = " . $test_slice->name());

ok($test_slice->name() eq '::test:1:3:1');


$slice = new Bio::EnsEMBL::Slice
  (-seq_region_name   => $CHR,
   -seq_region_length => $SEQ_REGION_LENGTH,
   -start             => $START,
   -end               => $END,
   -strand            => $STRAND,
   -coord_system      => $coord_system);


ok($slice->seq_region_name eq $CHR);

ok($slice->start == $START);
ok($slice->end == $END);
ok($slice->strand == $STRAND);
ok($slice->seq_region_length == $SEQ_REGION_LENGTH);

#
#Test - Slice::adaptor
#
$slice->adaptor($slice_adaptor);
ok($slice->adaptor == $slice_adaptor);


#
#1 Test Slice::name
#
#verify that chr_name start and end are contained in the name
my $name = $slice->name;
ok($name eq "chromosome:NCBI33:$CHR:$START:$END:$STRAND");

#
# Test Slice::length
#
ok($slice->length == ($END-$START + 1));

# Test exception name
is($slice->assembly_exception_type(), 'REF', 'Type of region is REF');

#
# Test get_attributes
#

my $clone = $slice_adaptor->fetch_by_region('clone','AL121583.25');

my @attrib = @{$clone->get_all_Attributes('htg_phase')};

ok(@attrib == 1 && $attrib[0]->value() == 4);

#
# Test expand
#
my $len = $clone->length();

$clone = $clone->expand(100,100);
ok(($clone->start == -99) && ($clone->end() == $len+100));

$clone = $clone->expand(-100,-100);
ok(($clone->start == 1) && ($clone->end() == $len));

$clone = $clone->expand(0,1000);
ok(($clone->start == 1) && ($clone->end() == $len + 1000));

$clone = $clone->expand(-1000, 0);
ok(($clone->start == 1001) && ($clone->end() == $len + 1000));


#
# Test Slice::invert
#
my $inverted_slice = $slice->invert;
ok($slice != $inverted_slice); #slice is not same object as inverted slice
#inverted slice on opposite strand
ok($slice->strand == ($inverted_slice->strand * -1)); 
#slice still on same strand
ok($slice->strand == $STRAND);


#
# Test Slice::seq
#
my $seq = uc $slice->seq;
my $invert_seq = uc $slice->invert->seq;

ok(length($seq) == $slice->length); #sequence is correct length

$seq = reverse $seq;  #reverse complement seq
$seq =~ tr/ACTG/TGAC/; 

ok($seq eq $invert_seq); #revcom same as seq on inverted slice

#
# Test Slice::subseq
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
# Test Slice::get_all_PredictionTranscripts
#
my $pts = $slice->get_all_PredictionTranscripts;
ok(@$pts == 24);

#
# Test Slice::get_seq_region_id
#
ok($slice->get_seq_region_id());

#
# Test Slice::get_all_DnaAlignFeatures
#
my $count = 0;
my $dafs = $slice->get_all_DnaAlignFeatures;
ok(@$dafs == 27081);
$count += scalar @$dafs;

#
# Test Slice::get_all_ProteinAlignFeatures
#
my $pafs = $slice->get_all_ProteinAlignFeatures;
ok(@$pafs == 7205);
$count += scalar @$pafs;

#
# Test Slice::get_all_SimilarityFeatures
#
ok($count == scalar @{$slice->get_all_SimilarityFeatures});

#
#  Test Slice::get_all_SimpleFeatures
#
ok(scalar @{$slice->get_all_SimpleFeatures});

#
#  Test Slice::get_all_RepeatFeatures
#
ok(scalar @{$slice->get_all_RepeatFeatures});

#
#  Test Slice::get_all_Genes
#
ok(scalar @{$slice->get_all_Genes});

#
#  Test Slice::get_all_Genes_by_type
#
ok(scalar @{$slice->get_all_Genes_by_type('protein_coding')});

#
#  Test Slice::get_all_Transcripts
#
ok(scalar @{$slice->get_all_Transcripts});



#
# Test Slice::get_all_KaryotypeBands
#
ok(scalar @{$slice->get_all_KaryotypeBands});


#
# Test Slice::get_RepeatMaskedSeq
#
$seq = $slice->seq;
ok(length($slice->get_repeatmasked_seq->seq) == length($seq));

my $softmasked_seq = $slice->get_repeatmasked_seq(['RepeatMask'], 1)->seq;

ok($softmasked_seq ne $seq);
ok(uc($softmasked_seq) eq $seq);

$softmasked_seq = $seq = undef;

#
# Test Slice::get_all_MiscFeatures
#
ok(scalar @{$slice->get_all_MiscFeatures()});

#
# Test Slice::project
#

my @segments = @{$slice->project( 'seqlevel' )};
ok(scalar @segments );

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


#my $super_slices = $slice->get_all_supercontig_Slices();


##
## get_all_supercontig_Slices()
##
#debug( "Supercontig starts at ".$super_slices->[0]->chr_start() );

#ok( $super_slices->[0]->chr_start() == 29591966 );

#debug( "Supercontig name ".$super_slices->[0]->name() );

#ok( $super_slices->[0]->name() eq "NT_028392" );

#
# get_base_count
#
$hash = $slice->get_base_count;
$a = $hash->{'a'};
$c = $hash->{'c'};
$t = $hash->{'t'};
$g = $hash->{'g'};
$n = $hash->{'n'};
$gc_content = $hash->{'%gc'};

debug( "Base count: a=$a c=$c t=$t g=$g n=$n \%gc=$gc_content");
ok($a == 234371 
   && $c == 224761 
   && $t == 243734 
   && $g == 227135 
   && $n == 0 
   && $gc_content == 48.59 
   && $a+$c+$t+$g+$n == $slice->length);


$slice = $slice_adaptor->fetch_by_region('chromosome', '20', 10, 30);

my $sr_slice = $slice->seq_region_Slice();

ok($sr_slice->start() == 1 &&
   $sr_slice->end()   == $slice->seq_region_length() &&
   $sr_slice->strand() == 1);








# synonym tests
debug("START syn test");
my $multi = $multi_db;
$multi->save("core", "seq_region_synonym");

debug("get slice");
$slice = $slice_adaptor->fetch_by_region('chromosome', 20, 1, 10);

my @alt_names = @{$slice->get_all_synonyms()};
foreach my $syn (@alt_names){
  debug("syn\t".$syn->name."\n");
}
debug("altnames ".scalar(@alt_names)."\n");
ok(@alt_names == 2);


$slice->add_synonym("20ish");
@alt_names = @{$slice->get_all_synonyms()};

ok(@alt_names == 3);


#slcie aleady stored so need to store syns
my $syn_adap =  $db->get_SeqRegionSynonymAdaptor; 
foreach my $syn (@alt_names){
 $syn_adap->store($syn);	
}

$slice = $slice_adaptor->fetch_by_region('chromosome', 20, 1, 10);

@alt_names = @{$slice->get_all_synonyms()};

ok(@alt_names == 3);

$multi->restore();



$multi->save("core", 'seq_region_synonym');

$slice = $slice_adaptor->fetch_by_region('chromosome', 1, 1, 10);

@alt_names = @{$slice->get_all_synonyms()};

ok(@alt_names == 0);

$slice->add_synonym("1ish");

@alt_names = @{$slice->get_all_synonyms()};

ok(@alt_names == 1);

foreach my $syn (@alt_names){
 $syn_adap->store($syn);	
}


$slice = $slice_adaptor->fetch_by_region('chromosome', 1, 1, 10);

@alt_names = @{$slice->get_all_synonyms()};

ok(@alt_names == 1);

$multi->restore();


#Test assembly exception type on HAP
my $hap_slice = $slice_adaptor->fetch_by_region(undef, '20_HAP1');
is($hap_slice->assembly_exception_type(), 'HAP', 'Ensuring haplotype regions are HAP');
my $chr_one_slice = $slice_adaptor->fetch_by_region('chromosome', '1', 1, 10);
is($chr_one_slice->assembly_exception_type(), 'REF', 'Ensuring reference regions are REF');

