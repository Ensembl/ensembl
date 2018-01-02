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
use IO::String;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::ProjectionSegment;
use Test::Exception;
use Test::Differences;

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
is($slice->seq_region_name, $CHR, "Slice name is $CHR");
is($slice->start, $START, "Slice start is $START");
is($slice->end, $END, "Slice end is $END");
is($slice->seq_region_length, 62842997, "Slice length is correct");
is($slice->adaptor, $slice_adaptor, "Slice has adaptor $slice_adaptor");


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




is($test_slice->length, 24, 'Tested slice length');

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
is($subseq, 'TGCAT', "Subseq is $subseq");

$subseq = $test_slice->subseq(2,6,-1);
is($subseq, 'ATGCA', "Subseq is $subseq");

# test that subslice works correctly with attached sequence

my $sub_slice = $test_slice->sub_Slice(2, 6);
is($sub_slice->seq, 'TGCAT', "Sub slice seq is correct");

# test that invert works correctly with attached sequence
is($sub_slice->invert()->seq(), 'ATGCA', "Inverted sub slice seq is correct");


# test that slice can be created without db, seq or coord system
{
  my $check = qr/MSG: Slice without coordinate system/;
  warns_like {
  $test_slice = Bio::EnsEMBL::Slice->new('-seq_region_name' => 'test', '-start'           => 1, '-end'             => 3)
  } qr/$check/, 'Checking we are still warning about lack of coordinate system';
}

ok($test_slice);
is($test_slice->seq(), 'NNN', "Test slice seq is only N's");

is($test_slice->name(), '::test:1:3:1', "Slice name is $test_slice->name");

# Test we can make a Slice with a seq_region_name of 0, odd but legal
$slice = new Bio::EnsEMBL::Slice
    (-seq_region_name   => 0,
     -seq_region_length => $SEQ_REGION_LENGTH,
     -start             => $START,
     -end               => $END,
     -strand            => $STRAND,
     -coord_system      => $coord_system);

is($slice->seq_region_name(), 0, 'Create Slice with seq_region_name of 0, odd but legal');

$slice = new Bio::EnsEMBL::Slice
  (-seq_region_name   => $CHR,
   -seq_region_length => $SEQ_REGION_LENGTH,
   -start             => $START,
   -end               => $END,
   -strand            => $STRAND,
   -coord_system      => $coord_system);


is($slice->seq_region_name, $CHR);

is($slice->start, $START, "Slice start is $START");
is($slice->seq_region_start, $START, "Slice seq_region_start is $START");
is($slice->end, $END, "Slice end is $END");
is($slice->seq_region_end, $END, "Slice seq_region_end is $END");
is($slice->strand, $STRAND, "Slice strand is $STRAND");
is($slice->seq_region_length, $SEQ_REGION_LENGTH, "Slice length is $SEQ_REGION_LENGTH");

#
#Test - Slice::adaptor
#
$slice->adaptor($slice_adaptor);
is($slice->adaptor, $slice_adaptor, "Slice has adaptor $slice_adaptor");


#
#1 Test Slice::name
#
#verify that chr_name start and end are contained in the name
my $name = $slice->name;
is($name, "chromosome:NCBI33:$CHR:$START:$END:$STRAND", "Chromosome name is $name");

#
# Test Slice::length
#
is($slice->length, ($END-$START + 1));

# Test exception name
is($slice->assembly_exception_type(), 'REF', 'Type of region is REF');

#
# Test get_attributes
#

my $clone = $slice_adaptor->fetch_by_region('clone','AL121583.25');

my @attrib = @{$clone->get_all_Attributes('htg_phase')};

is(@attrib, 1, "Attrib found");
is($attrib[0]->value(), 4, "First attrib is 4");

#
# Test get_all_DitagFeatures
#

my @ditags = @{ $slice->get_all_DitagFeatures() };

is(scalar(@ditags), 0, "Fetched ditag features from slice");

#
# Test expand
#
my $len = $clone->length();

$clone = $clone->expand(100,100);
is($clone->start, -99, "Clone start is correct");
is($clone->end(), $len+100, "Clone end is correct");

$clone = $clone->expand(-100,-100);
is($clone->start, 1, "Clone start matches 1");
is($clone->end(), $len, "Clone end matches $len");

$clone = $clone->expand(0,1000);
is($clone->start, 1, "Clone start matches 1");
is($clone->end(), $len + 1000, "Clone end matches $len + 1000");

$clone = $clone->expand(-1000, 0, 1);
is($clone->start, 1001, "Expanded clone start is correct if forced");
is($clone->end(), $len + 1000, "Expanded clone end is correct if forced");

#
# Test constrain_to_seq_region
#
my $tidy_clone = $clone->expand(1000000,10000000);
$tidy_clone = $tidy_clone->constrain_to_seq_region;
is($tidy_clone->start, 1, "Tidy clone is correct");
is($tidy_clone->end, 84710, 'Huge expand call truncates nicely');

$tidy_clone = $clone->expand(0,-1000);
$tidy_clone = $tidy_clone->constrain_to_seq_region;
is($tidy_clone->start, 1001, "Tidy clone $tidy_clone->start and $tidy_clone->end are correct");
is($tidy_clone->end(), 84710, 'constrain_to_seq_region does no harm');


#
# Test Slice::invert
#
my $inverted_slice = $slice->invert;
ok($slice != $inverted_slice); #slice is not same object as inverted slice
#inverted slice on opposite strand
is($slice->strand, ($inverted_slice->strand * -1), "Inverted slice on opposite strand is identical to initial slice"); 
#slice still on same strand
is($slice->strand, $STRAND, "Slice is still on the same strand");


#
# Test Slice::seq
#
my $seq = uc $slice->seq;
my $invert_seq = uc $slice->invert->seq;

is(length($seq), $slice->length, "Sequence is correct length");

$seq = reverse $seq;  #reverse complement seq
$seq =~ tr/ACTG/TGAC/; 

eq_or_diff($seq, $invert_seq, "revcom same as seq on inverted slice");

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

eq_or_diff($sub_seq, $invert_sub_seq);

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
is(scalar(@$dafs), 27081, 'Checking count of returned DnaAlignFeatures');
$count += scalar @$dafs;

#
# Test Slice::get_all_ProteinAlignFeatures
#
my $pafs = $slice->get_all_ProteinAlignFeatures;
is(scalar(@$pafs),7205, 'Checking count of returned ProteinAlignFeatures');
$count += scalar @$pafs;

#
# Test Slice::get_all_SimilarityFeatures
#
is($count, scalar @{$slice->get_all_SimilarityFeatures}, "Checking count of returned similarity features");

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
#  Test Slice::get_all_Genes_by_source
#
ok(scalar @{$slice->get_all_Genes_by_source('ensembl')});

#
#  Test Slice::get_all_Transcripts
#
is(scalar @{$slice->get_all_Transcripts}, 23, "Found 23 transcripts");

#
#  Test Slice::get_all_Transcripts_by_type
#
is(scalar @{$slice->get_all_Transcripts_by_type('protein_coding')}, 23, "Found 23 protein_coding transcripts");

#
#  Test Slice::get_all_Transcripts_by_source
#
is(scalar @{$slice->get_all_Transcripts_by_source('ensembl')}, 20, "Found 20 ensembl transcripts");

#
# Test Slice:get_all_Exons
#

my @exons = @{$slice->get_all_Exons};
is(scalar(@exons), 155, "Fetched all exons");


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
is(@alt_names, 2, "Got 2 altnames");


$slice->add_synonym("20ish");
@alt_names = @{$slice->get_all_synonyms()};

is(@alt_names, 3, "Got 3 alt names");

is($alt_names[1]->name, 'anoth_20', 'Correctly retrieved synonym name');
is($alt_names[1]->dbname, 'RFAM', 'Correctly retrieved synonym external db');


#slcie aleady stored so need to store syns
my $syn_adap =  $db->get_SeqRegionSynonymAdaptor; 
foreach my $syn (@alt_names){
 $syn_adap->store($syn);	
}

$slice = $slice_adaptor->fetch_by_region('chromosome', 20, 1, 10);

@alt_names = @{$slice->get_all_synonyms()};

is(@alt_names, 3, "Got 3 altnames");

$multi->restore();



$multi->save("core", 'seq_region_synonym');

$slice = $slice_adaptor->fetch_by_region('chromosome', 1, 1, 10);

@alt_names = @{$slice->get_all_synonyms()};

is(@alt_names, 0, "No altnames returned");


$slice->add_synonym("1ish");

@alt_names = @{$slice->get_all_synonyms()};

is(@alt_names, 1, "One synonym retrieved");

foreach my $syn (@alt_names){
 $syn_adap->store($syn);	
}


$slice = $slice_adaptor->fetch_by_region('chromosome', 1, 1, 10);

@alt_names = @{$slice->get_all_synonyms()};

is(@alt_names, 1, "One synonym found");

$multi->restore();

# Testing synonym searching
{
  my $chr_20 = $slice_adaptor->fetch_by_region('chromosome', 20);
  my ($syn) = @{$chr_20->get_all_synonyms('RFAM')};
  is($syn->name(), 'anoth_20', 'We have the right synonym');
  dies_ok { $chr_20->get_all_synonyms('RFAM', 'wibble') } 'Bad external DB version means dying code';
  dies_ok { $chr_20->get_all_synonyms('RFAMing', 'wibble') } 'Bad external DB name means dying code';
  ($syn) = @{$chr_20->get_all_synonyms('RFAM', 1)};
  is($syn->name(), 'anoth_20', 'We have the right synonym');
}

# test fetch_all on synonym adaptor
my $all_synonyms = $syn_adap->fetch_all();
is(@$all_synonyms, 3, 'fetch_all on synonym adaptor');


#Test assembly exception type on HAP
my $hap_slice = $slice_adaptor->fetch_by_region(undef, '20_HAP1');
is($hap_slice->assembly_exception_type(), 'HAP', 'Ensuring haplotype regions are HAP');
my $chr_one_slice = $slice_adaptor->fetch_by_region('chromosome', '1', 1, 10);
is($chr_one_slice->assembly_exception_type(), 'REF', 'Ensuring reference regions are REF');

#Test slice iterator
{
  my $large_slice = $slice_adaptor->fetch_by_region('chromosome', 1, 1, 21);
  my $map = sub { $_->length() };
  my $si = sub {
    my ($chunk) = @_;
    return $large_slice->sub_Slice_Iterator($chunk)->map($map)->to_arrayref();
  };
  is_deeply($si->(100), [21], 'Subslice larger than actual slice gives just 1 slice back');
  is_deeply($si->(10), [10,10,1], 'Subslice smaller than actual slice gives 3 slices back');
  is_deeply($si->(20), [20,1], 'Subslice just smaller than actual slice gives 2 slices back');
  is_deeply($si->(21), [21], 'Subslice equal to slice size gives 1 slice back');
  my $slice_count = $large_slice->sub_Slice_Iterator(1)->reduce(sub { $_[0]+1 }, 0);
  is($slice_count, 21, 'Giving a subslice size of 1 means 21 slices');
  
  {
    my $fake_slice = Bio::EnsEMBL::Slice->new(-SEQ => 'AAACCCTTTGGGA', 
      -START => 1, -END => 13, -SEQ_REGION_NAME => 'fake', 
      -COORD_SYSTEM => $coord_system);
    my $subseqs = $fake_slice->sub_Slice_Iterator(3)->map(sub { $_->seq() })->to_arrayref();
    my $expected = ['AAA','CCC','TTT','GGG','A'];
    is_deeply($subseqs, $expected, 'Calling seq on subslices returns only the sequence for the bounds');
  }
  
  {
    my $one_bp_slice = Bio::EnsEMBL::Slice->new(-SEQ => 'A', 
      -START => 1, -END => 1, -SEQ_REGION_NAME => 'fake', 
      -COORD_SYSTEM => $coord_system);
    my $subseqs = $one_bp_slice->sub_Slice_Iterator(1)->map(sub { $_->seq() })->to_arrayref();
    my $expected = ['A'];
    is_deeply($subseqs, $expected, 'Calling seq on subslices for 1bp slice returns only an A');
  }
}

# Test alternative region mappings from fetch_by_seq_region_id()
{
  my $current_slice = $slice_adaptor->fetch_by_region('chromosome', $CHR, $START, $END);
  my $alternative_seq_region_id = 1;

  #Get the alternative seq region id. Expect nothing
  ok(! defined $slice_adaptor->fetch_by_seq_region_id($alternative_seq_region_id), 'Asking for a non-existent ID results in no slice returned');

  #Save the tables and add the mapping values in
  my @tables = ('seq_region_mapping', 'mapping_set');
  $multi_db->save('core', @tables);
  my $ms_sql = 'insert into mapping_set (mapping_set_id, internal_schema_build, external_schema_build) values (?,?,?)';
  $db->dbc->sql_helper->execute_update(-SQL => $ms_sql, -PARAMS => [1, $db->_get_schema_build(), 'oldbuild']);
  my $srm_sql = 'insert into seq_region_mapping (mapping_set_id, internal_seq_region_id, external_seq_region_id) values (?,?,?)';
  $db->dbc->sql_helper->execute_update(-SQL => $srm_sql, -PARAMS => [1, $current_slice->get_seq_region_id(), $alternative_seq_region_id]);


  #Force a refresh in CSA
  delete $db->get_CoordSystemAdaptor->{$_} for qw/_internal_seq_region_mapping _external_seq_region_mapping/;
  $db->get_CoordSystemAdaptor->_cache_seq_region_mapping();
  my $alternative_slice = $slice_adaptor->fetch_by_seq_region_id($alternative_seq_region_id);
  ok(!defined $alternative_slice, 'Cannot retrieve the alternative slice without asking to look at alternatives');
  $alternative_slice = $slice_adaptor->fetch_by_seq_region_id($alternative_seq_region_id, undef, undef, undef, 1); #don't care about start,end,strand
  ok(defined $alternative_slice, 'Got a slice after asking for it');
  cmp_ok($current_slice->get_seq_region_id(), '==', $alternative_slice->get_seq_region_id(), 'Seq region IDs should be equivalent even though query seq_region_id was different');
  
  #Restore & force a refresh
  $multi_db->restore('core', @tables);
  delete $db->get_CoordSystemAdaptor->{$_} for qw/_internal_seq_region_mapping _external_seq_region_mapping/;
  $db->get_CoordSystemAdaptor->_cache_seq_region_mapping();

  # Just checking we are back to normal
  $alternative_slice = $slice_adaptor->fetch_by_seq_region_id($alternative_seq_region_id);
  ok(!defined $alternative_slice, 'Cannot retrieve the alternative slice post restore');
}


# Test slice attributes
my $current_slice = $slice_adaptor->fetch_by_region('chromosome', $CHR, $START, $END);
is($current_slice->is_chromosome, 1, "Slice is a chromosome");
is($current_slice->has_karyotype, 1, "Slice has a karyotype attribute");
is($current_slice->karyotype_rank, 20, "Karyotype rank is 20 could be found");

#
# Test get_genome_component
#
debug("Testing fetch_all_by_genome_component");
my $multi_polyploid = Bio::EnsEMBL::Test::MultiTestDB->new("polyploidy");
my $wheatdb = $multi_polyploid->get_DBAdaptor("core");
my $wheat_slice_adaptor = Bio::EnsEMBL::DBSQL::SliceAdaptor->new($wheatdb);
isa_ok($wheat_slice_adaptor, 'Bio::EnsEMBL::DBSQL::SliceAdaptor');

# should get an empty result for a slice on human chr (not polyploidy)
$slice = $slice_adaptor->fetch_by_region('chromosome', $CHR, $START, $END);
isa_ok($slice, 'Bio::EnsEMBL::Slice');
my $genome_component = $slice->get_genome_component();
ok(!$genome_component, "Genome component for human slice");

# test with slices on polyploid genome
$slice = $wheat_slice_adaptor->fetch_by_region('scaffold', 'IWGSC_CSS_5AL_scaff_2697823', 100, 10000);
isa_ok($slice, 'Bio::EnsEMBL::Slice');
$genome_component = $slice->get_genome_component();
is($genome_component, 'A', "Genome component from slice");

$slice = $wheat_slice_adaptor->fetch_by_region('scaffold', 'IWGSC_CSS_6BS_scaff_233977', 1000, 5000, -1);
isa_ok($slice, 'Bio::EnsEMBL::Slice');
$genome_component = $slice->get_genome_component();
is($genome_component, 'B', "Genome component from slice");

$slice = $wheat_slice_adaptor->fetch_by_region('scaffold', 'IWGSC_CSS_6DS_scaff_2121653', 1000, 3000);
isa_ok($slice, 'Bio::EnsEMBL::Slice');
$genome_component = $slice->get_genome_component($slice);
is($genome_component, 'D', "Genome component from slice");

# Deal with rare user case of requesting a region with end = ''

dies_ok(sub { $slice_adaptor->fetch_by_region('chromosome', 20, 1, '') },'Stringy null region end causes an exception');

  



done_testing();

