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

use Data::Dumper;
use Test::More;
use Test::Warnings;
use Test::Exception;
use Test::Differences;

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::ProjectionSegment;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);

#
# Slice Compiles
#
use_ok('Bio::EnsEMBL::CircularSlice');

my $circular_invert_broken = 1;

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new('circ');
my $db = $multi_db->get_DBAdaptor('core');
my $dbc = $db->dbc;
my $dbname = $dbc->dbname;

#
# Slice creation from adaptor
#
my $CHR = 'Chromosome';
my $START = 5_300_000;
my $END = 10_000;
my $STRAND = -1;
my $SEQ_REGION_LENGTH = 5_495_278;
my $COORD_SYSTEM = 'chromosome';
my $COORD_SYSTEM_VERSION = 'GCA_000292705.1';

my $slice_adaptor = $db->get_SliceAdaptor;
my $csa = $db->get_CoordSystemAdaptor();
my $coord_system = $csa->fetch_by_name($COORD_SYSTEM);

my $slice = $slice_adaptor->fetch_by_region($COORD_SYSTEM, $CHR, $START, $END);
isa_ok($slice, 'Bio::EnsEMBL::CircularSlice');
is($slice->is_circular, 1,"slice is circular");
is($slice->seq_region_name, $CHR,"seq region name $CHR");
is($slice->start, $START, "start == $START");
is($slice->end, $END, "end == $END");
is($slice->seq_region_length, $SEQ_REGION_LENGTH, 
   "seq_region_length == $SEQ_REGION_LENGTH");
is($slice->adaptor, $slice_adaptor, "adaptor is adaptor");

# centrepoint method
is($slice->centrepoint, 5402639, "slice centre point");

# _split private method
my ($sl1, $sl2) = $slice->_split;
isa_ok($sl1, 'Bio::EnsEMBL::CircularSlice');
is($sl1->is_circular, 1, "subslice is circular");
is($sl1->seq_region_name, $CHR,"subslice seq region name $CHR");
is($sl1->start, $START, "subslice start == $START");
is($sl1->end, $SEQ_REGION_LENGTH, "subslice end == $SEQ_REGION_LENGTH");

isa_ok($sl2, 'Bio::EnsEMBL::CircularSlice');
is($sl2->is_circular, 1, "subslice is circular");
is($sl2->seq_region_name, $CHR,"subslice seq region name $CHR");
is($sl2->start, 1, "subslice start == 1");
is($sl2->end, $END, "subslice end == $END");

#
# seq method
#
# - slice spanning the origin of replication
my $seq = $slice->seq;
is(length $seq, $slice->length, "sequence length");

# - slice not spanning the origin of replication
$slice = Bio::EnsEMBL::CircularSlice->new(-seq_region_name   => $CHR,
					  -seq_region_length => $SEQ_REGION_LENGTH,
					  -start             => 10,
					  -end               => 100,
					  -strand            => 1,
					  -coord_system      => $coord_system);
$seq = $slice->seq;
is(length $seq, $slice->length, "sequence length");



my $slstart = $slice_adaptor->fetch_by_region($COORD_SYSTEM, $CHR, 1, $END);
my $slend   = $slice_adaptor->fetch_by_region($COORD_SYSTEM, $CHR, $START, 
					      $SEQ_REGION_LENGTH);

#
# CircularSlice::new
#
my $test_seq = 'ATGCATGCATGCATGCATGCATGC';
my $test_slice = 
  Bio::EnsEMBL::CircularSlice->new(-seq_region_name   => 'misc',
				   -seq_region_length => 24,
				   -start             => 1,
				   -end               => 24,
				   -strand            => 1,
				   -coord_system      => $coord_system,
				   -seq               => $test_seq);
isa_ok($test_slice, 'Bio::EnsEMBL::CircularSlice');
is($test_slice->length, 24, 'slice length 24');
is(length $test_slice->seq, 24, 'slice seq length 24');


#
# Base count
#
my $hash = $test_slice->get_base_count;
my $a = $hash->{'a'};
my $c = $hash->{'c'};
my $t = $hash->{'t'};
my $g = $hash->{'g'};
my $n = $hash->{'n'};
my $gc_content = $hash->{'%gc'};

ok($a == 6 && $c == 6 && $t == 6 && $g == 6 && $n == 0 && $gc_content == 50 && 
   $a+$c+$t+$g+$n == $test_slice->length, 'counts of bases');

#
# - subseq works correctly with attached sequence
# - subslice works correctly with attached sequence
# - subslice works correctly with attached sequence
# - invert works correctly with attached sequence
#
my $subseq = $test_slice->subseq(2, 6);
is($subseq, 'TGCAT', 'subseq works correctly with attached sequence');

$subseq = $test_slice->subseq(2, 6, -1);
is($subseq, 'ATGCA', 'complement of the subseq');

my $sub_slice = $test_slice->sub_Slice(2, 6);
is($sub_slice->seq, 'TGCAT', 'subslice works with attached sequence');
is($sub_slice->invert()->seq(), 'ATGCA',
   'invert works correctly with attached sequence');

#
# Slice can be created without db, seq or coord system
#
my $check = qr/MSG: CircularSlice without coordinate system/;
warns_like{
  $test_slice = Bio::EnsEMBL::CircularSlice->new(-SEQ_REGION_NAME => 'test', -START => 1, -END => 3)
} qr/$check/, 'Checking we are still warning about lack of coordinate system'; 

isa_ok($test_slice, 'Bio::EnsEMBL::CircularSlice');
is($test_slice->seq(), 'NNN','sequence of created slice is NNN');
is($test_slice->name(),'::test:1:3:1','name of the created slice');

#
# - CircularSlice::adaptor
# - CircularSlice::name: chr_name start and end are contained in the name
# - CircularSlice::length 
# - CircularSlice::assembly_exception_type
#
$slice = Bio::EnsEMBL::CircularSlice->new(-seq_region_name   => $CHR,
					  -seq_region_length => $SEQ_REGION_LENGTH,
					  -start             => $START,
					  -end               => $END,
					  -strand            => 1,
					  -coord_system      => $coord_system);

isa_ok($slice, 'Bio::EnsEMBL::CircularSlice');
is($slice->seq_region_name, $CHR, "seq region name is $CHR");
is($slice->start, $START, "slice start $START");
is($slice->end, $END, "slice end $END");
is($slice->strand, 1, "slice strand");
is($slice->seq_region_length, $SEQ_REGION_LENGTH ,
   "slice seq region length $SEQ_REGION_LENGTH");

my $name = $slice->name; 
$slice->adaptor($slice_adaptor);
is($slice->adaptor, $slice_adaptor, 'slice adaptor method');
is($name, "chromosome:$COORD_SYSTEM_VERSION:$CHR:$START:$END:1", 
   '(chr_name, start, end) contained in name');
is($slice->length, 
   ($SEQ_REGION_LENGTH - ($START > $END ? $START - $END : $END-$START) + 1),
   'CircularSlice::length');
is($slice->assembly_exception_type(), 'REF', 'Type of region is REF');

#
# Test get_all_* methods with genes spanning the origin of replication
#
# consider all possible cases of relations between the slice and the
# spanning gene(s)
#
# Legend:
# - 'o': origin of replication
# - '|': gene boundary 
# - '[': slice boundary
#
my $slices = 
  {
   #
   # a. slice spanning the origin of replication:
   #    slice->start > slice->end
   #
   # - seq_region_start >= slice->start
   #
   #   ---[--|---|--o------]---
   #   ---[----|----o---|--]---
   #
   # - seq_region_end >= slice->start
   #
   #   ---[--|---|-o-------]---
   #
   a1 => { 
	  seq_region_id => 11, 
	  seq_region_name => $CHR, 
	  seq_region_length => $SEQ_REGION_LENGTH, 
	  start => 5493000, 
	  end => 400, 
	  coord_system => 'chromosome' },
   #
   # - seq_region_start <= slice->end
   #
   #   ---[-----o----|---|---]---
   #   ---[-----o----|-----]--|--
   #
   # - seq_region_end <= slice->end
   #
   #   ---[---|--o----|-----]---
   #
   a2 => { seq_region_id => 11, 
   	   seq_region_name => $CHR, 
   	   seq_region_length => $SEQ_REGION_LENGTH, 
   	   start => 5493000, 
   	   end => 2500, 
   	   coord_system => 'chromosome' },   
   #
   # - seq_region_start > seq_region_end
   #
   #   ---|---[---o---]---|---
   #
   a3 => { seq_region_id => 11, 
   	   seq_region_name => $CHR, 
   	   seq_region_length => $SEQ_REGION_LENGTH, 
   	   start => 5495000, 
   	   end => 100, 
   	   coord_system => 'chromosome' },
   #
   # b. slice not spanning the origin of replication:
   #    slice->start <= slice->end
   #
   # - seq_region_start <= slice->end AND seq_region_end >= slice->start
   #
   #   ---|--[--|---]--o--------
   #   ---[---|---|-]--o--------
   #   
   b1 => { seq_region_id => 11, 
	   seq_region_name => $CHR, 
   	   seq_region_length => $SEQ_REGION_LENGTH, 
   	   start => 5490000, 
   	   end => 5495000, 
   	   coord_system => 'chromosome' },
   #
   # - seq_region_start > seq_region_end AND seq_region_start <= slice->end
   #
   #   ---[--|---]--o----|----
   #   ---|--[---]--o----|----
   #
   b2 => { seq_region_id => 11, 
	   seq_region_name => $CHR, 
   	   seq_region_length => $SEQ_REGION_LENGTH, 
   	   start => 5494000, 
   	   end => 5495000, 
   	   coord_system => 'chromosome' },
   #
   # - seq_region_start > seq_region_end AND seq_region_end >= slice->start
   #
   #   ---|------o---[---|---]---
   #   ---|------o---[----]--|---
   #
   b2 => { seq_region_id => 11, 
	   seq_region_name => $CHR, 
   	   seq_region_length => $SEQ_REGION_LENGTH, 
   	   start => 100, 
   	   end => 400, 
   	   coord_system => 'chromosome' },
   
};

#
# dinamically generate data for testing various objects overlapping
# each slice
# data for each object (e.g. gene, transcript, exon etc) is queried
# directly from the db and compared (mainly the location definying 
# attributes) with the results of fetching by slice
#
my @query_templates = 
  (
   # slice->start > slice->end
   "SELECT %s_id as db_id, seq_region_start, seq_region_end, seq_region_strand%s FROM %s WHERE seq_region_id = %d AND (seq_region_start >= %d OR seq_region_start <= %d OR seq_region_end >= %d OR seq_region_end <= %d OR seq_region_start > seq_region_end)",
   # slice->start <= slice->end
   "SELECT %s_id as db_id, seq_region_start, seq_region_end, seq_region_strand%s FROM %s WHERE seq_region_id = %d AND ((seq_region_start <= %d AND seq_region_end >= %d) OR (seq_region_start > seq_region_end AND (seq_region_start <= %d OR seq_region_end >= %d)))"
);

# 
# NOTE
# Cannot test with the query above the following features:
# - AssemblyExceptionFeature
# - MarkerFeature (does not have strand attribute)
#
my @tables = qw(gene transcript exon simple_feature); # removed since there are no features of these types and BaseFeatureAdaptor would complain: prediction_transcript dna_align_feature protein_align_feature misc_feature ); 
my @stable_id_tables = qw(gene transcript exon);

# $dbc->do(sprintf "use $dbname");
my $sql_helper = $dbc->sql_helper;
foreach my $sl (sort keys %{$slices}) {
  note "Testing case $sl";
  my $s = $slices->{$sl};

  for my $strand (qw(1 -1)) {
  
    $slice = Bio::EnsEMBL::CircularSlice->new(-seq_region_name   => $s->{seq_region_name},
					      -seq_region_length => $s->{seq_region_length},
					      -start             => $s->{start},
					      -end               => $s->{end},
					      -strand            => $strand,
					      -coord_system      => $csa->fetch_by_name($s->{coord_system}),
					      -adaptor           => $slice_adaptor);

    isa_ok($slice, 'Bio::EnsEMBL::CircularSlice');
    is($slice->seq_region_name, $s->{seq_region_name}, sprintf "seq region name: %s", $s->{seq_region_name});
    is($slice->start, $s->{start}, sprintf "slice start: %d", $s->{start});
    is($slice->end, $s->{end}, sprintf "slice end: %d", $s->{end});
    is($slice->strand, $strand, sprintf "slice strand: %d", $strand);
    is($slice->seq_region_length, $s->{seq_region_length},
       sprintf "slice seq region length: %d", $s->{seq_region_length});

    my ($seq_id, $sstart, $send, $sstrand, $srl) = 
      ($s->{seq_region_id}, $s->{start}, $s->{end}, $strand, $s->{seq_region_length});

    foreach my $table (@tables) {
      my $query;
      
      my $table_in_stable_id_tables = 0;
      foreach my $stable_id_table (@stable_id_tables) {
        $table_in_stable_id_tables = 1 if $stable_id_table eq $table;
        last;
      }

      if ($sstart > $send) {
        $query = sprintf $query_templates[0], $table, $table_in_stable_id_tables ? ', stable_id':'', $table, $seq_id, $sstart, $send, $sstart, $send;
      } 
      else {
        $query = sprintf $query_templates[1], $table, $table_in_stable_id_tables ? ', stable_id':'', $table, $seq_id, $send, $sstart, $send, $sstart; 
      }

      my $expected_objects = $sql_helper->execute(
        -SQL      => $query,
        -USE_HASHREFS => 1,
        -CALLBACK => sub {
          my $row = shift @_;
          return feature_slice_boundaries($row, $sstart, $send, $sstrand, $srl);
        }
      );
    
      my $method = ucfirst($table); $method =~ s/_([a-z])/\u$1/g; $method = sprintf 'get_all_%ss', $method;
      # print "Calling $method\n";
      my $got_objects = $slice->$method;
      is(scalar @{$got_objects}, scalar @{$expected_objects}, sprintf "Number of $table objects on slice [%d, %d]", $sstart, $send);

      foreach my $expected (@{$expected_objects}) {
        my $got = (grep { $expected->{db_id} == $_->dbID } @{$got_objects})[0];
      
        SKIP: {
          skip 'Cannot test attributes with no returned object', 1 unless $got;
          ok($got, sprintf "Object of type $table (%d) retrieved", $got->dbID);	
          is($got->stable_id, $expected->{stable_id}, sprintf "%d: stable id", $got->dbID) if $table_in_stable_id_tables;
          is($got->start, $expected->{start}, sprintf "%d: start", $got->dbID);
          is($got->end, $expected->{end}, sprintf "%d: end", $got->dbID);
          is($got->strand, $expected->{strand}, sprintf "%d: strand", $got->dbID);
	  is($got->seq_region_start, $expected->{seq_region_start}, sprintf "%d: seq_region_start", $got->dbID);
	  is($got->seq_region_end, $expected->{seq_region_end}, sprintf "%d: seq_region_end", $got->dbID);
        }
      }
    }
  }
}

#
# Test retrieval of features (not genes, transcripts, exons)
#
# NOTE
# 
# If we exclude density features, simple features are the only feature 
# type stored in the test DB. 
#
$slice = $slice_adaptor->fetch_by_region($COORD_SYSTEM, $CHR, 5335000, 5337000); # 5336890, 5337706
is($slice->seq_region_name, $CHR,"seq region name $CHR");
is($slice->start, 5335000, "start == 5335000");
is($slice->end, 5337000, "end == 5337000");
is($slice->seq_region_length, $SEQ_REGION_LENGTH, 
   "seq_region_length == $SEQ_REGION_LENGTH");
is($slice->adaptor, $slice_adaptor, "adaptor is adaptor");

my $simple_features = $slice->get_all_SimpleFeatures;
is(scalar @{$simple_features}, 1, "Number of simple features on chromosome:11:5335000-5337000");
SKIP: {
  skip 'Cannot test simple feature location attributes', 1 unless scalar @{$simple_features} == 1;

  my $simple_feature = $simple_features->[0];
  is($simple_feature->start, 1891, "simple feature start");
  is($simple_feature->end, 2707, "simple feature end");
  is($simple_feature->strand, 1, "simple feature strand"); 
}

#
# Test get_attributes
#
my ($start, $end) = (30_000, 40_000);
$slice = Bio::EnsEMBL::CircularSlice->new(-seq_region_name   => $CHR,
					  -seq_region_length => $SEQ_REGION_LENGTH,
					  -start             => $start,
					  -end               => $end,
					  -strand            => 1,
					  -coord_system      => $coord_system,
					  -adaptor           => $slice_adaptor);

my @attrib = @{$slice->get_all_Attributes('circular_seq')};
is(@attrib == 1 && $attrib[0]->value() == 1, 1, 'get_all_Attributes using circular_seq attribute');

#
# Test expand
#
$slice = $slice->expand(100,100);
isa_ok($slice, 'Bio::EnsEMBL::CircularSlice');
is(($slice->start == $start - 100) && ($slice->end() == $end + 100), 1, 'expand: expand 100,100');

$slice = $slice->expand(-100,-100);
is(($slice->start == $start) && ($slice->end() == $end), 1, 'expand: contract -100,-100');

$slice = $slice->expand(0,1000);
is(($slice->start == $start) && ($slice->end() == $end + 1_000), 1, 'expand: extend end 0,1000');

$slice = $slice->expand(-1000, 0);
is(($slice->start == $start + 1_000) && ($slice->end() == $end + 1_000), 1, 'expand: contract start -1000,0');

#
# Test expand across origin
#
my $len = $slice_adaptor->fetch_by_region($COORD_SYSTEM, $CHR)->length;
($start, $end) = ($len - 10_000, 10_000); 
$slice = $slice_adaptor->fetch_by_region($COORD_SYSTEM, $CHR, $start, $end);
isa_ok($slice, 'Bio::EnsEMBL::CircularSlice');

$slice = $slice->expand(100,100);
isa_ok($slice, 'Bio::EnsEMBL::CircularSlice');
is(($slice->start == $start - 100) && ($slice->end() == $end + 100), 1, 'expand across origin: 100,100 ');

$slice = $slice->expand(-100,-100);
is(($slice->start == $start) && ($slice->end() == $end), 1, 'expand across origin: -100,-100');

$slice = $slice->expand(0,1000);
is(($slice->start == $start) && ($slice->end() == $end + 1_000), 1, 'expand across origin: 0,1000');

$slice = $slice->expand(-1000, 0);
is(($slice->start == $start + 1_000) && ($slice->end() == $end + 1_000), 1, 'expand across origin: -1000,0');
#
# Test expand from start and across
#
($start, $end) = (10, 10_000); 
$slice = $slice_adaptor->fetch_by_region($COORD_SYSTEM,$CHR, $start, $end);

$slice = $slice->expand(100,100);
isa_ok($slice, 'Bio::EnsEMBL::Slice');
is(($slice->start == $len - 90) && ($slice->end() == $end + 100),1,'expand near start: 100,100 ');

$slice = $slice->expand(-100,-100);
is(($slice->start == $start) && ($slice->end() == $end),1,'expand near start: -100,-100 ');
#
# Test expand with end at origin
#
($start, $end) = ($len - 10_000, $len); 
$slice = $slice_adaptor->fetch_by_region($COORD_SYSTEM, $CHR, $start, $end);
isa_ok($slice, 'Bio::EnsEMBL::Slice');

$slice = $slice->expand(0,1000);
isa_ok($slice, 'Bio::EnsEMBL::Slice');
is(($slice->start == $start) && ($slice->end() == 1000), 1,'expand near end: 0,1000');

$slice = $slice->expand(0, -1000);
isa_ok($slice, 'Bio::EnsEMBL::Slice');
is(($slice->start == $start) && ($slice->end() == $end), 1,'expand near end: 0,-1000');

#
# Test CircularSlice::invert
#
my $inverted_slice = $slice->invert;
isnt($slice, $inverted_slice,'slice is not same object as inverted slice');
#inverted slice on opposite strand
is($slice->strand,($inverted_slice->strand * -1),'inverted slice strand is correct'); 
#slice still on same strand
is($slice->strand, 1, 'slice still on same strand');


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
#
# $hash = $slice->get_base_count;
# $a = $hash->{'a'};
# $c = $hash->{'c'};
# $t = $hash->{'t'};
# $g = $hash->{'g'};
# $n = $hash->{'n'};
# $gc_content = $hash->{'%gc'};

# note( "Base count: a=$a c=$c t=$t g=$g n=$n \%gc=$gc_content");
# is($a, 10977, 'base count a');
# is($c, 8792, 'base count c');
# is($t, 13113, 'base count t');
# is($g, 9473, 'base count g');
# is($n, 43643, 'base count n');
# is($gc_content, 43.12, 'gc content'); 
# is($a+$c+$t+$g+$n, $slice->length, 'total bases == length');

#
# Test seq_region_Slice
#
$slice = $slice_adaptor->fetch_by_region($COORD_SYSTEM, $CHR, 10, 30);
my $sr_slice = $slice->seq_region_Slice();

is($sr_slice->start(), 1, 'seq_region_Slice: start == 1');
is($sr_slice->end(), $slice->seq_region_length, 'seq_region_Slice: end == length');
is($sr_slice->strand(), 1, 'seq_region_Slice: strand == 1');

#
# Test synonyms
#
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
is(scalar @alt_names, 0, 'number of alt names');

$slice->add_synonym("PrimaryChromosome");
@alt_names = @{$slice->get_all_synonyms()};
is(scalar @alt_names, 1, 'number of alt names after adding one');


# slice aleady stored so need to store syns
my $syn_adap =  $db->get_SeqRegionSynonymAdaptor; 
foreach my $syn (@alt_names){
 $syn_adap->store($syn);	
}

$slice = $slice_adaptor->fetch_by_region($COORD_SYSTEM, $CHR, 1, 10);
@alt_names = @{$slice->get_all_synonyms()};
is(scalar @alt_names, 1, 'number of alt names after saving and reloading');

$multi->restore();
$multi->save("core", 'seq_region_synonym');

$slice = $slice_adaptor->fetch_by_region($COORD_SYSTEM, $CHR, 1, 10);
@alt_names = @{$slice->get_all_synonyms()};
is(scalar @alt_names, 0, 'number of seq region synonyms');

$slice->add_synonym("1ish");
@alt_names = @{$slice->get_all_synonyms()};
is(scalar @alt_names, 1, 'number of seq region synonyms after adding one');

foreach my $syn (@alt_names){
 $syn_adap->store($syn);	
}

$slice = $slice_adaptor->fetch_by_region($COORD_SYSTEM, $CHR, 1, 10);
@alt_names = @{$slice->get_all_synonyms()};

is(scalar @alt_names, 1, 'number of seq region synonyms after saving and reloading');

$multi->restore();

#Test assembly exception type on HAP
my $pla_slice = $slice_adaptor->fetch_by_region($COORD_SYSTEM, $CHR, 1, 10);
is($pla_slice->assembly_exception_type(), 'REF', 'Ensuring reference regions are REF');

done_testing();

sub feature_slice_boundaries {
  my ($attrs, $sstart, $send, $sstrand, $srl) = @_;
  return unless defined $attrs;
  
  die "Undefined slice boundaries"
    unless defined $sstart and defined $send;
  die "Undefined slice strand"
    unless defined $sstrand;
  die "Undefined seq_region_length"
    unless defined $srl;

  my ($srs, $sre) =
    ($attrs->{seq_region_start}, $attrs->{seq_region_end});
  my ($start, $end);

  if ($sstrand == 1) {
    $start = $srs - $sstart + 1;
    $end = $sre - $sstart + 1;

    if ($sstart > $send) {
      if ($srs >= $sstart) {
	if ($end < 0) {
	  $end += $srl;
	}

      } elsif ($srs <= $send) {
	die "start > 0" unless $start < 0;
	die "end > 0" unless $end < 0;

	$start += $srl;
	$end += $srl;

      } elsif ($sre >= $sstart) {
	# do nothing

      } elsif ($sre <= $send) {
	if ($srs < $send) {
	  $start += $srl;
	}
	
	die "end > 0" unless $end < 0;
	$end += $srl;

      } elsif ($srs > $sre) {
	die "srs > sstart or sre < send"
	  unless $srs <= $sstart and $sre >= $send;

	$end += $srl;

      } else {
	die "Shouldn't be here";
      }
      
    } else {

      if ($srs <= $send and $sre >= $sstart) {
	# do nothing
      } elsif ($srs > $sre) {
	if ($srs <= $send) {
	  die "end > 0" unless $end < 0;

	  $end += $srl;

	} elsif ($sre >= $sstart) {

	  $start -= $srl;
	} else {
	  die "Shouldn't be here";
	}
      }
    }

  } else {
    $start = $send - $sre + 1;
    $end = $send - $srs + 1;

    if ($sstart > $send) {
      if ($srs >= $sstart) {
	$end += $srl;

	if ($sre > $sstart) {
	  $start += $srl;
	}

      } elsif ($srs <= $send) {
	# do nothing
      } elsif ($sre >= $sstart) {
	die "start > 0" unless $start < 0;
	die "end > 0" unless $end < 0;

	$start += $srl;
	$end += $srl;

      } elsif ($sre <= $send) {
	die "start < 0" unless $start > 0;

	if ($end < 0) {
	  $end += $srl;
	}

      } elsif ($srs > $sre) {
	die "srs > sstart or sre < send"
	  unless $srs <= $sstart and $sre >= $send;

	$end += $srl;

      } else {
	die "Shouldn't be here";
      }
      
    } else {

      if ($srs <= $send and $sre >= $sstart) {
	# do nothing
      } elsif ($srs > $sre) {
	if ($srs <= $send) {
	  
	  $start -= $srl;

	} elsif ($sre >= $sstart) {
	  $end += $srl;

	} else {
	  die "Shouldn't be here";
	}
      }
    }
    
  }
  
  return {
	  db_id => $attrs->{db_id},
	  stable_id => $attrs->{stable_id},
	  start => $start,
	  end => $end,
	  seq_region_start => $srs,
	  seq_region_end => $sre,
	  strand => $attrs->{seq_region_strand} * $sstrand
	 };

}
