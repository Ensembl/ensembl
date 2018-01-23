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

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Test::TestUtils;
use Test::Exception;

our $verbose = 0;

my ($CHR, $START, $END, $FLANKING) = ("20", 30_252_000, 31_252_001, 1000);

#
# slice adaptor compiles
#
ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db    = $multi->get_DBAdaptor('core');


#
# SliceAdaptor::new
#
my $slice_adaptor = Bio::EnsEMBL::DBSQL::SliceAdaptor->new($db);
ok($slice_adaptor->isa('Bio::EnsEMBL::DBSQL::SliceAdaptor'));
ok($slice_adaptor->db);

#
# fetch_by_region
#
my $slice = $slice_adaptor->fetch_by_region('chromosome',$CHR, $START, $END);
ok($slice->seq_region_name eq $CHR);
ok($slice->start == $START);
ok($slice->end   == $END);
ok($slice->seq_region_length == 62842997);
debug("slice seq_region length = " . $slice->seq_region_length());


#
# fetch_by_contig_name
#

my $projection = $slice->project('seqlevel');

#it is important to get a contig not cut off by slice start or end
unless(@$projection > 2) {
  warn("There aren't enough tiles in this path for this test to work");
}
my ($chr_start,$chr_end,$contig) = @{$projection->[1]};

ok($contig->length == ($chr_end - $chr_start + 1));

my $seq1 = $slice->subseq($chr_start, $chr_end);
my $seq2 = $contig->seq();

ok($seq1 eq $seq2);


#
# 16-17 fetch by transcript_stable_id
#
my $t_stable_id = 'ENST00000217315';
$slice = $slice_adaptor->fetch_by_transcript_stable_id($t_stable_id);
my $new_slice = $slice_adaptor->fetch_by_transcript_stable_id($t_stable_id,
                                                           $FLANKING);

ok($new_slice->start == $slice->start - $FLANKING);
ok($new_slice->end   == $slice->end   + $FLANKING);


#
# fetch_by_exon_stable_id
#

my $e_stable_id = 'ENSE00001048794';
$slice = $slice_adaptor->fetch_by_exon_stable_id($e_stable_id);
$new_slice = $slice_adaptor->fetch_by_exon_stable_id($e_stable_id, $FLANKING);
ok($new_slice->start == $slice->start - $FLANKING);
ok($new_slice->end   == $slice->end   + $FLANKING);


#
# 18-19 fetch by transcript_id
#
my $transcript = $db->get_TranscriptAdaptor->fetch_by_stable_id($t_stable_id);
my $tid = $transcript->dbID;
$slice = $slice_adaptor->fetch_by_transcript_id($tid);
$new_slice = $slice_adaptor->fetch_by_transcript_id($tid, $FLANKING);
ok($new_slice->start == $slice->start - $FLANKING);
ok($new_slice->end   == $slice->end   + $FLANKING);
ok($slice->seq_region_length == 62842997);
debug("new slice seq_region length = " . $new_slice->seq_region_length());

#
# 20-23 fetch_by_gene_stable_id
#
my $g_stable_id = 'ENSG00000125964';
$slice = $slice_adaptor->fetch_by_gene_stable_id($g_stable_id);
$new_slice = $slice_adaptor->fetch_by_gene_stable_id($g_stable_id, $FLANKING);
ok($new_slice->start == $slice->start - $FLANKING);
ok($new_slice->end   == $slice->end   + $FLANKING);

#verify we can retrieve the gene from this slice
my $gene_found = 0;
foreach my $g (@{$slice->get_all_Genes}) {
  if($g->stable_id eq $g->stable_id) {
    $gene_found = 1;
    last;
  }
}
ok($gene_found);

# same test for flanking slice
$gene_found = 0;
foreach my $g (@{$new_slice->get_all_Genes}) {
  if($g->stable_id eq $g->stable_id) {
    $gene_found = 1;
    last;
  }
}
ok($gene_found);



#
#  fetch_by_region (entire region)
#
$slice = $slice_adaptor->fetch_by_region('chromosome',$CHR);
ok($slice->seq_region_name eq $CHR);
ok($slice->start == 1);


#
# fetch_by_misc_feature_attribute
#
my $flanking= 1000;
$slice = $slice_adaptor->fetch_by_misc_feature_attribute('superctg',
                                                         'NT_030871',
                                                         $flanking);

ok($slice->seq_region_name eq '20');
ok($slice->start == 59707812 - $flanking);
ok($slice->end   == 60855021 + $flanking);

# 
# fetch_by_misc_feature_set
#
$slice = $slice_adaptor->fetch_by_misc_feature_set('superctg', 'NT_030871', 'tilepath');

ok($slice->seq_region_name eq '20');
ok($slice->start == 59707812);
ok($slice->end   == 60855021);


#
# normalized projected slice
#

#
# a slice with a PAR region
# 24,25
#
$slice = $slice_adaptor->fetch_by_region( "chromosome", "Y", 9_000_000, 11_000_000, 1 );

my $results = $slice_adaptor->fetch_normalized_slice_projection( $slice );

debug( "Pseudo autosomal region results" );
for my $projection ( @$results ) {
  debug( "Start: ".$projection->[0] );
  debug( "End: ".$projection->[1] );
  debug( "Slice ".$projection->[2] );
  debug( "-----------" );
}

ok( @$results == 3 );
ok( $results->[1]->[2]->seq_region_name() eq "20" );

#
# a slice with a haplotype 
# 26,27
#

$slice =  $slice_adaptor->fetch_by_region( "chromosome", "20_HAP1", 30_000_000, 31_000_000, 1 );
$results = $slice_adaptor->fetch_normalized_slice_projection( $slice );

debug( "Haplotype projection results" ); 
for my $projection ( @$results ) {
  debug( "Start: ".$projection->[0] );
  debug( "End: ".$projection->[1] );
  debug( "Slice ".$projection->[2] );
  debug( "-----------" );
}

ok( @$results == 3 );
ok( $results->[0]->[2]->seq_region_name() eq "20" );
ok( $results->[1]->[2]->seq_region_name() eq "20_HAP1" );
ok( $results->[2]->[2]->seq_region_name() eq "20" );


#try a projection from chromosome 20 to supercontigs
$slice = $slice_adaptor->fetch_by_region('chromosome', "20", 29_252_000, 
                                         31_252_001 );

debug("Projection from chromosome 20 to supercontig");
my @projection = @{$slice->project('supercontig')};
ok(@projection == 1);
ok($projection[0]->[2]->seq_region_name eq 'NT_028392');
foreach my $seg (@projection) {
  my ($start, $end, $slice) = @$seg;
  debug("$start-$end " . $slice->seq_region_name);
}

#try a projection from clone to supercontig
$slice = $slice_adaptor->fetch_by_region('clone', 'AL121583.25');

debug("Projection from clone AL121583.25 to supercontig");

@projection = @{$slice->project('supercontig')};
ok(@projection == 1);
ok($projection[0]->[2]->seq_region_name eq 'NT_028392');
foreach my $seg (@projection) {
  my ($start, $end, $slice) = @$seg;
  debug("$start-$end -> " . $slice->start . '-'. $slice->end . ' ' . $slice->seq_region_name);
}

#
# test storing a couple of different slices
#
my $csa = $db->get_CoordSystemAdaptor();
my $ctg_cs  = $csa->fetch_by_name('contig');

$multi->save('core', 'seq_region', 'dna', 'assembly');

my $ctg_len = 50;
my $name = 'testregion';

#
# Store a slice with sequence
#

my $ctg_slice = Bio::EnsEMBL::Slice->new(-COORD_SYSTEM    => $ctg_cs,
                                         -SEQ_REGION_NAME => $name,
                                         -SEQ_REGION_LENGTH => $ctg_len,
                                         -START           => 1,
                                         -END             => $ctg_len,
                                         -STRAND          => 1); 

my $seq   = 'A' x $ctg_len;



$slice_adaptor->store($ctg_slice, \$seq);

$ctg_slice = $slice_adaptor->fetch_by_region('contig', $name);

ok($ctg_slice->length == $ctg_len);
ok($ctg_slice->seq eq $seq);
ok($ctg_slice->seq_region_name eq $name);

#
# Store a slice without sequence
#

my $chr_cs  = $csa->fetch_by_name('chromosome');

my $chr_len = 50e6;
$name = 'testregion2';
my $chr_slice = Bio::EnsEMBL::Slice->new(-COORD_SYSTEM    => $chr_cs,
                                         -SEQ_REGION_NAME => $name,
                                         -SEQ_REGION_LENGTH => $chr_len,
                                         -START           => 1,
                                         -END             => $chr_len,
                                         -STRAND          => 1); 

$slice_adaptor->store($chr_slice);

$chr_slice = $slice_adaptor->fetch_by_region('chromosome', $name);
ok($chr_slice->length() == $chr_len);
ok($chr_slice->seq_region_length() == $chr_len);
ok($chr_slice->seq_region_name eq $name);

#
# Update a slice
#

$chr_slice->add_synonym('testregion3');
$slice_adaptor->update($chr_slice);

my $updated_slice = $slice_adaptor->fetch_by_region('chromosome', 'testregion3');
ok($updated_slice->length() == $chr_len);
ok($updated_slice->seq_region_length() == $chr_len);
ok($updated_slice->seq_region_name eq $name);



#
# Store an assembly between the slices
# Retrieve and remove said assembly
#
my $asm_start = 9999;
my $asm_slice = $chr_slice->sub_Slice( $asm_start, $asm_start + $ctg_len - 1 );
my $str = $slice_adaptor->store_assembly( $asm_slice, $ctg_slice );

ok( $str eq "chromosome:NCBI33:testregion2:9999:10048:1<>".
            "contig::testregion:1:50:1" );

my $new_mapper = $slice_adaptor->fetch_assembly($asm_slice, $ctg_slice);
is($new_mapper->from->name, 'chromosome:NCBI33:testregion2:9999:10048:1', 'Mapping from chromosome');
is($new_mapper->to->name, 'contig::testregion:1:50:1', 'Mapping to contig');

$slice_adaptor->remove_assembly($asm_slice, $ctg_slice);
my $empty_mapper = $slice_adaptor->fetch_assembly($asm_slice, $ctg_slice);
is($empty_mapper, undef, 'Mapper has been removed');

my $ctg_map = $chr_slice->project( $ctg_cs->name, $ctg_cs->version );
# Test currently fails as assembly cached somewhere.
#ok( @$ctg_map == 1 and
#    $ctg_map->[0]->[0] == $asm_slice->start and
#    $ctg_map->[0]->[1] == $asm_slice->end and
#    $ctg_map->[0]->[2]->name eq $ctg_slice->name );

my $chr_map = $ctg_slice->project( $chr_cs->name, $chr_cs->version );
# Test currently fails as assembly cached somewhere.
#ok( @$chr_map == 1 and
#    $chr_map->[0]->[0] == $ctg_slice->start and
#    $chr_map->[0]->[1] == $ctg_slice->end and
#    $chr_map->[0]->[2]->name eq $chr_slice->name );


$multi->restore('core', 'seq_region', 'dna', 'assembly');


#
# There was a bug such that features were not being retrieved
# from slices that had a start < 1.  This is a test for that case.
#
$slice = $slice_adaptor->fetch_by_region('chromosome', '20', 1,35_000_000);
debug("slice start = " . $slice->start);
debug("slice end   = " . $slice->end);

my $sfs1 = $slice->get_all_SimpleFeatures();
print_features($sfs1);

$slice = $slice_adaptor->fetch_by_region('chromosome', '20', -10, 35_000_000);

debug("slice start = " . $slice->start);
debug("slice end   = " . $slice->end);

my $sfs2 = $slice->get_all_SimpleFeatures();
print_features($sfs2);

ok(@$sfs1 == @$sfs2);


#
# test fetch_by_name
#
$slice = $slice_adaptor->fetch_by_name($slice->name());

ok($slice->coord_system->name eq 'chromosome');
ok($slice->seq_region_name eq '20');
ok($slice->start == -10);
ok($slice->strand == 1);
ok($slice->end == 35e6);

$slice = $slice_adaptor->fetch_by_name('clone::AL121583.25:1:10000:-1');

ok($slice->coord_system->name eq 'clone');
ok($slice->seq_region_name eq 'AL121583.25');
ok($slice->start == 1);
ok($slice->end == 10000);
ok($slice->strand == -1);


#
# test fetch_all
#

#default no duplicates and reference only
my $slices = $slice_adaptor->fetch_all('chromosome',undef);
print_slices($slices);
is(@$slices, 63, 'References slices for coord system chromosome');

# include duplicates
$slices = $slice_adaptor->fetch_all('chromosome', undef,0, 1);

print_slices($slices);
is(@$slices, 62, 'References slices for coord system chromosome when including duplicates (Y should become 1 region not 2)');


$slices = $slice_adaptor->fetch_all('contig', undef);

ok(@$slices == 13);

print_slices($slices);


$slices = $slice_adaptor->fetch_all('toplevel');

ok(@$slices == 1 && $slices->[0]->seq_region_name() eq '20');
print_slices($slices);

#
# test fetch_all_by_genome_component
#
debug("Testing fetch_all_by_genome_component");
my $multi_polyploid = Bio::EnsEMBL::Test::MultiTestDB->new("polyploidy");
my $wheatdb = $multi_polyploid->get_DBAdaptor("core");
my $wheat_slice_adaptor = Bio::EnsEMBL::DBSQL::SliceAdaptor->new($wheatdb);
isa_ok($wheat_slice_adaptor, 'Bio::EnsEMBL::DBSQL::SliceAdaptor');

# should throw if argument is not provided
throws_ok { $wheat_slice_adaptor->fetch_all_by_genome_component }
  qr/Undefined/, 'Call without argument';

# should throw if argument is an invalid genome component
throws_ok { $wheat_slice_adaptor->fetch_all_by_genome_component('C') }
  qr/Invalid/, 'Call with invalid argument';

# test with valid genome components
debug("Getting top level slices on A");
$slices = $wheat_slice_adaptor->fetch_all_by_genome_component('A');
ok(scalar @$slices == 1, "Number of top level slices on component A");
$slice = $slices->[0];
isa_ok($slice, 'Bio::EnsEMBL::Slice');
is($slice->seq_region_name, "IWGSC_CSS_5AL_scaff_2697823", "seq region name");
is($slice->start, 1, "slice start");
is($slice->end, 5428, "slice end");
is($slice->strand, 1, "slice strand");
is($slice->seq_region_length, 5428, "slice seq region length");
is($slice->adaptor, $wheat_slice_adaptor, 'slice adaptor');

debug("Getting top level slices on B");
$slices = $wheat_slice_adaptor->fetch_all_by_genome_component('B');
ok(scalar @$slices == 1, "Number of top level slices on component B");
$slice = $slices->[0];
isa_ok($slice, 'Bio::EnsEMBL::Slice');
is($slice->seq_region_name, "IWGSC_CSS_6BS_scaff_233977", "seq region name");
is($slice->start, 1, "slice start");
is($slice->end, 4562, "slice end");
is($slice->strand, 1, "slice strand");
is($slice->seq_region_length, 4562, "slice seq region length");
is($slice->adaptor, $wheat_slice_adaptor, 'slice adaptor');

debug("Getting top level slices on D");
$slices = $wheat_slice_adaptor->fetch_all_by_genome_component('D');
ok(scalar @$slices == 1, "Number of top level slices on component D");
$slice = $slices->[0];
isa_ok($slice, 'Bio::EnsEMBL::Slice');
is($slice->seq_region_name, "IWGSC_CSS_6DS_scaff_2121653", "seq region name");
is($slice->start, 1, "slice start");
is($slice->end, 18301, "slice end");
is($slice->strand, 1, "slice strand");
is($slice->seq_region_length, 18301, "slice seq region length");
is($slice->adaptor, $wheat_slice_adaptor, 'slice adaptor');

# 
# test get_genome_component_for_slice
#
debug("Testing get_genome_component_for_slice");

# should throw if argument is not provided
throws_ok { $slice_adaptor->get_genome_component_for_slice }
  qr/Undefined/, 'Call with undefined argument';

# should get an empty result for a slice on human chr (not polyploidy)
$slice = $slice_adaptor->fetch_by_region('chromosome',$CHR, $START, $END);
isa_ok($slice, 'Bio::EnsEMBL::Slice');

my $genome_component = $slice_adaptor->get_genome_component_for_slice($slice);
ok(!$genome_component, "Genome component for human slice");

# test with slices on polyploid genome
$slice = $wheat_slice_adaptor->fetch_by_region('scaffold', 'IWGSC_CSS_5AL_scaff_2697823', 100, 10000);
isa_ok($slice, 'Bio::EnsEMBL::Slice');
$genome_component = $wheat_slice_adaptor->get_genome_component_for_slice($slice);
is($genome_component, 'A', "Genome component from slice");

$slice = $wheat_slice_adaptor->fetch_by_region('scaffold', 'IWGSC_CSS_6BS_scaff_233977', 1000, 5000, -1);
isa_ok($slice, 'Bio::EnsEMBL::Slice');
$genome_component = $wheat_slice_adaptor->get_genome_component_for_slice($slice);
is($genome_component, 'B', "Genome component from slice");

$slice = $wheat_slice_adaptor->fetch_by_region('scaffold', 'IWGSC_CSS_6DS_scaff_2121653', 1000, 3000);
isa_ok($slice, 'Bio::EnsEMBL::Slice');
$genome_component = $wheat_slice_adaptor->get_genome_component_for_slice($slice);
is($genome_component, 'D', "Genome component from slice");

#
# test the fuzzy matching of clone accessions
#
my $clone_name = 'AL031658';
$slice = $slice_adaptor->fetch_by_region('clone', $clone_name);

debug("Fuzzy matched clone name $clone_name Got " . 
     $slice->seq_region_name);

ok($slice->seq_region_name =~ /$clone_name\.\d+/);

#make sure that it does not fuzzy match too much
$slice = $slice_adaptor->fetch_by_region('contig', $clone_name);
ok(!defined($slice));
print_slices([$slice]);

#make sure that you can fetch a seq_region without knowing its version
$slice = $slice_adaptor->fetch_by_region(undef, '20');
ok(defined($slice) && $slice->seq_region_name eq '20');

$slice = $slice_adaptor->fetch_by_region('toplevel', '20');
ok(defined($slice) && $slice->seq_region_name eq '20');

$slice = $slice_adaptor->fetch_by_region('toplevel', '20', 10, 20);
ok(defined($slice) && $slice->start == 10 && $slice->end == 20);

$slice = $slice_adaptor->fetch_by_region(undef, '20', 10, 20, 1, 'NCBI33');
ok(defined($slice) && $slice->seq_region_name eq '20');

$slice = $slice_adaptor->fetch_by_region(undef, '20', 10, 20, 1, 'bogus');
ok(!defined($slice));


$slice = $slice_adaptor->fetch_by_region('toplevel', '20', 10, 20, 1, 'bogus');
ok(defined($slice) && $slice->seq_region_name eq '20');

# try fuzzy matching in conjunction with coord system guessing
$clone_name = 'AL031658';
$slice = $slice_adaptor->fetch_by_region(undef, $clone_name);
ok($slice->seq_region_name =~ /$clone_name\.\d+/);

# Testing synonym fetching
{
  my $syn_slice = $slice_adaptor->fetch_by_region(undef, 'anoth_20');
  is($syn_slice->seq_region_name(), '20', 'Ensuring slice is Chr20 as expected');
  my $chr_syn_slice = $slice_adaptor->fetch_by_region('chromosome', 'anoth_20');
  is($chr_syn_slice->seq_region_name(), '20', 'Ensuring slice is Chr20 as expected');
  $chr_syn_slice = $slice_adaptor->fetch_by_region('toplevel', 'chrx');
  is($chr_syn_slice->seq_region_name(), 'X', 'Ensuring slice is ChrX as expected');
}

#{
#  my @slices = @{$slice_adaptor->fetch_all_by_synonym('anoth_20', 'UniGene')};
#  is(scalar(@slices), 0, 'Checking querying with a bad external name means no Slices');
#  
#  @slices = @{$slice_adaptor->fetch_all_by_synonym('anoth_20', 'RFAM')}; #Yeah ... RFAM
#  is(scalar(@slices), 1, 'Checking querying with a good external name means Slices');
#  is($slices[0]->seq_region_name(), '20', 'Ensuring slice is Chr20 as expected');
#}

# test that with multiple sequence regions with the same name, the
# highest (lowest-numbered) ranked comes out first
$multi->hide('core', 'seq_region');

my $sth = $db->dbc->prepare(qq{INSERT INTO seq_region (coord_system_id, name,
                                                  length)
                SELECT cs.coord_system_id, 'TESTREGION', 1000000
                FROM coord_system cs
                WHERE cs.name in ('supercontig', 'chromosome')});

$sth->execute();
$sth->finish();

$slice = $slice_adaptor->fetch_by_region('toplevel', 'TESTREGION');

ok($slice->seq_region_name() eq 'TESTREGION');
ok($slice->coord_system()->name() eq 'chromosome');


$multi->restore('core', 'seq_region');

###### FETCH BY LOCATION
test_toplevel_location('1:1-1000', 'chromosome', '1', 1, 1000);
test_toplevel_location('1:1-', 'chromosome', '1', 1, 246874334);
test_toplevel_location('1:-10', 'chromosome', '1', 1, 10);
test_toplevel_location('1:100', 'chromosome', '1', 100, 246874334);
test_toplevel_location('1:', 'chromosome', '1', 1, 246874334);
test_toplevel_location('1', 'chromosome', '1', 1, 246874334);

test_toplevel_location('1:1..1000', 'chromosome', '1', 1, 1000);
test_toplevel_location('1:1..', 'chromosome', '1', 1, 246874334);
test_toplevel_location('1:..10', 'chromosome', '1', 1, 10);
test_toplevel_location('1:100', 'chromosome', '1', 100, 246874334);
test_toplevel_location('1:', 'chromosome', '1', 1, 246874334);
test_toplevel_location('1', 'chromosome', '1', 1, 246874334);

test_toplevel_location('1:1_1000', 'chromosome', '1', 1, 1000);
test_toplevel_location('1:1_', 'chromosome', '1', 1, 246874334);
test_toplevel_location('1:_10', 'chromosome', '1', 1, 10);
test_toplevel_location('1:100', 'chromosome', '1', 100, 246874334);
test_toplevel_location('1:', 'chromosome', '1', 1, 246874334);
test_toplevel_location('1', 'chromosome', '1', 1, 246874334);

test_toplevel_location('1: 1-1,000', 'chromosome', '1', 1, 1000);
test_toplevel_location('1: 1-1,000,000', 'chromosome', '1', 1, 1000000);
test_toplevel_location('1: 1-1 000 000', 'chromosome', '1', 1, 1000000);
test_toplevel_location('1: 1', 'chromosome', '1', 1, 246874334);
test_toplevel_location('1: -10', 'chromosome', '1', 1, 10);
test_toplevel_location('1: 100', 'chromosome', '1', 100, 246874334);
test_toplevel_location('1:100..2000000000', 'chromosome', '1', 100, 246874334);
test_toplevel_location('1:100..2E9', 'chromosome', '1', 100, 246874334);

# Try with optional g. addition
test_toplevel_location('1:g.1_1000', 'chromosome', '1', 1, 1000);

# Try chr
my $ucsc = 1;
test_toplevel_location('chr1: 1-1,000', 'chromosome', '1', 1, 1000, 1, $ucsc);
test_toplevel_location('chr1: 1-1,000,000', 'chromosome', '1', 1, 1000000, 1, $ucsc);
test_toplevel_location('chr1: 1-1 000 000', 'chromosome', '1', 1, 1000000, 1, $ucsc);
test_toplevel_location('chr1: 1', 'chromosome', '1', 1, 246874334, 1, $ucsc);
test_toplevel_location('chr1: -10', 'chromosome', '1', 1, 10, 1, $ucsc);
test_toplevel_location('chr1: 100', 'chromosome', '1', 100, 246874334, 1, $ucsc);
test_toplevel_location('chr1:100..2000000000', 'chromosome', '1', 100, 246874334, 1, $ucsc);
test_toplevel_location('chr1:100..2E9', 'chromosome', '1', 100, 246874334, 1, $ucsc);

# Try negative locations
test_toplevel_location('chr1: -10-1,000', 'chromosome', '1', 1, 1000, 1, $ucsc);
test_toplevel_location('chr1: -10..1,000', 'chromosome', '1', 1, 1000, 1, $ucsc);
test_toplevel_location('chr1: -10_1,000', 'chromosome', '1', 1, 1000, 1, $ucsc);
ok(!defined $slice_adaptor->fetch_by_toplevel_location('1:-1000--10', 1), 'Checking with a bogus region with negative coords returns undef');

#Try strands
test_toplevel_location('1:1-1000:1', 'chromosome', '1', 1, 1000, 1);
test_toplevel_location('1:1-1000:-1', 'chromosome', '1', 1, 1000, -1);
test_toplevel_location('1:1-1000:+', 'chromosome', '1', 1, 1000, 1);
test_toplevel_location('1:1-1000:-', 'chromosome', '1', 1, 1000, -1);
test_toplevel_location('1:1-1000..1', 'chromosome', '1', 1, 1000, 1);
test_toplevel_location('1:1-1000--1', 'chromosome', '1', 1, 1000, -1);

dies_ok { $slice_adaptor->fetch_by_toplevel_location(); } 'Checking calling without a location fails';
dies_ok { $slice_adaptor->fetch_by_toplevel_location('', 1); } 'Checking calling with a blank location fails';
dies_ok { $slice_adaptor->fetch_by_toplevel_location('1:1000000000..100', 1); } 'Checking calling with an excessive start throws an error';
ok(!defined $slice_adaptor->fetch_by_toplevel_location('wibble', 1), 'Checking with a bogus region returns undef');
ok(!defined $slice_adaptor->fetch_by_toplevel_location('1:-100--50', 1), 'Checking with a bogus region with negative coords returns undef');

# Try without toplevel_location
{
  my $location = 'AL359765.6.1.13780:2-100';
  
  note "Testing $location by asking for seqlevel";
  my $seqlevel_slice = $slice_adaptor->fetch_by_location($location, 'seqlevel');
  test_slice($location, $seqlevel_slice, 'contig', 'AL359765.6.1.13780', 2, 100, 1);
  
  note "Testing $location by asking for contig";
  my $contig_slice = $slice_adaptor->fetch_by_location($location, 'contig');
  test_slice($location, $contig_slice, 'contig', 'AL359765.6.1.13780', 2, 100, 1);
}

{
  #Non standard name check
  my ($name, $start, $end, $strand) = $slice_adaptor->parse_location_to_values('GL21446.1');
  is($name, 'GL21446.1', 'Name parses out');
  ok(!defined $start, 'Start is undefined');
  ok(!defined $end, 'End is undefined');
  ok(!defined $strand, 'Strand is undefined');
}


## Test patch data

my $patch_db    = $multi->get_DBAdaptor('patch');
$slice_adaptor = $patch_db->get_SliceAdaptor();

my $unique_slices = $slice_adaptor->fetch_by_region_unique("chromosome", "HG1304_PATCH");
my $unique_slice = $unique_slices->[0];
my $initial_slice = $slice_adaptor->fetch_by_region("chromosome", "HG1304_PATCH");
my $exceptions = $initial_slice->get_all_AssemblyExceptionFeatures();
is($unique_slice->start, $exceptions->[0]->start, "Start of patch");
is($unique_slice->end, $exceptions->[0]->end, "End of patch");


## Test karyotype data

my $band_slice = $slice_adaptor->fetch_by_chr_band('6', 'p21.33');
is($band_slice->seq_region_name, '6', 'Fetched chromosome 6');


############# METHODS BELOW HERE 

sub test_toplevel_location {
  my ($location, $cs_name, $seq_region_name, $start, $end, $strand, $ucsc) = @_;
  my $no_warnings = 1;
  my $no_fuzz = undef;
  my $incoming_slice = $slice_adaptor->fetch_by_toplevel_location($location, $no_warnings, $no_fuzz, $ucsc);
  test_slice($location, $incoming_slice, $cs_name, $seq_region_name, $start, $end, $strand);
  return;
}

sub test_slice {
  my ($location, $incoming_slice, $cs_name, $seq_region_name, $start, $end, $strand) = @_;
  $strand ||= 1;
  my $def = ok(defined $incoming_slice, "We expect a defined Slice for $location");
  SKIP : {
    skip 'Incoming slice is undefined', 5 if ! $def;
    is($incoming_slice->coord_system_name(), $cs_name, "Checking coord system name for $location");
    is($incoming_slice->seq_region_name(), $seq_region_name, "Checking seq region name for $location");
    is($incoming_slice->start(), $start, "Checking start for $location");
    is($incoming_slice->end(), $end, "Checking end for $location");
    is($incoming_slice->strand(), $strand, "Checking strand for $location");
  }
  return;
}

sub print_slices {
  my $slices = shift;
  foreach my $slice (@$slices) {
    debug(($slice) ? $slice->name() : "UNDEF");
  } 
  debug( "Got ". scalar(@$slices));
}

sub print_features {
  my $fs = shift;
  foreach my $f (@$fs) {
    my $start  = $f->start();
    my $end    = $f->end();
    my $strand = $f->strand();
    debug("  $start-$end($strand)");
  }
}

done_testing();
