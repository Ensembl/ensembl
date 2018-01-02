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
use Test::Warnings qw( allow_warnings );
use Test::Exception;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;

our $verbose = 0; #set to 1 to turn on debug printouts

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

ok(1);



my $db = $multi->get_DBAdaptor( 'core' );

ok($db);


# Exon specific tests

my $exonad = $db->get_ExonAdaptor();
my $slice_adaptor = $db->get_SliceAdaptor();

my $slice = $slice_adaptor->fetch_by_region('chromosome', '20',
                                            30_811_000,
                                            32_000_000);
ok($exonad);

my $exon = Bio::EnsEMBL::Exon->new();


$exon->start(1000);
ok(&test_getter_setter($exon, 'start', 200));

$exon->end(1400);
ok(&test_getter_setter($exon, 'end', 400));

$exon->strand(1);
ok(&test_getter_setter($exon, 'strand', -1));

$exon->phase(0);
ok(&test_getter_setter($exon, 'phase', -1));

$exon->slice( $slice );
ok(&test_getter_setter($exon, 'slice', $slice));

# should try to store (!)
$exon->end_phase( -1 );
ok(&test_getter_setter($exon, 'end_phase', 1));

ok( test_getter_setter( $exon, "created_date", time() ));
ok( test_getter_setter( $exon, "modified_date", time() ));

#
# find supporting evidence for the exon
#
my @evidence = ();
my @fs = ();
push @fs, @{$db->get_DnaAlignFeatureAdaptor->fetch_all_by_Slice($slice)};
push @fs, @{$db->get_ProteinAlignFeatureAdaptor->fetch_all_by_Slice($slice)};

while(my $f = shift @fs) {
  #debug("feature at: " . $f->start . "-" . $f->end);
  next if $f->start > $exon->end || $f->end < $exon->start;
  push(@evidence, $f);
  # cheat it into storing it again
  $f->dbID( undef );
  $f->adaptor( undef );
}

my $count = scalar(@evidence);
debug("adding $count supporting features");
$exon->add_supporting_features(@evidence);

$multi->hide( "core", "exon", "supporting_feature", 
	      "protein_align_feature", "dna_align_feature");

# We get some 'datatype mismatch: bind param (11) 3.2e-42 as float' warnings
#
allow_warnings(1) if $db->dbc->driver() eq 'SQLite';

$exonad->store($exon);

allow_warnings(0) if $db->dbc->driver() eq 'SQLite';

ok($exon->dbID() && $exon->adaptor == $exonad);

# now test fetch_by_dbID

my $newexon = $exonad->fetch_by_dbID($exon->dbID);

ok($newexon);



debug("exon chr start  = " . $exon->start);
debug("exon chr end    = " . $exon->end);
debug("exon chr strand = " . $exon->strand);

debug("newexon start  = " . $newexon->start());
debug("newexon end    = " . $newexon->end());
debug("newexon strand = " . $newexon->strand());

ok($newexon->start == 30811999 &&
   $newexon->end == 30812399 &&
   $newexon->strand==1);


#
# Test transform to another slice
#
$slice = $exon->slice();
$slice = $db->get_SliceAdaptor->fetch_by_region('chromosome',
                                         $slice->seq_region_name,
                                         $slice->start + $exon->start - 11,
                                         $slice->start + $exon->end + 9);
$exon = $exon->transfer($slice);
debug("exon start  = " . $exon->start);
debug("exon end    = " . $exon->end);
debug("exon strand = " . $exon->strand);
ok($exon->start == 11 && $exon->end == 411 && $exon->strand==1);


#
# Test Transform to contig coord system
#
$exon = $exon->transform('contig');

debug("exon start  = " . $exon->start);
debug("exon end    = " . $exon->end);
debug("exon strand = " . $exon->strand);
debug("exon seq_region = " . $exon->slice->seq_region_name);

ok($exon->start == 913);
ok($exon->end   == 1313);
ok($exon->strand == 1);
ok($exon->slice->seq_region_name eq 'AL034550.31.1.118873');


#regression test, supporting evidence was lost post transform before...
my $se_count = scalar(@{$exon->get_all_supporting_features});

debug("Got $se_count supporting feature after transform");
ok($se_count == $count);
$exon->flush_supporting_features;
is(scalar(@{$exon->get_all_supporting_features}), 0, "All supporting features flushed");

#make sure that supporting evidencewas stored correctly
$se_count = scalar(@{$newexon->get_all_supporting_features});
debug("Got $se_count from newly stored exon");
ok($se_count == $count);


# list_ functions
debug ("Exon->list_dbIDs");
my $ids = $exonad->list_dbIDs();
ok (@{$ids});

debug ("Exon->list_stable_ids");
my $stable_ids = $exonad->list_stable_ids();
ok (@{$stable_ids});

#hashkey
my $hashkey = $exon->hashkey();
debug($hashkey);

ok($hashkey eq $exon->slice->name . '-' . $exon->start . '-' .
               $exon->end . '-' . $exon->strand . '-' . $exon->phase .
               '-' . $exon->end_phase);

$multi->restore();


# regression test
# make sure that sequence fetching and caching is not broken
$exon->stable_id('TestID');
my $first_seq = $exon->seq();
my $second_seq = $exon->seq();

ok($first_seq->seq() && $first_seq->seq() eq $second_seq->seq());
ok($first_seq->display_id()  && $first_seq->display_id() eq $second_seq->display_id());


#
# test the remove method works
#

$multi->save("core", "exon", "supporting_feature",
  "dna_align_feature", "protein_align_feature", 'meta_coord');

my $ex_count = count_rows($db, 'exon');
my $supfeat_count = count_rows($db, 'supporting_feature');

$exon = $exonad->fetch_by_stable_id('ENSE00000859937');

# check the created and modified times
my @date_time = localtime( $exon->created_date());
ok( $date_time[3] == 6 && $date_time[4] == 11 && $date_time[5] == 104 );

@date_time = localtime( $exon->modified_date());
ok( $date_time[3] == 6 && $date_time[4] == 11 && $date_time[5] == 104 );


my $supfeat_minus = @{$exon->get_all_supporting_features()};

$exonad->remove($exon);

ok(count_rows($db, 'exon') == $ex_count - 1);
ok(count_rows($db, 'supporting_feature') == $supfeat_count - $supfeat_minus);

$multi->restore();

#
# tests for multiple versions of transcripts in a database
#

$exon = $exonad->fetch_by_stable_id('ENSE00001109603');
debug("fetch_by_stable_id");
ok( $exon->dbID == 162033 );

$exon->stable_id_version('ENSE00000171455.4');
is($exon->stable_id, 'ENSE00000171455', 'Stable id set with stable_id_version');
is($exon->version, 4, 'Version set with stable_id_version');
is($exon->stable_id_version, 'ENSE00000171455.4', 'Stable id and version from stable_id_version');

$exon->stable_id_version('ENSE00000171456');
is($exon->stable_id, 'ENSE00000171456', 'Stable id set with stable_id_version');
is($exon->version, undef, 'Version undef from stable_id_version');
is($exon->stable_id_version, 'ENSE00000171456', 'Stable id and no version from stable_id_version');

$exon = $exonad->fetch_by_stable_id('ENSE00001109603.1');
ok($exon->dbID == 162033, 'fetch_by_stable_id with version');

$exon = $exonad->fetch_by_stable_id('ENSE00001109603.1a');
ok(!defined($exon), 'fetch_by_stable_id with bad version');

$exon = $exonad->fetch_by_stable_id_version('ENSE00001109603', 1);
ok($exon->dbID == 162033, 'fetch_by_stable_id_version with version');

$exon = $exonad->fetch_by_stable_id_version('ENSE00001109603', '1a');
ok(!defined($exon), 'fetch_by_stable_id_version with bad version');

my @exons = @{ $exonad->fetch_all_versions_by_stable_id('ENSE00001109603') };
debug("fetch_all_versions_by_stable_id");
ok( scalar(@exons) == 1 );

# store/update tests

$multi->hide( "core", "exon", "supporting_feature", 
	      "protein_align_feature", "dna_align_feature", 'meta_coord');

my $e1 = Bio::EnsEMBL::Exon->new(
  -start => 10,
  -end => 1000,
  -strand => 1,
  -slice => $slice,
  -phase => 0,
  -end_phase => 0,
  -stable_id => 'ENSE0001'
);

my $e2 = Bio::EnsEMBL::Exon->new(
  -start => 10,
  -end => 1000,
  -strand => 1,
  -slice => $slice,
  -phase => 0,
  -end_phase => 0,
  -stable_id => 'ENSE0001',
  -version => 2,
  -is_current => 0
);

my $e3 = Bio::EnsEMBL::Exon->new(
  -start => 10,
  -end => 1000,
  -strand => 1,
  -slice => $slice,
  -phase => 0,
  -end_phase => 0,
  -stable_id => 'ENSE0001',
  -is_current => 0
);
$e3->version(undef);

$exonad->store($e1);
$exonad->store($e2);
$exonad->store($e3);

$exon = $exonad->fetch_by_stable_id('ENSE0001');
ok( $exon->is_current == 1);
ok( $exon->version == 1);

@exons = @{ $exonad->fetch_all_versions_by_stable_id('ENSE0001') };
foreach my $e (@exons) {
  if (defined $e->version && $e->version == 2) {
    ok($e->is_current == 0);
  }
}

my $null_versions = 0;
foreach my $e (@exons) {
  if (! defined $e->version) {
    $null_versions++;
  }
}
is ( $null_versions, 1, "Null/undef version stored and retrieved");

$multi->restore();

# TESTS 36-47: Tests for cdna_start(), cdna_end(), cdna_coding_start(),
# cdna_coding_end(), coding_region_start(), and coding_region_end().

my $transcriptad = $db->get_TranscriptAdaptor();
my $transcript   = $transcriptad->fetch_by_stable_id('ENST00000246229');

@exons = @{ $transcript->get_all_Exons() };

$exon = shift @exons;    # First exon is non-coding.

ok( $exon->cdna_start($transcript) == 1 );
ok( $exon->cdna_end($transcript) == 88 );
ok( !defined $exon->cdna_coding_start($transcript) );
ok( !defined $exon->cdna_coding_end($transcript) );
ok( !defined $exon->coding_region_start($transcript) );
ok( !defined $exon->coding_region_end($transcript) );

is ( $exon->rank($transcript), 1, "First exon has rank 1");

$exon = shift @exons;    # Second exon is coding.

ok( $exon->cdna_start($transcript) == 89 );
ok( $exon->cdna_end($transcript) == 462 );
ok( $exon->cdna_coding_start($transcript) == 203 );
ok( $exon->cdna_coding_end($transcript) == 462 );
ok( $exon->coding_region_start($transcript) == 30577779 );
ok( $exon->coding_region_end($transcript) == 30578038 );

is ( $exon->rank($transcript), 2, "Second exon has rank 2");

my $pep = $exon->peptide($transcript);
is($pep->seq, 'MTTFFTSVPPWIQDAKQEEEVGWKLVPRPRGREAESQVKCQCEISGTPFSNGEKLRPHSLPQPEQRPYSCPQLHCGKAFASKYKLYR', 'Retrieved peptide sequence');

SKIP: {
  skip 'No registry support for SQLite yet', 1 if $db->dbc->driver() eq 'SQLite';

  #test the get_species_and_object_type method from the Registry
  my $registry = 'Bio::EnsEMBL::Registry';
  my ( $species, $object_type, $db_type ) = $registry->get_species_and_object_type('ENSE00000859937');
  ok( $species eq 'homo_sapiens' && $object_type eq 'Exon');
}

# UTR and coding region tests. Only testing simple +ve orientation transcript ATMO but it is a start
# tests are based on offsetted coordinates from ENST00000000233 in release 67
{
  my $base_cs = Bio::EnsEMBL::CoordSystem->new(-NAME => 'chromosome', -RANK => 1);
  my $base_slice = Bio::EnsEMBL::Slice->new(-COORD_SYSTEM => $base_cs, -SEQ_REGION_NAME => 'a', -STRAND => 1, -START => 1, -END => 2000, -SEQ => 'A'x2000, -SEQ_REGION_LENGTH => 2000);
  my $base_transcript = Bio::EnsEMBL::Transcript->new(
    -START => 99,
    -END => 1759,
    -STRAND => 1,
    -SLICE => $base_slice
  );
  
  my $start_exon = Bio::EnsEMBL::Exon->new(-START => 99, -END => 319, -STRAND => 1, -STABLE_ID => 'Exon1', -SLICE => $base_slice);
  throws_ok { $start_exon->rank($base_transcript) } qr/does not have/, "No exons in transcript";
  $base_transcript->add_Exon($start_exon);
  is ($start_exon->rank($base_transcript), 1, "Start exon in position 1");
  my $end_exon = Bio::EnsEMBL::Exon->new(-START => 1267, -END => 1759, -STRAND => 1, -STABLE_ID => 'Exon2', -SLICE => $base_slice);
  throws_ok { $end_exon->rank($base_transcript) } qr/does not belong/, "Exon does not belong to transcript";
  $base_transcript->add_Exon($end_exon);

  $base_transcript->translation(Bio::EnsEMBL::Translation->new(
    -START_EXON => $start_exon,
    -SEQ_START => 155,
    -END_EXON => $end_exon,
    -SEQ_END => 87
  ));

  is($start_exon->cdna_coding_start($base_transcript), 155, 'Coding starts at 155bp into the first exon');
  is($start_exon->cdna_coding_end($base_transcript), 221, 'Coding ends at 221bp in the first exon (at the exon end)');
  is($start_exon->coding_region_start($base_transcript), (99+155)-1, 'Coding starts at an offset of 99bp plus coding start');
  is($start_exon->coding_region_end($base_transcript), (99+221)-1, 'Coding region end is the offset of 99bp with the exon length');
  
  is($end_exon->cdna_coding_start($base_transcript), 222, 'CDNA coding start in last exon should be first exon + 1bp');
  is($end_exon->cdna_coding_end($base_transcript), (222+87)-1, 'CDNA coding end should be 86 plus first exon length');
  is($end_exon->coding_region_start($base_transcript), 1267, 'Seq region location start is same as exon start');
  is($end_exon->coding_region_end($base_transcript), (1267+87)-1, 'Seq region location end is offsetted by exon coding length');
}

# Checking the reverse strand. Taken from ENST00000321407 in E! 66 with 66
# API to check for regressions
{
  my $base_cs = Bio::EnsEMBL::CoordSystem->new(-NAME => 'chromosome', -RANK => 1);
  my $base_slice = Bio::EnsEMBL::Slice->new(
    -COORD_SYSTEM => $base_cs, -SEQ_REGION_NAME => 'a', 
    -STRAND => 1, -START => 1, -END => 6000, -SEQ => 'A'x6000, 
    -SEQ_REGION_LENGTH => 6000
  );
  my $base_transcript = Bio::EnsEMBL::Transcript->new(
    -START => 672,
    -END => 5661,
    -SLICE => $base_slice,
    -STRAND => -1
  );
  
  
  my $start_exon = Bio::EnsEMBL::Exon->new(-START => 4205, -END => 5661, -STRAND => -1, -SLICE => $base_slice);
  $base_transcript->add_Exon($start_exon);
  my $end_exon = Bio::EnsEMBL::Exon->new(-START => 672, -END => 3363, -STRAND => -1, -SLICE => $base_slice);
  $base_transcript->add_Exon($end_exon);
  
  $base_transcript->translation(Bio::EnsEMBL::Translation->new(
    -START_EXON => $start_exon,
    -SEQ_START => 426,
    -END_EXON => $end_exon,
    -SEQ_END => 1296
  ));

  
  is($start_exon->cdna_coding_start($base_transcript), 426, 'CDNA start equals SEQ_START');
  is($start_exon->cdna_coding_end($base_transcript), 1457, 'CDNA end equals SEQ_START plus exon length');
  is($start_exon->coding_region_start($base_transcript), 4205, 'Coding region start is start of first exon');
  is($start_exon->coding_region_end($base_transcript), 5236, 'Coding region end is start of first exon plus its length');
  is($start_exon->frame, 2, 'Coding start has frame 2');
  
  is($end_exon->cdna_coding_start($base_transcript), 1458, 'CDNA coding start equals END of first Exon + 1');
  is($end_exon->cdna_coding_end($base_transcript), 2753, 'CDNA coding end equals start plus its length into the last exon');
  is($end_exon->coding_region_start($base_transcript), 2068, 'Start is the end of the last exon minus the coding length');
  is($end_exon->coding_region_end($base_transcript), 3363, 'End is the same as the end exon end');

}

# Testing exon is_coding
{
  # Create a transcript with start and end same as exon and the coding regions falls within their boundries - POSITIVE STRAND
  my $sa = $db->get_SliceAdaptor();
  my $local_slice = $sa->fetch_by_region('chromosome', "20");
  #create a transcript instance
  my $tr = Bio::EnsEMBL::Transcript->new(-SLICE => $local_slice, -START => 2000, -END => 3000, -STRAND => 1); 
  #create an exon instance
  my $start_exon = Bio::EnsEMBL::Exon->new(-START => 2000, -END => 3000, -STRAND => 1, -SLICE => $local_slice);
  $tr->add_Exon($start_exon);
  $tr->translation(Bio::EnsEMBL::Translation->new(
      -SEQ_START => 100,
      -SEQ_END => 500,
      -START_EXON => $start_exon,
      -END_EXON => $start_exon,
    ));


  my $exon_one = $tr->get_all_Exons()->[0];

  ok($tr->translate, "Transcript can translate");
  is($exon_one->start, $tr->start, 'Exon start equals Transcript start');
  is($exon_one->end, $tr->end, 'Exon end equals Transcript end');

  is($exon_one->cdna_coding_start($tr), 100, 'CDNA coding start equals SEQ_START');
  is($exon_one->cdna_coding_end($tr), 500, 'CDNA coding end equals SEQ_END');

  is($exon_one->coding_region_start($tr), $tr->coding_region_start, 'coding_region_start is 2099'); # $exon_one->start + 100 -1 = 2099
  is($exon_one->coding_region_end($tr), $tr->coding_region_end, 'coding_region_end is 2499'); # $tr->coding_region_start + 500 -100 = 2499

  ok($tr->coding_region_start > $exon_one->start,  'coding_region_start > exon_start');
  ok($tr->coding_region_end < $exon_one->end,  'coding_region_end < exon_end');

  my $is_coding = $exon_one->is_coding($tr);
  is($is_coding, 1, "Exon is coding");

  # REPEAT WITH NEGATIVE STRAND
  #create a transcript instance
  $tr = Bio::EnsEMBL::Transcript->new(-SLICE => $local_slice, -START => 2000, -END => 3000, -STRAND => -1); 
  #create an exon instance
  $start_exon = Bio::EnsEMBL::Exon->new(-START => 2000, -END => 3000, -STRAND => -1, -SLICE => $local_slice);
  $tr->add_Exon($start_exon);
  $tr->translation(Bio::EnsEMBL::Translation->new(
      -SEQ_START => 100,
      -SEQ_END => 500,
      -START_EXON => $start_exon,
      -END_EXON => $start_exon,
    ));


  $exon_one = $tr->get_all_Exons()->[0];

  ok($tr->translate, "Transcript can translate");
  is($exon_one->start, $tr->start, 'Exon start equals Transcript start');
  is($exon_one->end, $tr->end, 'Exon end equals Transcript end');

  is($exon_one->cdna_coding_start($tr), 100, 'CDNA coding start equals SEQ_START');
  is($exon_one->cdna_coding_end($tr), 500, 'CDNA coding end equals SEQ_END');

  is($exon_one->coding_region_start($tr), $tr->coding_region_start, 'coding_region_start is 2099');
  is($exon_one->coding_region_end($tr), $tr->coding_region_end, 'coding_region_end is 2499');

  ok($tr->coding_region_start > $exon_one->start,  'coding_region_start > exon_start');
  ok($tr->coding_region_end < $exon_one->end,  'coding_region_end < exon_end');

  $is_coding = $exon_one->is_coding($tr);
  is($is_coding, 1, "Exon is coding");

}

done_testing();
