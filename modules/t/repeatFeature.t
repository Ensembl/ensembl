use strict;

use Bio::EnsEMBL::RepeatFeature;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::RepeatConsensus;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::RepeatFeature;


BEGIN { $| = 1;
	use Test;
	plan tests => 16;
}


use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0;


my $coord_system = Bio::EnsEMBL::CoordSystem->new
  (-NAME    => 'chromosome',
   -VERSION => 'NCBI34',
   -DBID    => 123,
   -RANK    => 1);

my $analysis = Bio::EnsEMBL::Analysis->new(-LOGIC_NAME => 'test');

my $slice = Bio::EnsEMBL::Slice->new(-COORD_SYSTEM    => $coord_system,
                                     -SEQ_REGION_NAME => 'X',
                                     -SEQ_REGION_LENGTH => 15e6,
                                     -START           => 1_000_000,
                                     -END             => 2_000_000);

my $repeat_consensus = Bio::EnsEMBL::RepeatConsensus->new
  (-REPEAT_CONSENSUS => 'ACTG',
   -NAME             => 'ACTG(n)',
   -LENGTH           => 4,
   -REPEAT_CLASS    => 'Simple_repeat');


#
# Test new and getters
#

my $start  = 10;
my $end    = 100;
my $strand = -1;
my $hstart = 1;
my $hend = 90;
my $dbID = 123;
my $score = 12.5;

my $rf = Bio::EnsEMBL::RepeatFeature->new
  (-START   => $start,
   -END     => $end,
   -STRAND  => $strand,
   -ANALYSIS => $analysis,
   -SLICE   => $slice,
   -HSTART  => $hstart,
   -HEND    => $hend,
   -SCORE   => $score,
   -REPEAT_CONSENSUS => $repeat_consensus);


ok($rf && $rf->isa('Bio::EnsEMBL::RepeatFeature'));

ok($rf->start == $start);
ok($rf->end == $end);
ok($rf->strand == $strand);
ok($rf->analysis == $analysis);
ok($rf->slice == $slice);
ok($rf->hstart == $hstart);
ok($rf->hend == $hend);
ok($rf->score == $score);
ok($rf->repeat_consensus == $repeat_consensus);


#
# Test Getter/Setters
#

$repeat_consensus = Bio::EnsEMBL::RepeatConsensus->new
  (-REPEAT_CONSENSUS => 'ACTG',
   -NAME             => 'ACTG(n)',
   -LENGTH           => 4,
   -REPEAT_CLASS    => 'Simple_repeat');


ok(test_getter_setter($rf, 'hstart', 120));
ok(test_getter_setter($rf, 'hend', 200));
ok(test_getter_setter($rf, 'repeat_consensus', $repeat_consensus));
ok(test_getter_setter($rf, 'score', '45.5'));

ok($rf->display_id eq 'ACTG(n)');

ok($rf->hstrand == 1);
