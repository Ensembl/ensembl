use strict;
use warnings;

use lib 't';

BEGIN { $| = 1;
	use Test;
	plan tests => 11;
}

use TestUtils qw( debug test_getter_setter );

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::CoordSystem;

our $verbose= 0;

my $coord_system = Bio::EnsEMBL::CoordSystem->new
  (-NAME    => 'chromosome',
   -VERSION => 'NCBI34',
   -DBID    => 123,
   -TOP_LEVEL => 1);

my $analysis = Bio::EnsEMBL::Analysis->new(-LOGIC_NAME => 'test');

my $slice = Bio::EnsEMBL::Slice->new(-COORD_SYSTEM    => $coord_system,
                                     -SEQ_REGION_NAME => 'X',
                                     -START           => 1_000_000,
                                     -END             => 2_000_000);

#
# 1-6 Test new and getters
#
my $start  = 10;
my $end    = 100;
my $strand = -1;
my $feature = Bio::EnsEMBL::Feature->new(-START => 10,
                                         -END   => 100,
                                         -STRAND => -1,
                                         -ANALYSIS => $analysis,
                                         -SLICE => $slice);


ok($feature && $feature->isa('Bio::EnsEMBL::Feature'));

ok($feature->start == $start);
ok($feature->end == $end);
ok($feature->strand == $strand);
ok($feature->analysis == $analysis);
ok($feature->slice == $slice);

#
# 7-11 Test setters
#
$analysis = Bio::EnsEMBL::Analysis->new(-LOGIC_NAME => 'new analysis');
$slice = Bio::EnsEMBL::Slice->new(-COORD_SYSTEM    => $coord_system,
                                  -SEQ_REGION_NAME => 'Y',
                                  -START           => 1_000_000,
                                  -END             => 2_000_000);
ok(&test_getter_setter($feature, 'start', 1000));
ok(&test_getter_setter($feature, 'end', 2000));
ok(&test_getter_setter($feature, 'strand', 1));
ok(&test_getter_setter($feature, 'analysis', $analysis));
ok(&test_getter_setter($feature, 'slice', $slice));


#
# 8 Test length
#

ok($feature->length == $start-$end+1);


#
# 9-11 Test move
#

$feature->move(300, 500, 1);

ok($feature->start == 300);
ok($feature->end == 500);
ok($feature->strand == 1);

