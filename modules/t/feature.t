use strict;
use warnings;

use lib 't';

BEGIN { $| = 1;
	use Test;
	plan tests => 75;
}

use TestUtils qw( debug test_getter_setter );

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::CoordSystem;
use MultiTestDB;

our $verbose= 0; #turn on or off debug statements


my $multi_db = MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');

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
# Test new and getters
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
# Test setters
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
# Test length
#

ok($feature->length == ($end-$start+1));


#
# Test move
#

$feature->move(300, 500, 1);

ok($feature->start == 300);
ok($feature->end == 500);
ok($feature->strand == 1);


#################
# Test transform
#################

$slice = $db->get_SliceAdaptor->fetch_by_region('chromosome',
                                                 '20',
                                                 30_249_935,
                                                 31_000_000);

$feature->slice($slice);

#
# Test Transform chromosome -> contig
#

$feature = $feature->transform('contig');

debug("\nchromosome -> contig");
debug("start  = " . $feature->start());
debug("end    = " . $feature->end());
debug("strand = " . $feature->strand());
debug("seq_region = " . $feature->slice->seq_region_name());

ok($feature->start() == 400);
ok($feature->end() == 600);
ok($feature->strand() == 1);
ok($feature->slice->seq_region_name() eq 'AL359765.6.1.13780');
ok($feature->slice->coord_system->name() eq 'contig');


#
# Test Transform contig -> supercontig
#

$feature = $feature->transform('supercontig');

debug("\ncontig -> supercontig");
debug("start  = " . $feature->start());
debug("end    = " . $feature->end());
debug("strand = " . $feature->strand());
debug("seq_region = " . $feature->slice->seq_region_name());

ok($feature->start()  == 658269);
ok($feature->end()    == 658469);
ok($feature->strand() == 1);
ok($feature->slice->seq_region_name() eq 'NT_028392');
ok($feature->slice->coord_system->name() eq 'supercontig');


#
# Test Transform supercontig -> contig
#

$feature = $feature->transform('contig');
debug("\nsupercontig -> contig");
debug("start  = " . $feature->start());
debug("end    = " . $feature->end());
debug("strand = " . $feature->strand());
debug("seq_region = " . $feature->slice->seq_region_name());

ok($feature->start() == 400);
ok($feature->end() == 600);
ok($feature->strand() == 1);
ok($feature->slice->seq_region_name() eq 'AL359765.6.1.13780');
ok($feature->slice->coord_system->name() eq 'contig');


#
# Test Transform contig -> clone
#

$feature = $feature->transform('clone');
debug("\ncontig -> clone");
debug("start  = " . $feature->start());
debug("end    = " . $feature->end());
debug("strand = " . $feature->strand());
debug("seq_region = " . $feature->slice->seq_region_name());

ok($feature->start() == 400);
ok($feature->end() == 600);
ok($feature->strand() == 1);
ok($feature->slice->seq_region_name() eq 'AL359765.6');
ok($feature->slice->coord_system->name() eq 'clone');


#
# Test transform to into gap
#
$slice = $db->get_SliceAdaptor->fetch_by_region('chromosome',
                                                 '20',
                                                 1_000_000,
                                                 2_000_000);
$feature->slice($slice);

ok(!defined($feature->transform('contig')));



###############
# Test transfer
###############
$slice = $db->get_SliceAdaptor->fetch_by_region('chromosome',
                                                 '20',
                                                 30_249_935,
                                                 31_000_000);

$feature->slice($slice);


#
# Transfer to expanded inverted chr slice
#

#get larger slice on other strand
$slice = $slice->invert()->expand(100,100);

$feature = $feature->transfer($slice);
debug("\ntransfer to inverted chromosomal slice");
debug("start  = " . $feature->start());
debug("end    = " . $feature->end());
debug("strand = " . $feature->strand());
debug("seq_region = " . $feature->slice->seq_region_name());

ok($feature->start()  == 749567);
ok($feature->end()    == 749767);
ok($feature->strand() == -1);


#
# Transfer to contig slice
#

$slice = $db->get_SliceAdaptor->fetch_by_region('contig',
                                                'AL359765.6.1.13780',
                                                '30', '3000');

$feature = $feature->transfer($slice);
debug("\ntransfer to contig slice");
debug("start  = " . $feature->start());
debug("end    = " . $feature->end());
debug("strand = " . $feature->strand());
debug("seq_region = " . $feature->slice->seq_region_name());


ok($feature->start()  == 471);
ok($feature->end()    == 671);
ok($feature->strand() == 1);


#
# Transfer to supercontig slice
#
$slice  = $db->get_SliceAdaptor->fetch_by_region('supercontig',
                                                 'NT_028392');
debug( "----------------------------" );

$feature = $feature->transfer($slice);
debug("\ntransfer to supercontig slice");
debug("start  = " . $feature->start());
debug("end    = " . $feature->end());
debug("strand = " . $feature->strand());
debug("seq_region = " . $feature->slice->seq_region_name());

ok($feature->start()  == 658369);
ok($feature->end()    == 658569);
ok($feature->strand() == 1);


#
# Project over a contig boundary
#

$feature->move(671600,671800);

debug( "Before projection to contig" );
debug("start  = " . $feature->start());
debug("end    = " . $feature->end());
debug("strand = " . $feature->strand());
debug("seq_region = " . $feature->slice->seq_region_name());

my @projection = @{$feature->project('contig')};

debug( "After project to contig" );
foreach my $segment (@projection) {
  ($start, $end, $slice) = @$segment;
  debug("[$start-$end] -> " . $slice->seq_region_name
        . ' ' . $slice->start . '-' . $slice->end . '('.$slice->strand().')');
}
debug( "-----------------------" );

ok(@projection == 2);

my ($seg1, $seg2) = @projection;

ok($seg1 && $seg1->[0] == 1);
ok($seg1 && $seg1->[1]   == 50);
ok($seg1 && $seg1->[2]->start == 13731);
ok($seg1 && $seg1->[2]->end   == 13780);
ok($seg1 && $seg1->[2]->seq_region_name eq 'AL359765.6.1.13780');

ok($seg2 && $seg2->[0] == 51);
ok($seg2 && $seg2->[1] == 201);
ok($seg2 && $seg2->[2]->start == 101);
ok($seg2 && $seg2->[2]->end   == 251);
ok($seg2 && $seg2->[2]->seq_region_name eq 'AL031658.11.1.162976');

debug('forward strand slice feature projection');
foreach my $segment (@projection) {
  ($start, $end, $slice) = @$segment;
  debug("[$start-$end] -> " . $slice->seq_region_name
        . ' ' . $slice->start . '-' . $slice->end . '('.$slice->strand().')');
}

#
# Try a projection using an inverted slice
#

$slice = $feature->slice()->invert();
$feature->slice($slice);
$feature->move(4_391_950, 4_392_150);

@projection = @{$feature->project('contig')};

ok(@projection == 2);

($seg1, $seg2) = @projection;

ok($seg1 && $seg1->[0] == 1 );
ok($seg1 && $seg1->[1]   ==  8 );
ok($seg1 && $seg1->[2]->start == 101);
ok($seg1 && $seg1->[2]->end   == 108);
ok($seg1 && $seg1->[2]->seq_region_name eq 'AL031658.11.1.162976');

ok($seg2 && $seg2->[0] == 9);
ok($seg2 && $seg2->[1] == 201);
ok($seg2 && $seg2->[2]->start == 13588);
ok($seg2 && $seg2->[2]->end   == 13780);
ok($seg2 && $seg2->[2]->seq_region_name eq 'AL359765.6.1.13780');

debug('negative strand slice feature projection');
foreach my $segment (@projection) {
  ($start, $end, $slice) = @$segment;
  debug("[$start-$end] -> " . $slice->seq_region_name
        . ' ' . $slice->start . '-' . $slice->end . '('.$slice->strand().')');
}

#
# Try project to same coord system
#


@projection = @{$feature->project('supercontig')};

ok(@projection == 1);

($seg1) = @projection;
$slice = $feature->slice();
ok($seg1 && $seg1->[0] == 1 );
ok($seg1 && $seg1->[1]   == $feature->length() );
ok($seg1 && $seg1->[2]->start == $slice->end - $feature->end() + 1);
ok($seg1 && $seg1->[2]->end   == $slice->end - $feature->start() + 1);
ok($seg1 && $seg1->[2]->strand == $slice->strand * $feature->strand());
ok($seg1 && $seg1->[2]->seq_region_name eq 'NT_028392');


debug('projection to same coord system');
foreach my $segment (@projection) {
  ($start, $end, $slice) = @$segment;
  debug("[$start-$end] -> " . $slice->seq_region_name
        . ' ' . $slice->start . '-' . $slice->end . '('.$slice->strand().')');
}

ok($feature->display_id eq '');