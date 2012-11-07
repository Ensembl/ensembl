use strict;
use warnings;
use Test::More;

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::KaryotypeBand;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::CoordSystem;

#
# Test constructor
#

my $coord_system = Bio::EnsEMBL::CoordSystem->new
  (-NAME    => 'chromosome',
   -VERSION => 'NCBI34',
   -DBID    => 123,
   -RANK => 1);


my $slice = Bio::EnsEMBL::Slice->new(-COORD_SYSTEM    => $coord_system,
                                     -SEQ_REGION_NAME => 'X',
                                     -SEQ_REGION_LENGTH => 15e6,
                                     -START           => 1,
                                     -END             => 2e6);

my $start  = 1;
my $end    = 1e6;
my $name   = 'q.11';
my $stain  = 'gpos50';
my $kb = Bio::EnsEMBL::KaryotypeBand->new(-START => $start,
                                          -END   => $end,
                                          -STAIN => $stain,
                                          -NAME  => $name,
                                          -SLICE => $slice);


ok($kb->start() == $start);
ok($kb->end()   == $end);
ok($kb->stain() eq $stain);
ok($kb->name() eq $name);
ok($kb->slice == $slice);
ok($kb->display_id eq $name);

#
# test getter/setters
#
ok(test_getter_setter($kb, 'name', 'p.31'));
ok(test_getter_setter($kb, 'start', 12_200_000));
ok(test_getter_setter($kb, 'end',   13_000_000));
ok(test_getter_setter($kb, 'stain', 'gpos60'));

done_testing();