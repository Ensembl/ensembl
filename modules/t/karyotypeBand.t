use strict;
use warnings;

use lib 't';

BEGIN { $| = 1;
	use Test;
	plan tests => 9;
}


use TestUtils qw(debug test_getter_setter);
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
   -TOP_LEVEL => 1);


my $slice = Bio::EnsEMBL::Slice->new(-COORD_SYSTEM    => $coord_system,
                                     -SEQ_REGION_NAME => 'X',
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

#
# test getter/setters
#
ok(test_getter_setter($kb, 'name', 'p.31'));
ok(test_getter_setter($kb, 'start', 12_200_000));
ok(test_getter_setter($kb, 'end',   13_000_000));
ok(test_getter_setter($kb, 'stain', 'gpos60'));

