use strict;
use warnings;

use lib 't';
use Bio::EnsEMBL::CoordSystem;


BEGIN { $| = 1;
	use Test;
	plan tests => 9;
}


my $name    = 'chromosome';
my $version = 'NCBI33';
my $dbID    = 1;
my $top_level = 1;
my $sequence_level = 0;


#
# 1 Test constructor
#
my $coord_system = Bio::EnsEMBL::CoordSystem->new
  (-NAME    => $name,
   -VERSION => $version,
   -DBID    => $dbID,
   -TOP_LEVEL => $top_level,
   -SEQUENCE_LEVEL => $sequence_level);


ok($coord_system && $coord_system->isa('Bio::EnsEMBL::CoordSystem'));


#
# 2 Test name()
#

ok($coord_system->name() eq $name);


#
# 3 Test version()
#

ok($coord_system->version() eq $version);


#
# 4 Test is_top_level()
#
ok($coord_system->is_top_level());

#
# 5 Test is_sequence_level()
#
ok(!$coord_system->is_sequence_level());

#
# 6-7 Test equals()
#

my $coord_system2 = Bio::EnsEMBL::CoordSystem->new
  (-NAME    => $name,
   -VERSION => $version,
   -DBID    => 123,
   -TOP_LEVEL => $top_level);

ok($coord_system->equals($coord_system2));

$coord_system2 = Bio::EnsEMBL::CoordSystem->new
  (-NAME    => 'chromosome',
   -VERSION => 'NCBI34',
   -DBID    => 123,
   -TOP_LEVEL => $top_level);

ok(!$coord_system->equals($coord_system2));




