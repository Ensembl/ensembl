
use lib 't';

use strict;
use warnings;

use Bio::EnsEMBL::DensityFeature;
use Bio::EnsEMBL::Analysis;

use TestUtils qw(debug test_getter_setter);


our $verbose = 0; #set to 1 to turn on debug printouts



BEGIN { $| = 1;
	use Test;
	plan tests => 6;
}

use TestUtils qw( debug );

my $analysis = Bio::EnsEMBL::Analysis->new(-DBID => 1,
                                           -LOGIC_NAME => 'test');


my $start = 10;
my $end   = 102;
my $density_value = 123;

my $df = Bio::EnsEMBL::DensityFeature->new
  (-start    => $start,
   -end      => $end,
   -analysis => $analysis,
   -density_value => $density_value);

ok($df->start == $start && $df->analysis == $analysis && $df->end == $end);
ok($df->strand == 0);

ok($df->density_value == $density_value);

$df = Bio::EnsEMBL::DensityFeature->new_fast
  ({'start'    => $start,
    'end'      => $end,
    'analysis' => $analysis,
    'density_value' => $density_value});

ok($df->start == $start && $df->analysis == $analysis && $df->end == $end);
ok($df->strand == 0);
ok($df->density_value == $density_value);



