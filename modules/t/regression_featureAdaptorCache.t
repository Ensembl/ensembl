use strict;
use warnings;
use Test::More;
use Bio::EnsEMBL::Test::MultiTestDB;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $multi->get_DBAdaptor( "core" );

my $sa = $db->get_SliceAdaptor();
my $ga = $db->get_GeneAdaptor();

ok(!$ga->db()->no_cache(), 'Checking cache is on');

my $cache_assert = sub {
  my ($expected) = @_;
  is(scalar(keys %{$ga->{_slice_feature_cache}}), $expected, sprintf('Asserting cache has %d element(s)', $expected));
};

my $run = sub {
  my $start = 30_249_935;
  my $end = 31_254_640;
  my $offset = 0; 
  my @regions = (
    [$start, $end + $offset++],
    [$start, $end + $offset++],
    [$start, $end + $offset++],
    [$start, $end + $offset++],
    [$start, $end + $offset++]
  );
  $ga->fetch_all_by_Slice($sa->fetch_by_region( "chromosome", "20", @{$regions[0]} ));
  
  $cache_assert->(1);  
  foreach my $region (@regions) {
    my $slice = $sa->fetch_by_region( "chromosome", "20", @{$region} );
    my $features = $ga->fetch_all_by_Slice($slice);
  }
  $cache_assert->(4);
};

$run->();
$ga->clear_cache();
$run->();

done_testing();