
use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(verbose);
use TestUtils qw( debug );
use MultiTestDB;

use lib 't';

BEGIN { $| = 1;
	use Test;
	plan tests => 43;
}



our $verbose = 1;
verbose('WARNING');

my $multi_db = MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');


my $dfa = $db->get_DensityFeatureAdaptor();

ok(ref($dfa) && $dfa->isa('Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor'));

my $slice_adaptor = $db->get_SliceAdaptor();
my $slice = $slice_adaptor->fetch_by_region('chromosome', '20', 1, 600);

my $dfs = $dfa->fetch_all_by_Slice($slice, 'SNPDensity');
print_features($dfs);

$dfs = $dfa->fetch_all_by_Slice($slice, 'SNPDensity', 10, 1);
print_features($dfs);

$dfs = $dfa->fetch_all_by_Slice($slice, 'SNPDensity', 50, 1);
print_features($dfs);

$dfs = $dfa->fetch_all_by_Slice($slice, 'SNPDensity', 2, 1);
print_features($dfs);



#
# helper to draw an ascii representation of the density features
#
sub print_features {
  my $features = shift;

  my $sum = 0;
  my $length = 0;

  debug("\n");
  foreach my $f (@$features) {
    my $draw_width = 1;
    my $density_value = $f->density_value();
    my $draw_height = int(0.75 * $density_value);
    $sum += $density_value;
    $length += $f->length();
    for(my $i = 0; $i < $draw_width; $i++) {
      debug(('*'x$draw_height));
    }
  }
  my $avg = undef;
  $avg = $sum/$length if($length < 0);
  debug("\nSum=$sum, Length=$length, Avg/Base=$sum");
}





