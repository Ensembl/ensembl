
use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(verbose);
use TestUtils qw( debug );
use MultiTestDB;

use lib 't';

BEGIN { $| = 1;
	use Test;
	plan tests => 11;
}



our $verbose = 0;
verbose('WARNING');

my $multi_db = MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');


my $dfa = $db->get_DensityFeatureAdaptor();

ok(ref($dfa) && $dfa->isa('Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor'));

my $slice_adaptor = $db->get_SliceAdaptor();
my $slice = $slice_adaptor->fetch_by_region('chromosome', '20', 1, 600);

my $dfs = $dfa->fetch_all_by_Slice($slice, 'SNPDensity');
ok(@$dfs);
print_features($dfs);
$dfs = $dfa->fetch_all_by_Slice($slice, 'RepeatCoverage');
ok(@$dfs);
print_features($dfs);


$dfs = $dfa->fetch_all_by_Slice($slice, 'SNPDensity', 10, 1);
ok(@$dfs);
print_features($dfs);
$dfs = $dfa->fetch_all_by_Slice($slice, 'RepeatCoverage', 10, 1);
ok(@$dfs);
print_features($dfs);

$dfs = $dfa->fetch_all_by_Slice($slice, 'SNPDensity', 50, 1);
ok(@$dfs);
print_features($dfs);
$dfs = $dfa->fetch_all_by_Slice($slice, 'RepeatCoverage', 50, 1);
ok(@$dfs);
print_features($dfs);

$dfs = $dfa->fetch_all_by_Slice($slice, 'SNPDensity', 2, 1);
ok(@$dfs);
print_features($dfs);
$dfs = $dfa->fetch_all_by_Slice($slice, 'RepeatCoverage', 2, 1);
ok(@$dfs);
print_features($dfs);



$slice_adaptor->fetch_by_region('chromosome', '20', 50, 600);

$dfs = $dfa->fetch_all_by_Slice($slice, 'SNPDensity', 50, 1);
ok(@$dfs);
print_features($dfs);
$dfs = $dfa->fetch_all_by_Slice($slice, 'RepeatCoverage', 50, 1);
ok(@$dfs);
print_features($dfs);



#
# helper to draw an ascii representation of the density features
#
sub print_features {
  my $features = shift;

  return if(!@$features);

  my $sum = 0;
  my $length = 0;
  my $type = $features->[0]->density_value_type();

  debug("\n");
  foreach my $f (@$features) {
    my $draw_width = 1;
    my $density_value = $f->density_value();
    my $draw_height = int(0.75 * $density_value);
    $sum += $density_value;
    $length += $f->length();
    for(my $i = 0; $i < $draw_width; $i++) {
      debug(('*'x$draw_height) . "($density_value)");
    }
  }
  my $avg = undef;
  $avg = $sum/$length if($length < 0);
  debug("Type=$type, Sum=$sum, Length=$length, Avg/Base=$sum");
}





