
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
verbose('INFO');

my $multi_db = MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');


my $dfa = $db->get_DensityFeatureAdaptor();

ok(ref($dfa) && $dfa->isa('Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor'));

my $slice_adaptor = $db->get_SliceAdaptor();
my $slice = $slice_adaptor->fetch_by_region('chromosome', '20', 1, 600);

my $dfs = $dfa->fetch_all_by_Slice($slice, 'SNPDensity');
print_features($dfs);

my $dfs = $dfa->fetch_all_by_Slice($slice, 'SNPDensity', 10, 1);
print_features($dfs);


sub print_features {
  my $features = shift;

  debug("\n");
  foreach my $f (@$features) {
    my $start = $f->start();
    my $end   = $f->end();
    my $density_value = $f->density_value();
    my $n = int(0.75 * $density_value);
    debug(('*'x$n)."($start-$end)");
  }
  debug("\n");
}





