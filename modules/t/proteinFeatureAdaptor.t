use strict;

use Bio::EnsEMBL::Test::TestUtils;

BEGIN { $| = 1;
	use Test;
	plan tests => 2;
}

use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;

# get a core DBAdaptor
my $dba = $multi->get_DBAdaptor("core");


#
# Test get_ProteinFeatureAdaptor works
#
my $pfa = $dba->get_ProteinFeatureAdaptor();

ok($pfa && ref($pfa) && $pfa->isa('Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor'));


my $pfs = $pfa->fetch_all_by_translation_id(21724);

print_features($pfs);

ok(@$pfs == 15);

sub print_features {
  my $features = shift;
  foreach my $f (@$features) {
    if(defined($f)) {
      debug($f->start.'-'.$f->end. ' -> '.$f->hseqname. ':'.
            $f->hstart. '-'.$f->hend); 
    } 
  }
}
