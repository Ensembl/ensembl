use lib 't';

BEGIN { $| = 1;
	use Test;
	plan tests => 2;
}


use MultiTestDB;
use TestUtils qw(test_getter_setter debug);

our $verbose = 1;

my $multi = MultiTestDB->new();
ok(1);

my $db = $multi->get_DBAdaptor('core');
my $dafa = $db->get_DnaAlignFeatureAdaptor();

# list_dbIDs
debug("DnaAlignFeatureAdaptor->list_dbIDs");
my $ids = $dafa->list_dbIDs();
ok (@{$ids});


#
# Test retrieval via Slice
#

my $slice = $db->get_SliceAdaptor->fetch_by_region('chromosome','20',
                                                   30_500_000, 31_000_000);

my $feats = $dafa->fetch_all_by_Slice($slice);
debug('---fetching by chromosomal slice---');
debug("Got " . scalar(@$feats) . " features back");
ok(@$feats == 11069);
print_features($feats);



sub print_features {
  return if(!$verbose);
  my $features = shift;
  foreach my $f (@$features) {
    if(defined($f)) {
      my $seqname = $f->slice->seq_region_name();
      my $analysis = $f->analysis->logic_name();
      debug($seqname . ' ' . $f->start().'-'.$f->end().'('.$f->strand().
            ') ['. $f->dbID.'] '. $f->cigar_string . ' ' .
            $f->hstart .'-'.$f->hend.' ('.$f->hstrand.')'.$f->score() .
            " ($analysis)");
    } else {
      debug('UNDEF');
    }
  }
}
