use lib 't';
use strict;

BEGIN { $| = 1;
	use Test;
	plan tests => 2;
}


use MultiTestDB;
use TestUtils qw(test_getter_setter debug);

our $verbose = 0;

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

my $feats;

my $chr_slice = $db->get_SliceAdaptor->fetch_by_region('chromosome','20',
                                                   30_500_000, 30_700_000);

$feats = $dafa->fetch_all_by_Slice($chr_slice);
debug('---fetching by chromosomal slice---');
debug("Got " . scalar(@$feats) . " features back");
ok(@$feats == 1354);
print_features($feats);


my $ctg_slice;
my $ctg_slice  = $db->get_SliceAdaptor->fetch_by_region('contig',
                                                       'AL031658.11.1.162976',
                                                       1,
                                                       50_000);
$feats = $dafa->fetch_all_by_Slice($ctg_slice);
debug('--- contig AL031658.11.1.162976 (1-50000) features ---');
debug("Got " . scalar(@$feats));
ok(@$feats == 709);
print_features($feats);


#
# Test fetch_by_dbID
#
my $feat = $dafa->fetch_by_dbID(22171863, 'contig');
debug('--- fetching by dbID in contig coords ---');
ok($feat);
print_features([$feat]);


$feat = $dafa->fetch_by_dbID(22171863, 'supercontig');
debug('--- fetching by dbID in supercontig coords ---');
ok($feat);
print_features([$feat]);


#
# Test fetch_by_Slice_and_pid
#
$feats = $dafa->fetch_all_by_Slice_and_pid($chr_slice, '90');
debug('--- fetching by chr Slice and pid (90) ---');
debug("Got " . scalar(@$feats));
ok(@$feats == 253);
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
            " ($analysis) " . $f->percent_id);
    } else {
      debug('UNDEF');
    }
  }
}
