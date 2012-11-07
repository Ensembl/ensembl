
use Test::More;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::DnaPepAlignFeature;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
ok(1);

my $db = $multi->get_DBAdaptor('core');
my $pafa = $db->get_ProteinAlignFeatureAdaptor();

# list_dbIDs
debug("ProteinAlignFeatureAdaptor->list_dbIDs");
my $ids = $pafa->list_dbIDs();
ok (@{$ids});




#
# Test retrieval via Slice
#

my $feats;

my $chr_slice = $db->get_SliceAdaptor->fetch_by_region('chromosome','20',
                                                   30_500_000, 30_700_000);

$feats = $pafa->fetch_all_by_Slice($chr_slice);
debug('---fetching by chromosomal slice---');
debug("Got " . scalar(@$feats) . " features back");
ok(@$feats == 2429);
print_features($feats);

my $ctg_slice;
$ctg_slice  = $db->get_SliceAdaptor->fetch_by_region('contig',
                                                     'AL031658.11.1.162976',
                                                     1,
                                                     50_000);
$feats = $pafa->fetch_all_by_Slice($ctg_slice);
debug('--- contig AL031658.11.1.162976 (1-50000) features ---');
debug("Got " . scalar(@$feats));
ok(@$feats == 357);
      print_features($feats);


#
# Test fetch_by_dbID
#
my $feat = $pafa->fetch_by_dbID(5339568);
debug('--- fetching by dbID ---');
ok($feat);
print_features([$feat]);
ok($feat->db_name eq 'EMBL');
ok($feat->db_display_name eq 'EMBL');

$feat = $feat->transform('supercontig');
debug('--- transforming to supercontig coords ---');
ok($feat);
print_features([$feat]);


#
# Test fetch_by_Slice_and_pid
#
$feats = $pafa->fetch_all_by_Slice_and_pid($chr_slice, '90');
debug('--- fetching by chr Slice and pid (90) ---');
debug("Got " . scalar(@$feats));
ok(@$feats == 64);
print_features($feats);

#
# Test store
#

$multi->save('core', 'protein_align_feature');

my $analysis   = $feat->analysis;
my $slice      =
  $db->get_SliceAdaptor->fetch_by_region('contig','AL031658.11.1.162976');
my $start      = 100;
my $end        = 200;
my $strand     = 1;
my $hstart     = 10;
my $hend       = 110;
my $hstrand    = -1;
my $hseqname   = 'test';
my $score      = 80;
my $percent_id = 90;
my $evalue     = 23.2;
my $cigar_string = '100M';
my $hcoverage  = 99.5;
my $external_db_id = 2200;

$feat = Bio::EnsEMBL::DnaPepAlignFeature->new
  (-START  => $start,
   -END    => $end,
   -STRAND => $strand,
   -SLICE  => $slice,
   -HSTART => $hstart,
   -HEND   => $hend,
   -HSTRAND => 1,
   -HSEQNAME => $hseqname,
   -CIGAR_STRING => '100M',
   -PERCENT_ID => $percent_id,
   -SCORE    => $score,
   -P_VALUE => $evalue,
   -ANALYSIS => $analysis,
   -HCOVERAGE => $hcoverage,
   -EXTERNAL_DB_ID => $external_db_id );

ok(!$feat->is_stored($db));

$pafa->store($feat);

ok($feat->is_stored($db));

my $dbID = $feat->dbID();

$feat = $pafa->fetch_by_dbID($dbID);

ok($feat->dbID == $dbID);
ok($feat->start == $start);
ok($feat->end  == $end);
ok($feat->strand == $strand);
ok($feat->slice->name eq $slice->name);
ok($feat->hstart == $hstart);
ok($feat->hend   == $hend);
ok($feat->hseqname eq $hseqname);
ok($feat->cigar_string eq $cigar_string);
ok($feat->percent_id == $percent_id);
ok($feat->score == $score);
ok($feat->p_value == $evalue);
ok($feat->analysis->logic_name eq $analysis->logic_name);
ok($feat->external_db_id == $external_db_id);
ok($feat->hcoverage == $hcoverage);

$multi->restore('core', 'protein_align_feature');




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

done_testing();