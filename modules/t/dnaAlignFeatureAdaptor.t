# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use Test::More;
use Test::Warnings qw(warning);
use Bio::EnsEMBL::Test::MultiTestDB;

use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
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
is(@$feats, 6188, "Found 6188 features");
print_features($feats);


my $ctg_slice;
$ctg_slice  = $db->get_SliceAdaptor->fetch_by_region('contig',
                                                     'AL031658.11.1.162976',
                                                     1,
                                                     50_000);
$feats = $dafa->fetch_all_by_Slice($ctg_slice);
debug('--- contig AL031658.11.1.162976 (1-50000) features ---');
debug("Got " . scalar(@$feats));
is(@$feats, 709, "found 709 features");
print_features($feats);


#
# Test fetch_by_dbID
#
my $feat = $dafa->fetch_by_dbID(22171863);
debug('--- fetching by dbID ---');
ok($feat);
print_features([$feat]);
is($feat->db_name, 'EMBL', "Correct feature db_name");
is($feat->db_display_name, 'EMBL', "Correct feature db_display_name");

$feat = $feat->transform('supercontig');
debug('--- transform to supercontig coords ---');
ok($feat);
print_features([$feat]);


#
# Test fetch_by_Slice_and_pid
#
$feats = $dafa->fetch_all_by_Slice_and_pid($chr_slice, '90');
debug('--- fetching by chr Slice and pid (90) ---');
debug("Got " . scalar(@$feats));
is(@$feats, 452, "Found 452 features");
#print_features($feats);

#
# Test store
#
$multi->save('core', 'dna_align_feature', 'attrib_type', 'dna_align_feature_attrib');

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
my $align_type = 'ensembl';
my $hcoverage  = 99.5;
my $external_db_id = 2200;

warning {
$feat = Bio::EnsEMBL::DnaDnaAlignFeature->new
  (-START  => $start,
   -END    => $end,
   -STRAND => $strand,
   -SLICE  => $slice,
   -HSTART => $hstart,
   -HEND   => $hend,
   -HSTRAND => $hstrand,
   -HSEQNAME => $hseqname,
   -CIGAR_STRING => '100M',
   -PERCENT_ID => $percent_id,
   -SCORE    => $score,
   -P_VALUE => $evalue,
   -ANALYSIS => $analysis,
   -HCOVERAGE => $hcoverage,
   -EXTERNAL_DB_ID => $external_db_id);
};

my $attrib = Bio::EnsEMBL::Attribute->new(
  -NAME        => 'test_daf_name',
  -CODE        => 'test_daf_code',
  -DESCRIPTION => 'test_daf_desc',
  -VALUE       => 'test_daf_value');

$feat->add_Attributes($attrib);

ok(!$feat->is_stored($db));

$dafa->store($feat);

ok($feat->is_stored($db));

my $dbID = $feat->dbID();

$feat = $dafa->fetch_by_dbID($dbID, 'contig');

is($feat->dbID, $dbID, "Correct dbID");
is($feat->start, $start, "Correct start");
is($feat->end, $end, "Correct end");
is($feat->strand, $strand, "Correct strand");
is($feat->slice->name, $slice->name, "Correct name");
is($feat->hstart, $hstart, "Correct hstart");
is($feat->hend, $hend, "Correct hend");
is($feat->hseqname, $hseqname, "Correct hseqname");
is($feat->cigar_string, $cigar_string, "Correct cigar string");
is($feat->percent_id, $percent_id, "Correct percent id");
is($feat->score, $score, "Correct score");
is(sprintf("%.6f", $feat->p_value), sprintf("%.6f", $evalue), "Correct evalue");
is($feat->analysis->logic_name, $analysis->logic_name, "Correct logic_name");
is($feat->external_db_id, $external_db_id, "Correct external_db_id");
is($feat->hcoverage, $hcoverage, "Correct hcoverage");

my @attribs = @{$feat->get_all_Attributes};
is(@attribs, 1, "Stored 1 attribute");

my $attrib = $attribs[0];
ok($attrib->code eq 'test_daf_code', 'Stored attribute');

$dafa->remove($feat);

is_rows(0, $db, "dna_align_feature", "where dna_align_feature_id = ? ", [$dbID]);
is_rows(0, $db, "dna_align_feature_attrib", "where dna_align_feature_id = ? ", [$dbID]);


$multi->restore('core', 'dna_align_feature', 'attrib_type', 'dna_align_feature_attrib');


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
