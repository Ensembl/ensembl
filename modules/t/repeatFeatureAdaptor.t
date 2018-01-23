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
use warnings;

use Test::More;
use Test::Warnings;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0;


my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $db = $multi->get_DBAdaptor( 'core' );

my $slice_adaptor = $db->get_SliceAdaptor();
my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name('RepeatMask');
my $slice = $slice_adaptor->fetch_by_seq_region_id(319456);

my $repeat_f_ad = $db->get_RepeatFeatureAdaptor();
my $repeat_c_ad = $db->get_RepeatConsensusAdaptor();


#
# Test storing and simple retrieval
#

$multi->hide( "core", "repeat_feature", "repeat_consensus" );
my $repeat_consensus = Bio::EnsEMBL::RepeatConsensus->new();

$repeat_consensus->length(10);
$repeat_consensus->repeat_class('dummy');
$repeat_consensus->name('dummy');
$repeat_consensus->repeat_consensus('ATGCATGCAT');

ok($repeat_consensus);

$repeat_c_ad->store($repeat_consensus);

ok(1);

my $repeat_feature = Bio::EnsEMBL::RepeatFeature->new();

$repeat_feature->start(26);
$repeat_feature->end(65);
$repeat_feature->strand(1);
$repeat_feature->hstart(6);
$repeat_feature->hend(45);
$repeat_feature->score(100);
$repeat_feature->analysis($analysis);
$repeat_feature->repeat_consensus($repeat_consensus);
$repeat_feature->slice( $slice );

ok($repeat_feature);

$repeat_f_ad->store( $repeat_feature );


ok(1);

my $repeats = $repeat_f_ad->fetch_all_by_Slice($slice);

my $repeat = $repeats->[0];

ok($repeat);

ok($repeat->start == 26);
ok($repeat->hend == 45);

my $dbID = $repeat->dbID;

my $r = $repeat_f_ad->fetch_by_dbID($dbID);

ok($r->dbID == $dbID && $r->start == 26 && $r->hend == 45);


#
# Test  list_dbIDs
#
my $ids = $repeat_f_ad->list_dbIDs();
ok (@{$ids});

$multi->restore('core' );



#
# Test retrieval via Slice
#
$slice = $db->get_SliceAdaptor->fetch_by_region('chromosome',
                                                '20', 30_000_000,
                                                40_000_000);

my $feats = $repeat_f_ad->fetch_all_by_Slice($slice);
debug('fetching by chromosomal slice---');
debug("Got " . scalar(@$feats) . " features back");
# print_features($feats);

$feats = $repeat_f_ad->fetch_all_by_Slice( $slice, undef, "LTR" );
debug( "fetching by type LTR" );
debug( "Got ".scalar( @$feats ). " back" );

{
  my $f = $repeat_f_ad->fetch_all_by_Slice( $slice, undef, ["LTR", 'LTR'] );
  is(scalar(@{$f}), 126, 'Checking LTRs returned if querying by array');
}

# print_features( $feats );

$r = $repeat_f_ad->fetch_by_dbID(518841);
$r = $r->transform('supercontig');
debug('---fetching by dbID and transform to supercontig coords---');
# print_features([$r]);

done_testing();



sub print_features {
  my $features = shift;
  foreach my $f (@$features) {
    if(defined($f)) {
      my $seqname = $f->slice->seq_region_name();
      my $analysis = $f->analysis->logic_name();
      debug($seqname . ' ' . $f->start().'-'.$f->end().'('.$f->strand().
            ') ['. $f->dbID.'] '. $f->repeat_consensus->name() . ' ' .
            $f->hstart .'-'.$f->hend.' '.$f->score() .
            " ($analysis)");
    } else {
      debug('UNDEF');
    }
  }
}
