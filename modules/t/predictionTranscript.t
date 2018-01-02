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
use Bio::EnsEMBL::PredictionTranscript;

our $verbose = 0;

#
# 1 PredictionTranscript compiles
#
ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $db = $multi->get_DBAdaptor( "core" );

my $sa = $db->get_SliceAdaptor();

my $slice = $sa->fetch_by_region("chromosome", "20", 30_252_000, 31_252_001 );

my $p_transs = $slice->get_all_PredictionTranscripts();

#
# 2 Verify prediciton transcripts can be obtained
#
debug( "Retrieved ".scalar( @$p_transs )." PredictionTranscripts." );
ok( scalar( @$p_transs ) );


foreach my $pt (@$p_transs) {
  debug( "Stable_id: ".$pt->stable_id()." Start: ".$pt->start." End: ".$pt->end() , "\n" );
}

my $pt = $p_transs->[0];

my $exons = $pt->get_all_Exons;

#
# 3 test get all exons
#
debug( "First PT had ".scalar( @$exons ). " exons." );
ok(scalar @$exons);

#
#  4 test new
#
my $new_pt = new Bio::EnsEMBL::PredictionTranscript( -exons => $exons,
                                                    -display_label => 'test');
ok(scalar @{$new_pt->get_all_Exons});
ok($new_pt->display_label() eq $new_pt->stable_id() && 
   $new_pt->display_label eq 'test');

#
# 5 test stable_id
#
debug( "stable_id: ".$pt->stable_id() );
ok($pt->stable_id =~ /.*/ );


#
# 6 test start
#
ok(&Bio::EnsEMBL::Test::TestUtils::test_getter_setter($pt, 'start', 8));

#
# 7 test end
#
ok(&Bio::EnsEMBL::Test::TestUtils::test_getter_setter($pt, 'end', 9));


#
# 8 test analysis 
#
my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name('Vertrna');
ok(&Bio::EnsEMBL::Test::TestUtils::test_getter_setter($pt, 'analysis', $analysis));

#
# 9 test dbID
#
ok(&Bio::EnsEMBL::Test::TestUtils::test_getter_setter($pt, 'dbID', 11));

#
# 10 test adaptor
#
my $pta = $db->get_PredictionTranscriptAdaptor;
ok(&Bio::EnsEMBL::Test::TestUtils::test_getter_setter($pt, 'adaptor', $pta));

#
# 11-14 test add Exon
#
$pt = new Bio::EnsEMBL::PredictionTranscript();

my $exon = new Bio::EnsEMBL::PredictionExon;
$exon->start(40);
$exon->end(50);
$exon->slice( $slice );
$exon->strand( 1 );
$pt->add_Exon($exon);

$exon = new Bio::EnsEMBL::PredictionExon;
$exon->start(20);
$exon->end(30);
$exon->slice( $slice );
$exon->strand( 1 );
$pt->add_Exon($exon);

$exon = new Bio::EnsEMBL::PredictionExon;
$exon->start( 1 );
$exon->end(10);
$exon->slice( $slice );
$exon->strand( 1 );
$pt->add_Exon($exon);


#check that transcript start + end updated
ok( $pt->end() == 50 );
ok( $pt->start() == 1 );

my $all_exons = $pt->get_all_Exons();
ok( $all_exons->[0]->start() == 1 );
ok( $all_exons->[2]->end() == 50 );

#
# 15-18 -1 strand checks for add_Exon
#

$pt = new Bio::EnsEMBL::PredictionTranscript();

$exon = new Bio::EnsEMBL::PredictionExon;
$exon->start(40);
$exon->end(50);
$exon->slice( $slice );
$exon->strand( -1 );
$pt->add_Exon($exon);

$exon = new Bio::EnsEMBL::PredictionExon;
$exon->start( 1 );
$exon->end(10);
$exon->slice( $slice );
$exon->strand( -1 );
$pt->add_Exon($exon);

$exon = new Bio::EnsEMBL::PredictionExon;
$exon->start(20);
$exon->end(30);
$exon->slice( $slice );
$exon->strand( -1 );
$pt->add_Exon($exon);

ok( $pt->end() == 50 );
ok( $pt->start() == 1 );

$all_exons = $pt->get_all_Exons();
ok( $all_exons->[0]->start() == 40 );
ok( $all_exons->[2]->end() == 10 );






#
# 19-23 test flush exons
#
debug( "Flush exons effect test" );
$pt->flush_Exons;
ok(scalar @{$pt->get_all_Exons} == 0);
debug( "pt->start ".($pt->start()||"undef") );
ok(!defined $pt->start);
debug( "pt->end ".($pt->end()||"undef") );
ok(!defined $pt->end);
debug( "pt->coding_start ".($pt->coding_region_start()||"undef") );
ok(!defined $pt->coding_region_start);
debug( "pt->coding_end ".($pt->coding_region_end()||"undef") );
ok(!defined $pt->coding_region_end);


#
# 24 test get_all_translateable_Exons
#
$pt = $p_transs->[0];

ok(scalar @{$pt->get_all_translateable_Exons} == scalar @{$pt->get_all_Exons});


#
# 25 test length
#
my $len = 0;
foreach my $ex (@{$pt->get_all_Exons}) {
  $len += $ex->length();
}
ok($len == $pt->length);

#
# 26 test translate
#
my $translated = $pt->translate->seq();
debug( "Translated sequence: $translated" );
ok( $translated );

#
# 27 test spliced_seq()
#
my $spliced_seq = $pt->spliced_seq();
debug( "Spliced seq: ".$spliced_seq );
ok( $spliced_seq );

#
# 28 test pep2genomic
#
my $pend = $pt->length / 3;
my $pstart = 2;

my $defined_exons_count = 0;
foreach my $e (@{$pt->get_all_Exons}) {
    $defined_exons_count++;
}
# should return genomic coords for each exon since covers entire peptide
ok($defined_exons_count == $pt->pep2genomic($pstart, $pend));


#
# 29 test cdna2genomic
#

ok($defined_exons_count == $pt->cdna2genomic( 1, $pt->length()));


#
# 30 test type
#
ok(&Bio::EnsEMBL::Test::TestUtils::test_getter_setter($pt, 'biotype', 'test'));

#
# 31 test fetch_by_stable_id
#

my $stable_id = 'GENSCAN00000011401';

$pt = $pta->fetch_by_stable_id($stable_id);
ok($pt->stable_id() eq $stable_id);

# 32 list_dbIDs
my $ids = $pta->list_dbIDs();
ok (@{$ids});


ok($pt->display_id eq $pt->stable_id);


#
# test PredictionTranscriptAdaptor::store
#

my $analysis_adaptor = $db->get_AnalysisAdaptor();

$analysis = $analysis_adaptor->fetch_by_logic_name('Genscan');

$multi->hide('core', 'prediction_transcript', 'prediction_exon');

my @exons;
push @exons, Bio::EnsEMBL::PredictionExon->new
  (-START  => 100,
   -END    => 200,
   -STRAND => 1,
   -P_VALUE => 0.98,
   -SCORE  => 50,
   -SLICE  => $slice,
   -PHASE  => 0);

push @exons, Bio::EnsEMBL::PredictionExon->new
  (-START  => 300,
   -END    => 400,
   -STRAND => 1,
   -P_VALUE => 0.99,
   -SCORE  => 75,
   -SLICE  => $slice,
   -PHASE  => $exons[0]->length % 3);
   
$pt = Bio::EnsEMBL::PredictionTranscript->new
  (-EXONS => \@exons,
   -SLICE => $slice,
   -ANALYSIS => $analysis);

$pta->store($pt);

$pt = $pta->fetch_by_dbID($pt->dbID);

ok($pt && $pt->dbID);

# check that display label is automatically generated if it was unspecified
ok($pt->display_label() eq $pt->stable_id() 
   && $pt->display_label() eq 'GENSCAN' . ('0'x(11-length($pt->dbID()))) . $pt->dbID());

@exons = @{$pt->get_all_Exons};
ok(@exons == 2);

ok($exons[0]->start == $slice->start + 100 - 1);
ok($exons[1]->start == $slice->start + 300 - 1);

ok($exons[0]->score == 50);
ok($exons[1]->score == 75);

$multi->restore();

done_testing();
