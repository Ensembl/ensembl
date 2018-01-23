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
use Test::Exception;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

my $db = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $db->get_DBAdaptor('core');

require_ok('Bio::EnsEMBL::Exon');
require_ok('Bio::EnsEMBL::Intron');
require_ok('Bio::EnsEMBL::IntronSupportingEvidence');
require_ok('Bio::EnsEMBL::Transcript');

my $feature_id = sub {
  my ($f) = @_;
  no warnings 'uninitialized';
  return join(q{|}, 
    $f->seq_region_name(), $f->start, $f->end(), $f->strand(),
    ($f->dbID()||q{?})
  ); 
};

my @basic_ise_args = (
  -HIT_NAME => 'wibble',
  -SCORE_TYPE => 'DEPTH',
  -SCORE => 3.000,
  -ANALYSIS => $dba->get_AnalysisAdaptor()->fetch_by_dbID(1282),,
  -IS_SPLICE_CANONICAL => 0
);

my $ise_adaptor = $dba->get_IntronSupportingEvidenceAdaptor();
my $transcript_adaptor = $dba->get_TranscriptAdaptor();
my $exon_adaptor = $dba->get_ExonAdaptor();
my @intron_tables = ('intron_supporting_evidence', 'transcript_intron_supporting_evidence'); 
  
my $assert_ise_vs_intron = sub {
  my ($transcript_id, $start, $end, $strand) = @_;
  
  $db->hide('core', @intron_tables);
    
  my $t = $transcript_adaptor->fetch_by_dbID($transcript_id);
  my $exons = $t->get_all_Exons();
  my ($e1, $e2) = ($exons->[0], $exons->[1]);
  my $intron_from_exons = Bio::EnsEMBL::Intron->new($e1, $e2);
  
  note 'Asserting Intron tests';
  warns_like 
    { Bio::EnsEMBL::Intron->new($e1, $exon_adaptor->fetch_by_dbID(162033)); }
    qr/Exons have different slice references/, 
    'Intron must warn if a different reference slice is used';

  ok($intron_from_exons->is_splice_canonical(), 'Checking Intron is canonical in its splicing');
  
  note 'Starting IntronSupportingEvidence tests';
  
  my $ise = Bio::EnsEMBL::IntronSupportingEvidence->new(
    -SLICE => $t->slice(),
    -STRAND => $strand,
    -START => $start,
    -END => $end,
    @basic_ise_args,
  );
  my $ise_from_intron = Bio::EnsEMBL::IntronSupportingEvidence->new(
    -INTRON => $intron_from_exons, @basic_ise_args,
  );
  
  ok($ise_from_intron->is_splice_canonical(), 'Checking IntronSupportingEvidence is canonical in its splicing');
  is($feature_id->($intron_from_exons), $feature_id->($ise_from_intron), 'IntronSupportingEvidence returns the equivalent Intron object');
  
  throws_ok { $ise->score_type('RUBBISH') } qr/not allowed/, 'IntronSupportingEvidence refuses unsupported ENUMs';
  
  throws_ok { $ise_adaptor->update($ise) } qr/Cannot/, 'Checking IntronSupportingEvidenceAdaptor refuses to update a non-persisted object';
  
  note 'Checking our find Exon methods operate';
  ok($e1 eq $ise_from_intron->find_previous_Exon($t), 'Checking find_previous_exon() finds exon 1 (5 prime) by coordinates');
  ok($e2 eq $ise_from_intron->find_next_Exon($t), 'Checking find_next_exon() finds exon 2 (3 prime) by coordinates');
  
  my $strandedness = ($strand == 1) ? '+ve' : '-ve';
  is($feature_id->($ise_from_intron), $feature_id->($ise), "Checking both ISEs on the $strandedness strand are equivalent");
  
  is($t->add_IntronSupportingEvidence($ise), 1, 'First addition of duplciate ISE works');
  is($t->add_IntronSupportingEvidence($ise_from_intron), 0, 'Second addition of duplciate ISE fails to add');
  
  is(scalar(@{$t->get_all_IntronSupportingEvidence()}), 1, 'Checking Transcript has only 1 ISE');
  
  note('Storing and removing linkage before storing the IntronSupportingEvidence object');
  throws_ok { $ise_adaptor->remove($ise) } qr/Cannot delete .+ supporting evidence/, 'Checking API refuses to remove a non-persisted object';
  throws_ok { $ise_adaptor->store_transcript_linkage($ise, $t) } qr/Cannot perform the link/, 'Checking API refuses to store a non-persisted linkage';
  throws_ok { $ise_adaptor->remove_transcript_linkage($ise, $t) } qr/evidence.+has not/, 'Checking API refuses to remove if the object has not been persisted';
  $t->{dbID} = undef;
  $t->{adaptor} = undef;
  throws_ok { $ise_adaptor->remove_transcript_linkage($ise, $t) } qr/transcript.+has not/, 'Checking API refuses to remove if the object has not been persisted';
  $t->{dbID} = $transcript_id;
  $t->{adaptor} = $dba->get_TranscriptAdaptor();
  
  note('Multiple calls to store for IntronSupportingEvidence and Transcript');
  $ise_adaptor->store($ise_from_intron);
  is($ise_adaptor->store($ise_from_intron), $ise_from_intron->dbID(), 'Same ID should be returned');
  $ise_adaptor->store($ise);
  $ise_adaptor->store_transcript_linkage($ise, $t);
  $ise_adaptor->store_transcript_linkage($ise_from_intron, $t);
  $ise_adaptor->store_transcript_linkage($ise, $t);
  $ise_adaptor->store_transcript_linkage($ise, $t);
  
  #Slice searches
  note 'Performing Slice based searches';
  my $assert_slice_search = sub {
    my ($start, $end, $count) = @_;
    my $slice = $dba->get_SliceAdaptor()->fetch_by_region('toplevel', 20, $start, $end);
    my $features = $ise_adaptor->fetch_all_by_Slice($slice);
    is(scalar(@{$features}), $count, sprintf(q{Checking region '%s'}, $slice->name()));
  };
  
  $assert_slice_search->($start, $end, 1);
  $assert_slice_search->($start-6000, $start-3000, 0);
  $assert_slice_search->($end+3000, $end+6000, 0);
  $assert_slice_search->($start-3000, $end+3000, 1);
  
  is($ise->dbID(), $ise_from_intron->dbID(), 'Checking both IntronSupportingEvidence objects are given the same dbID');
  is_rows(1, $dba, 'intron_supporting_evidence', 'where seq_region_start =?', [$ise->start()]);
  is_rows(1, $dba, 'transcript_intron_supporting_evidence', 'where intron_supporting_evidence_id =?', [$ise->dbID()]);
  is_rows(1, $dba, 'transcript_intron_supporting_evidence', 'where previous_exon_id =? and next_exon_id =?', [$e1->dbID(), $e2->dbID()]);
  
  $ise_adaptor->remove_transcript_linkage($ise, $t);
  is_rows(0, $dba, 'transcript_intron_supporting_evidence', 'where intron_supporting_evidence_id =?', [$ise->dbID()]);
  $ise_adaptor->store_transcript_linkage($ise_from_intron, $t);
  is_rows(1, $dba, 'transcript_intron_supporting_evidence', 'where intron_supporting_evidence_id =?', [$ise->dbID()]);
  
  note 'Checking our find Exon methods operate with DBIDs';
  ok($e1 eq $ise_from_intron->find_previous_Exon($t), 'Checking find_previous_exon() finds exon 1 (5 prime) with DBIDs');
  ok($e2 eq $ise_from_intron->find_next_Exon($t), 'Checking find_next_exon() finds exon 2 (3 prime) DBIDs');
  
  note 'Fetching a new Transcript and getting the linked IntronSupportingEvidence';
  my $new_transcript = $transcript_adaptor->fetch_by_dbID($t->dbID());
  my $new_evidence = $new_transcript->get_all_IntronSupportingEvidence();
  is(scalar(@{$new_evidence}), 1, 'Checking we have 1 ISE');
  is($new_evidence->[0]->score()*1, $ise->score(), 'Checking scores are the same');
  ok($new_evidence->[0]->equals($ise), 'Checking retrieved feature is the same');
  
  note 'Using a Transcript with no official links to the ISE';
  my $bad_transcript = $transcript_adaptor->fetch_by_dbID(21727);
  ok(!defined $ise->find_next_Exon($bad_transcript), 'Checking we cannot find the flanking exon with an incorrect exon');
  
  note 'Updating IntronSupportingEvidence';
  $ise->hit_name('Woobly');
  $ise_adaptor->update($ise);
  is_rows(1, $dba, 'intron_supporting_evidence', 'where hit_name=?', [$ise->hit_name()]);
  
  note('Removing IntronSupportingEvidence but still has Transcripts attached');
  throws_ok { $ise_adaptor->remove($ise) } qr/transcripts attached/, 'Cannot remove ISE if we have not detached all Transcripts';
  
  note('Removing non-existent transcript linkage');
  $ise_adaptor->remove_transcript_linkage($ise, $transcript_adaptor->fetch_by_dbID(21717));
  is_rows(1, $dba, 'transcript_intron_supporting_evidence', 'where intron_supporting_evidence_id =?', [$ise->dbID()]);
  
  note('Removing all transcript linkages');
  $ise_adaptor->remove_all_transcript_linkages($ise);
  is_rows(0, $dba, 'transcript_intron_supporting_evidence', 'where intron_supporting_evidence_id =?', [$ise->dbID()]);
  
  note('Removed IntronSupportingEvidence');
  $ise_adaptor->remove($ise);
  is_rows(0, $dba, 'intron_supporting_evidence', 'where seq_region_start =?', [$ise->start()]);
  $db->restore('core', @intron_tables);
};

############# Basic tests

#Forward work
$assert_ise_vs_intron->(21716, 30274426, 30284450, 1);
#Reverse work
$assert_ise_vs_intron->(21719, 30326248, 30327734, -1);

############# Transcript centric tests (note only written on forward Transcripts; prior tests cover reverse strands)

{
  note 'Starting Transcript based tests. Hiding transcript and intron tables. Saving exon state';
  $db->hide('core', qw/transcript exon_transcript transcript_attrib transcript_supporting_feature meta_coord/, @intron_tables);
  $db->save('core', qw/exon dna_align_feature protein_align_feature/);
  
  my $gene = $dba->get_GeneAdaptor()->fetch_by_dbID(18258);
  my $slice = $gene->slice();
  my $analysis = $gene->analysis();
  my @exons = map { $exon_adaptor->fetch_by_dbID($_) } (161880, 161882, 161879);
  $_->slice($slice) for @exons;
  my $intron_one = Bio::EnsEMBL::Intron->new(@exons[0..1]);
  my $ise_one = Bio::EnsEMBL::IntronSupportingEvidence->new(
    -INTRON => $intron_one, -HIT_NAME => 'transwib', -SCORE => '1', -SCORE_TYPE => 'DEPTH', -ANALYSIS => $analysis
  );
  my $intron_two = Bio::EnsEMBL::Intron->new(@exons[1..2]);
  my $transcript_one = Bio::EnsEMBL::Transcript->new(-EXONS => \@exons, -ANALYSIS => $analysis, -SLICE => $slice);
  $transcript_one->add_IntronSupportingEvidence($ise_one);
  $transcript_one->add_IntronSupportingEvidence(Bio::EnsEMBL::IntronSupportingEvidence->new(
    -INTRON => $intron_two, -HIT_NAME => 'transwob', -SCORE => '2', -SCORE_TYPE => 'DEPTH', -ANALYSIS => $analysis
  ));
  $transcript_one->add_IntronSupportingEvidence($ise_one);
  is(scalar @{$transcript_one->get_all_IntronSupportingEvidence()}, 2, 'Checking additonal adds of the same data are ignored');
  
  #Transcript two shares an Intron with Transcript one
  my $transcript_two = Bio::EnsEMBL::Transcript->new(-EXONS => [@exons[0..1]], -ANALYSIS => $analysis, -SLICE => $slice);
  $transcript_two->add_IntronSupportingEvidence(Bio::EnsEMBL::IntronSupportingEvidence->new(
    -INTRON => $intron_one,
    -HIT_NAME => 'transwib', -SCORE => '1', -SCORE_TYPE => 'DEPTH', -ANALYSIS => $analysis
  ));
  
  my $assert_table_counts = sub {
    my ($t_count, $et_count, $ise_count, $tise_count) = @_;
    is_rows($t_count, $dba, 'transcript');
    is_rows($et_count, $dba, 'exon_transcript');
    is_rows($ise_count, $dba, 'intron_supporting_evidence');
    is_rows($tise_count, $dba, 'transcript_intron_supporting_evidence');
  };
  
  $transcript_adaptor->store($transcript_one);
  $transcript_adaptor->store($transcript_two);
  note 'Stored both transcripts';
  $assert_table_counts->(2,5,2,3);
  
  #Start removing
  $transcript_adaptor->remove($transcript_one);
  note 'Removed transcript one';
  $assert_table_counts->(1,2,1,1);
  note 'Removed transcript two';
  $transcript_adaptor->remove($transcript_two);
  $assert_table_counts->(0,0,0,0);
  
  $db->restore('core');
  note 'All tables restored';
  $assert_table_counts->(27,173,0,0);
}

############# Feature transformation (asserting what occurs when we try to project to something else)
{
  my $transcript = $transcript_adaptor->fetch_by_dbID(21716);
  my $intron = Bio::EnsEMBL::Intron->new(@{$transcript->get_all_Exons()}[0..1]);
  my $ise = Bio::EnsEMBL::IntronSupportingEvidence->new(
   -INTRON => $intron, -HIT_NAME => 'transmogrification', -SCORE => '5', -SCORE_TYPE => 'DEPTH', -ANALYSIS => $transcript->analysis()
  );
  $transcript->add_IntronSupportingEvidence($ise);
  
  my $assert_trans = sub {
    my ($t, $type) = @_;
    note 'Attempted a '.$type.' at the Transcript level';
    my @ise = @{$t->get_all_IntronSupportingEvidence()};
    is(scalar(@ise), 1, 'Still only has 1 IntronSupportingEvidence');
    is($ise[0]->start(), 10912, 'Checking '.$type.' start');
    is($ise[0]->end(), 20936, 'Checking '.$type.' end');
    is($ise[0]->strand(), 1, 'Checking '.$type.' strand');
  };
  
  $assert_trans->($transcript->transform('seqlevel'), 'transform');
  $assert_trans->($transcript->transfer($dba->get_SliceAdaptor()->fetch_by_region('contig', 'AL031658.11.1.162976')), 'transfer');
}

############# Data movement tests (asserting what occurs when we migrate data from one DB to another)
{
  my $new_dba = $db->get_DBAdaptor('empty');
  $db->hide('core', @intron_tables);
  
  my $empty_obj = sub {
    $_[0]->dbID(undef);
    $_[0]->adaptor(undef);
    return;
  };
  
  my $assert_table_counts = sub {
    my ($g_count, $t_count, $et_count, $e_count, $ise_count, $tise_count) = @_;
    is_rows($g_count, $new_dba, 'gene');
    is_rows($t_count, $new_dba, 'transcript');
    is_rows($et_count, $new_dba, 'exon_transcript');
    is_rows($e_count, $new_dba, 'exon');
    is_rows($ise_count, $new_dba, 'intron_supporting_evidence');
    is_rows($tise_count, $new_dba, 'transcript_intron_supporting_evidence');
  };
  
  $db->hide('empty', @intron_tables, qw/gene transcript exon_transcript exon/);
  
  $assert_table_counts->(0,0,0,0,0,0);
  
  my $gene = $dba->get_GeneAdaptor->fetch_by_dbID(18256);
  my ($transcript) = grep { $_->dbID() == 21716 } @{$gene->get_all_Transcripts()};
  my $intron = Bio::EnsEMBL::Intron->new(@{$transcript->get_all_Exons()}[0..1]);
  my $ise = Bio::EnsEMBL::IntronSupportingEvidence->new(
   -INTRON => $intron, -HIT_NAME => 'datamovement', -SCORE => '1', -SCORE_TYPE => 'DEPTH', -ANALYSIS => $transcript->analysis()
  );
  $transcript->add_IntronSupportingEvidence($ise);
  $ise_adaptor->store($ise);
  $ise_adaptor->store_transcript_linkage($ise, $transcript);
   
  foreach my $t (@{$gene->get_all_Transcripts()}) {
    foreach my $e (@{$t->get_all_Exons()}) {
      $empty_obj->($e);
    }
    $empty_obj->($t);
    $t->display_xref(undef);
  }
  $empty_obj->($gene);
  $gene->display_xref(undef);
  
  $new_dba->get_GeneAdaptor()->store($gene);
  note 'Checking all objects have been copied across as expected';
  $assert_table_counts->(1,2,10,9,1,1);
  $db->restore('empty');
  $db->restore('core');
}

done_testing(); 
