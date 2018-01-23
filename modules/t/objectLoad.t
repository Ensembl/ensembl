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
#use Test::Differences qw/eq_or_diff/;
use Bio::EnsEMBL::Test::MultiTestDB;

my $db = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $db->get_DBAdaptor('core');

$db->save('core', qw/intron_supporting_evidence transcript_intron_supporting_evidence/);

my $ga = $dba->get_GeneAdaptor();

#Test the keys held by a fully populated gene change
my $gene = $ga->fetch_by_stable_id('ENSG00000131044');
$gene->load();
my @generic_keys = qw/
  dbID adaptor stable_id created_date modified_date version
/;
my @location_keys = qw/start end strand slice/;

my $expected_keys = {
  'gene' => [sort @generic_keys, @location_keys, qw/
    analysis
    attributes
    is_current biotype source
    description
    canonical_transcript canonical_transcript_id _transcript_array
    dbentries display_xref external_name external_db external_status
  /],
  'transcript' => [sort @generic_keys, @location_keys, qw/
    analysis
    attributes
    _ise_array
    _supporting_evidence
    _trans_exon_array
    
    cdna_coding_start cdna_coding_end
    alternative_translations translation edits_enabled
    
    is_current biotype description source
    
    dbentries display_xref external_name external_db external_status external_display_name
  /],
  'translation' => [sort @generic_keys, qw/
    attributes
    start_exon end_exon start end transcript
    seq
    protein_features
    dbentries
  /],
  'exon' => [sort @generic_keys, @location_keys, qw/
    _seq_cache
    _supporting_evidence
    phase end_phase is_constitutive is_current 
  /],
};
$expected_keys->{transcript_canonical} = [sort @{$expected_keys->{transcript}}, 'is_canonical'];

assert_keys($gene, 'gene', 'Fully populated gene should have all elements filled');
my $transcripts = $gene->get_all_Transcripts();
is(scalar(@{$transcripts}), 2, 'Gene has 2 transcripts');
assert_transcript($transcripts->[0], 'transcript_canonical');
assert_transcript($transcripts->[1], 'transcript');

$db->restore();

done_testing();

sub assert_transcript {
  my ($transcript, $key) = @_;
  assert_keys($transcript, $key, 'Fully populated transcript');
  assert_keys($transcript->translation(), 'translation', 'Fully populated translation');
  foreach my $exon (@{$transcript->get_all_Exons()}) {
    assert_keys($exon, 'exon', 'Fully populated exon');
  }
}

sub assert_keys {
  my ($feature, $type_key, $msg) = @_;
  is_deeply(get_keys($feature), $expected_keys->{$type_key}, $msg);
#  eq_or_diff(get_keys($feature), $expected_keys->{$type_key}, $msg);
}

sub get_keys {
  return [sort keys %{$_[0]}];
}
