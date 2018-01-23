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
use Test::Builder;
use Test::Exception;
use Test::MockObject;
use Test::MockObject::Extends;

use Bio::EnsEMBL::Utils::TranscriptSelector;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Test::TestUtils;

my $verbose = Test::Builder->new()->no_diag() ? 0 : 1;

my $transcript_selector;

warns_like 
  { $transcript_selector = Bio::EnsEMBL::Utils::TranscriptSelector->new(undef, $verbose) } 
  qr/Running without CCDS DB/, 'Checking we warn about the lack of a CCDS DB';

# Test the sorting algorithm 
# encoded arrays as follows:
#[transcript dbID, source , biotype, translation length, transcript length, stable ID]
my $sortables = [
    [qw( a 1 1 1 500 250 ENST7 )],
    [qw( b 1 2 1 500 250 ENST6 )],
    [qw( c 1 3 1 500 250 ENST5 )],
    [qw( d 1 1 1 450 250 ENST4 )],
    [qw( e 1 1 1 500 250 ENST3 )],
    [qw( f 0 3 3 0   700 ENST2 )],
    [qw( g 0 3 3 0   700 ENST1 )],
]; 

my $sorted = $transcript_selector->sort_into_canonical_order($sortables);

note "Sorted order";
note  join(',',@$sorted);

my $correct_order = [qw(e a d b c g f)];
is_deeply($sorted,$correct_order,'Canonical sort order');

$sortables = [];
$sorted = $transcript_selector->sort_into_canonical_order($sortables);
note join(',',@$sorted);
is(scalar(@$sorted), 0,'Null data into sort.');

# create a mock CCDS dba to test db-dependent code
# returns a pretend Slice object containing some test genes
# This tests check_Ens_trans_against_CCDS() and shortcuts
# the need for a slice adaptor.

my $mock_sa = Test::MockObject->new();
$mock_sa->set_isa('Bio::EnsEMBL::DBSQL::SliceAdaptor');
$mock_sa->mock('get_seq_region_id', sub {
   return 1; 
});

my $coord_system = Bio::EnsEMBL::CoordSystem->new(
    -NAME => 'landofgiants',
    -TOP_LEVEL => 0,
    -RANK => 1,
    -DBID => 1,
);

my $slice = Bio::EnsEMBL::Slice->new(
    -START => 1,
    -END => 10000,
    -STRAND => 1,
    -SEQ_REGION_LENGTH => 1e4,
    -SEQ_REGION_NAME => '1',
    -COORD_SYSTEM => $coord_system,
    -ADAPTOR => $mock_sa,
    -SEQ => 'N' x 10000,
);

my $other_slice = Bio::EnsEMBL::Slice->new(
    -START => 1,
    -END => 4000,
    -STRAND => -1,
    -SEQ_REGION_LENGTH => 4e3,
    -SEQ_REGION_NAME => '1',
    -COORD_SYSTEM => $coord_system,
    -ADAPTOR => $mock_sa,
    -SEQ => 'A' x 4000,
);

$slice = Test::MockObject::Extends->new($slice);
$slice->mock('get_all_Attributes', sub {return [];});

my $exon = Bio::EnsEMBL::Exon->new(
    -START => 1000,
    -END => 2000,
    -STRAND => 1,
    -DBID => 15,
    -STABLE_ID => 'ENSE01',
    -SLICE => $slice,
);
    
my $analysis = Bio::EnsEMBL::Analysis->new(
    -id => 1,
    -logic_name => 'bananas_are_nice',
);
    
my $transcript1 = Bio::EnsEMBL::Transcript->new(
    -DBID => 1,
    -STABLE_ID => 'ENST01',
    -BIOTYPE => 'protein_coding',
    -IS_CURRENT => 1,
    -SLICE => $slice,
    -ANALYSIS => $analysis,
);
$transcript1 = Test::MockObject::Extends->new($transcript1);
$transcript1->mock('translate', sub {return $slice});

my $transcript2 = Bio::EnsEMBL::Transcript->new(
    -DBID => 2,
    -STABLE_ID => 'ENST02',
    -BIOTYPE => 'nonsense_mediated_decay',
    -IS_CURRENT => 1,
    -SLICE => $slice,
    -ANALYSIS => $analysis,
);
$transcript2 = Test::MockObject::Extends->new($transcript2);
$transcript2->mock('translate', sub {return $slice});

my $transcript3 =  Bio::EnsEMBL::Transcript->new(
    -DBID => 3,
    -STABLE_ID => 'ENST03',
    -BIOTYPE => 'flying_poofish',
    -IS_CURRENT => 1,
    -SLICE => $slice,
    -ANALYSIS => $analysis,
);

$transcript3 = Test::MockObject::Extends->new($transcript3);
$transcript3->mock('translate', sub {return $other_slice});

$transcript2->add_Exon($exon);
$transcript2->translation(Bio::EnsEMBL::Translation->new(
    -START_EXON => $exon,
    -END_EXON   => $exon,
    -SEQ_START  => 1,
    -SEQ_END    => 1001,
    
  )
);
$transcript1->add_Exon($exon);
$transcript3->add_Exon($exon);

my $transcripts = [ $transcript1, $transcript2, $transcript3 ];

my $gene = Bio::EnsEMBL::Gene->new(
    -START  => 123,
    -END    => 2045,
    -STRAND => 1,
    -BIOTYPE => 'protein_coding',
    -TRANSCRIPTS => $transcripts,
    -SLICE => $slice,
    -STABLE_ID => 'CCDS01'
);
    
my $mock_slice = Test::MockObject->new();
$mock_slice->mock('get_seq_region_id', sub {
    return 1;   
});
$mock_slice->mock('is_circular', sub { return 0;});
$mock_slice->mock('get_all_Genes', sub {
    return [$gene];
});

$mock_sa->mock('fetch_by_region', sub {
    return $mock_slice;
});
$mock_sa->mock('is_reference', sub {
    return 1;
});
$mock_sa->mock('is_circular', sub { return 0;});
my $fake_seq_adaptor = Test::MockObject->new();
$fake_seq_adaptor->mock('fetch_by_Slice_start_end_strand', sub {my $seq = 'A'x 20; return \$seq;});
my $fake_db = Test::MockObject->new();
$fake_db->mock('get_SequenceAdaptor', sub {return $fake_seq_adaptor;});
$mock_sa->mock('db', sub {return $fake_db;});

my $mock_dba = Test::MockObject->new();
$mock_dba->mock('get_SliceAdaptor', sub {
    return $mock_sa;
});

$transcript_selector = Bio::EnsEMBL::Utils::TranscriptSelector->new($mock_dba, $verbose);
ok($transcript_selector->check_Ens_trans_against_CCDS($transcript2),'CCDS transcript lookup with good data');
ok($transcript_selector->check_Ens_trans_against_CCDS($transcript3) != 1,'CCDS transcript lookup with non-coding');

my $canonical_transcript = $transcript_selector->select_canonical_transcript_for_Gene($gene);
is($canonical_transcript->stable_id, 'ENST02','Full select canonical transcript');

my $transcript4 = Bio::EnsEMBL::Transcript->new(
    -DBID => 5,
    -STABLE_ID => 'ENST04',
    -BIOTYPE => 'protein_coding',
    -IS_CURRENT => 1,
    -SLICE => $slice,
    -ANALYSIS => $analysis,
);
$transcript4->add_Exon($exon);

my $transcript5 = Bio::EnsEMBL::Transcript->new(
    -DBID => 4,
    -STABLE_ID => 'ENST05',
    -BIOTYPE => 'protein_coding',
    -IS_CURRENT => 1,
    -SLICE => $slice,
    -ANALYSIS => $analysis,
);
$transcript5->add_Exon($exon);

$transcripts = [ $transcript1, $transcript3, $transcript4, $transcript5 ];

$gene = Bio::EnsEMBL::Gene->new(
    -START  => 123,
    -END    => 2045,
    -STRAND => 1,
    -BIOTYPE => 'protein_coding',
    -TRANSCRIPTS => $transcripts,
    -SLICE => $slice,
);
$canonical_transcript = $transcript_selector->select_canonical_transcript_for_Gene($gene);
is($canonical_transcript->stable_id, 'ENST01', 'Sorting with no CCDS option and equal lengths');


$gene = Bio::EnsEMBL::Gene->new(
    -STABLE_ID => 'ENSGFAKE',
    -START  => 123,
    -END    => 2045,
    -STRAND => 1,
    -BIOTYPE => 'protein_coding',
    -SLICE => $slice,
);
warns_like { $canonical_transcript = $transcript_selector->select_canonical_transcript_for_Gene($gene)} qr/No transcripts attached/, 'Checking we warn about lack of transcripts against a gene';
note ($canonical_transcript);
is($canonical_transcript, undef, "Gene with no transcripts, fault tolerance");

done_testing();
