# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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
use Test::Deep;
use Test::Warnings qw(warning);
use Test::Exception;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::BaseAlignFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::ProteinFeature;
use Bio::EnsEMBL::Analysis;


my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');
my $pfa = $db->get_ProteinFeatureAdaptor();

my $analysis_db = 'test_db';

my $analysis = Bio::EnsEMBL::Analysis->new(-LOGIC_NAME => 'gifts_import', -DB => $analysis_db);

my $seq_start  = 1;
my $seq_end    = 111;
my $hit_start = 0;
my $hit_end   = 0;
my $hit_name = "P20366";
my $cigar_string = "MD:Z:0G105";
my $align_type = "mdtag";
my $translation_id = 21739;

my $pf = Bio::EnsEMBL::ProteinFeature->new
  (-START       => $seq_start,
   -END         => $seq_end,
   -ANALYSIS    => $analysis,
   -HSTART      => $hit_start,
   -HEND        => $hit_end,
   -HSEQNAME    => $hit_name,
   -CIGAR_STRING => $cigar_string,
   -ALIGN_TYPE => $align_type,
   -TRANSLATION_ID     => $translation_id,
   -ADAPTOR           => $pfa
   );


ok($pf && $pf->isa('Bio::EnsEMBL::ProteinFeature'));

my $alignment_string = $pf->alignment_strings();
is($alignment_string->[0], "DNSSLSGEERLKCKLGKSFLLEKSLGKGMLIHCSLGVSMGKGKPPSPLTLTSFPPFCDLAKSAFHVVLTTTGVKLTMIPYSRSRLMSSEDLAEIPQLQKLSIPHGF", "Got query string right");
is($alignment_string->[1], "GNSSLSGEERLKCKLGKSFLLEKSLGKGMLIHCSLGVSMGKGKPPSPLTLTSFPPFCDLAKSAFHVVLTTTGVKLTMIPYSRSRLMSSEDLAEIPQLQKLSIPHGF", "Got target string right");


my $chunks = $pf->_get_mdz_chunks("MD:Z:0G105");
my $expected = ['0', 'G', '105'];
cmp_deeply($chunks, $expected, "Got the right chunks ['0', 'G', '105']");

$chunks = $pf->_get_mdz_chunks("MD:Z:96^RHKTDSFVGLMGKRALNS0V14");
$expected = ['96', '^', 'RHKTDSFVGLMGKRALNS', '0', 'V', '14'];
cmp_deeply($chunks, $expected, "Got the right chunks ['96', '^', 'RHKTDSFVGLMGKRALNS', '0', 'V', '14']");

$chunks = $pf->_get_mdz_chunks("MD:Z:35^VIVALE31^GRPLIQPRRKKAYQLEHTFQGLLGKRSLFTE10");
$expected = ['35', '^', 'VIVALE', '31', '^', 'GRPLIQPRRKKAYQLEHTFQGLLGKRSLFTE', '10'];
cmp_deeply($chunks, $expected, "Got the right chunks  ['35', '^', 'VIVALE', '31', '^', 'GRPLIQPRRKKAYQLEHTFQGLLGKRSLFTE', '10']");



done_testing();
