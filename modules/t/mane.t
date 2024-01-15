# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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
use Test::Warnings;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

require_ok('Bio::EnsEMBL::MANE');

our $verbose = 0; #set to 1 to turn on debug printouts

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

ok(1);

my $db = $multi->get_DBAdaptor( 'core' );

ok($db);

my $stable_id = 'ENST00000202017';
my $transcript_adaptor = $db->get_TranscriptAdaptor();
my $transcript = 
  $transcript_adaptor->fetch_by_stable_id($stable_id);

print "Check if transcript is part of MANE\n";
is($transcript->is_mane, 1, "$stable_id is part of mane");

my $mane_transcript = $transcript->mane_transcript();
print "Retrieve MANE attributes\n";
is($mane_transcript->stable_id, $transcript->stable_id, "MANE transcript has same stable ID as the original transcript");
is($mane_transcript->refseq, "NM_030815.3", "MANE transcript matches RefSeq accession");
is($mane_transcript->type, "MANE_Select", "MANE transcript belongs to MANE_Select");

my $stable_id2 = "ENST00000252021";
my $transcript2 = $transcript_adaptor->fetch_by_stable_id($stable_id2);

print "Check transcript is not part of MANE\n";
is($transcript2->is_mane, 0, "$stable_id2 is not part of MANE");
my $mane_transcript2 = $transcript2->mane_transcript();
is($mane_transcript2, undef, "No result when the transcript is not part of MANE");

done_testing();
