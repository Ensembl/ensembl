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
use vars qw( $verbose );

use Test::More;
use Test::Warnings;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

$verbose = 0;

ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

ok($multi);

my $db = $multi->get_DBAdaptor( "core" );

my $gene = $db->get_GeneAdaptor->fetch_by_transcript_stable_id( "ENST00000217347" );
my $geneid = $gene->stable_id;
ok( $geneid );

$gene = $db->get_GeneAdaptor->fetch_by_translation_stable_id( "ENSP00000278995" );
$geneid = $gene->stable_id;

$gene = $db->get_GeneAdaptor()->fetch_by_stable_id( "ENSG00000101321" );
ok( $gene );


done_testing();
