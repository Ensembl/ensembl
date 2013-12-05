# Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

use Bio::EnsEMBL::DBSQL::GenomeContainer;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Test::MultiTestDB;


our $verbose = 0;

use Test::More;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();


#
# Test constructor
#

my $db = $multi->get_DBAdaptor("patch");
debug("Test database instatiated");
ok($db);
my $sql_helper = $db->dbc->sql_helper;
my $genome = $db->get_adaptor('GenomeContainer');

ok($genome && $genome->isa('Bio::EnsEMBL::DBSQL::GenomeContainer'));

# 
# Test version()
#

my $sql = "SELECT DISTINCT version FROM coord_system
        WHERE attrib like '%default_version%'
          AND version is not null";
my $version = $sql_helper->execute_single_result(-SQL => $sql);
is($genome->get_version(), $version, "Genome version is correct");

$sql = "SELECT meta_value FROM meta
         WHERE meta_key = 'assembly.accession'";
my $accession = $sql_helper->execute_single_result(-SQL => $sql);
is($genome->get_accession(), $accession, "Genome accession is correct");


#
# Test feature counts
#

$sql = "select sum(value) from seq_region_attrib sa, attrib_type at, seq_region s
        where at.attrib_type_id = sa.attrib_type_id
        and s.seq_region_id = sa.seq_region_id
        and code = ?";
my $coding_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['coding_cnt']);
is($coding_count, $genome->get_coding_count, "Coding count is correct");

my $lnoncoding_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['lnoncoding_cnt']);
is($lnoncoding_count, $genome->get_lnoncoding_count, "Long non coding count is correct");

my $snoncoding_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['snoncoding_cnt']);
is($snoncoding_count, $genome->get_snoncoding_count, "Short non coding count is correct");

my $pseudogene_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['pseudogene_cnt']);
is($pseudogene_count, $genome->get_pseudogene_count, "Pseudogene count is correct");

my $alt_coding_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['coding_acnt']);
is($alt_coding_count, $genome->get_alt_coding_count, "Coding count on alternate sequences is correct");

my $alt_lnoncoding_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['lnoncoding_acnt']);
is($alt_lnoncoding_count, $genome->get_alt_lnoncoding_count, "Long non coding count on alternate sequences is correct");

my $alt_snoncoding_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['snoncoding_acnt']);
is($alt_snoncoding_count, $genome->get_alt_snoncoding_count, "Short non coding count on alternate sequences is correct");

my $short_variation_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['SNPCount'], -NO_ERROR => 1);
is($short_variation_count, $genome->get_short_variation_count, "Short variants count is correct");

is_rows($genome->get_prediction_count, $db, "prediction_transcript");
is_rows($genome->get_prediction_count('genscan'), $db, "prediction_transcript", "where analysis_id = ?", [8440]);

#
# Test genome length
#

$sql = "SELECT sum(length) FROM seq_region sr, seq_region_attrib sra, attrib_type at, coord_system cs
        WHERE sr.seq_region_id = sra.seq_region_id
          AND sra.attrib_type_id = at.attrib_type_id
          AND sr.coord_system_id = cs.coord_system_id 
          AND at.code = 'toplevel' 
          AND cs.name != 'lrg' 
          AND sr.seq_region_id NOT IN 
            (SELECT DISTINCT seq_region_id FROM assembly_exception ae WHERE ae.exc_type != 'par' )";
my $ref_length = $sql_helper->execute_single_result(-SQL => $sql); 
is($ref_length, $genome->get_ref_length, "Reference length is correct");

$sql = "SELECT sum(length(sequence)) FROM dna";
my $total_length = $sql_helper->execute_single_result(-SQL => $sql);    
print "$total_length is the calculated sql length\n";
is($total_length, $genome->get_total_length, "Total length is correct");



done_testing();
