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

use Test::Warnings;
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

my $rcoding_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['coding_rcnt'], -NO_ERROR => 1);
is($rcoding_count, $genome->get_rcoding_count, "Readthough coding count is correct");

my $lnoncoding_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['noncoding_cnt_l']);
is($lnoncoding_count, $genome->get_lnoncoding_count, "Long non coding count is correct");

my $rlnoncoding_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['noncoding_rcnt_l'], -NO_ERROR => 1);
is($rlnoncoding_count, $genome->get_rlnoncoding_count, "Readthrough long non coding count is correct");

my $snoncoding_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['noncoding_cnt_s']);
is($snoncoding_count, $genome->get_snoncoding_count, "Short non coding count is correct");

my $rsnoncoding_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['noncoding_rcnt_s'], -NO_ERROR => 1);
is($rsnoncoding_count, $genome->get_rsnoncoding_count, "Readthrough short non coding count is correct");

my $mnoncoding_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['noncoding_cnt_m'], -NO_ERROR => 1);
is($mnoncoding_count, $genome->get_mnoncoding_count, "Miscellaneous non coding count is correct");

my $rmnoncoding_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['noncoding_rcnt_m'], -NO_ERROR => 1);
is($rmnoncoding_count, $genome->get_rmnoncoding_count, "Readthrough miscellaneous non coding count is correct");

my $pseudogene_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['pseudogene_cnt']);
is($pseudogene_count, $genome->get_pseudogene_count, "Pseudogene count is correct");

my $rpseudogene_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['pseudogene_rcnt'], -NO_ERROR => 1);
is($rpseudogene_count, $genome->get_rpseudogene_count, "Readthrough pseudogene count is correct");

my $alt_rpseudogene_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['pseudogene_racnt'], -NO_ERROR => 1);
is($alt_rpseudogene_count, $genome->get_alt_rpseudogene_count, "Readthrough pseudogene count on alternate sequences is correct");

my $alt_coding_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['coding_acnt']);
is($alt_coding_count, $genome->get_alt_coding_count, "Coding count on alternate sequences is correct");

my $alt_rcoding_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['coding_racnt'], -NO_ERROR => 1);
is($alt_rcoding_count, $genome->get_alt_rcoding_count, "Readthrough coding count on alternate sequences is correct");

my $alt_lnoncoding_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['lnoncoding_acnt']);
is($alt_lnoncoding_count, $genome->get_alt_lnoncoding_count, "Long non coding count on alternate sequences is correct");

my $alt_rlnoncoding_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['lnoncoding_racnt'], -NO_ERROR => 1);
is($alt_rlnoncoding_count, $genome->get_alt_rlnoncoding_count, "Readthrough long non coding count on alternate sequences is correct");

my $alt_snoncoding_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['snoncoding_acnt']);
is($alt_snoncoding_count, $genome->get_alt_snoncoding_count, "Short non coding count on alternate sequences is correct");

my $alt_rsnoncoding_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['snoncoding_racnt'], -NO_ERROR => 1);
is($alt_rsnoncoding_count, $genome->get_alt_rsnoncoding_count, "Readthrough short non coding count on alternate sequences is correct");

my $alt_mnoncoding_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['mnoncoding_acnt'], -NO_ERROR => 1);
is($alt_mnoncoding_count, $genome->get_alt_mnoncoding_count, "Short non coding count on alternate sequences is correct");

my $alt_rmnoncoding_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['mnoncoding_racnt'], -NO_ERROR => 1);
is($alt_rmnoncoding_count, $genome->get_alt_rmnoncoding_count, "Readthrough short non coding count on alternate sequences is correct");

my $short_variation_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['SNPCount'], -NO_ERROR => 1);
is($short_variation_count, $genome->get_short_variation_count, "Short variants count is correct");

my $structural_variation_count = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => ['StructuralVariation'], -NO_ERROR => 1);
is($structural_variation_count, $genome->get_structural_variation_count, "Structural variants count is correct");

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
is($total_length, $genome->get_total_length, "Total length is correct");

#
# Test transcript counts
#

my $transcript_sql = "select count(*) from transcript t, seq_region s
where t.seq_region_id = s.seq_region_id
and t.seq_region_id not in (
select sa.seq_region_id from seq_region_attrib sa, attrib_type at
where at.attrib_type_id = sa.attrib_type_id
and at.code = 'non_ref')";
my $transcript_count = $sql_helper->execute_single_result(-SQL => $transcript_sql);
is($transcript_count, $genome->get_transcript_count(), "Number of transcripts is correct");
my $alt_transcript_sql = "select count(*) from transcript t, seq_region s, seq_region_attrib sa, attrib_type at
where s.seq_region_id = t.seq_region_id
and s.seq_region_id = sa.seq_region_id
and sa.attrib_type_id = at.attrib_type_id
and at.code = 'non_ref'
and biotype not in ('LRG_gene')";
my $alt_transcript_count = $sql_helper->execute_single_result(-SQL => $alt_transcript_sql);
is($alt_transcript_count, $genome->get_alt_transcript_count(), "Number of alt transcripts is correct");

my $counts = $genome->fetch_all_statistics();
is(scalar(@$counts), 16, "All separate statistics retrieved");
foreach my $genome_object (@$counts) {
  debug($genome_object->statistic);
  debug($genome_object->code);
  debug($genome_object->name);
  debug($genome_object->description);
  debug($genome_object->dbID);
  debug($genome_object->species_id);
}

#
# Test region methods
#

my $toplevel = $genome->get_toplevel();
is(scalar(@$toplevel), 2, "2 toplevel regions");
my $karyotype = $genome->get_karyotype();
is(scalar(@$karyotype), 2, "2 karyotype regions");
my $coord_systems = $genome->get_coord_systems();
is(scalar(@$coord_systems), 2, "2 coord systems");

#
# Test genebuild methods
#

my $gb_start = $genome->get_genebuild_start_date();
is($gb_start, undef, 'Genebuild start date');
my $gb_method = $genome->get_genebuild_method();
is($gb_method, undef, 'Genebuild method');
my $gb_release = $genome->get_genebuild_initial_release_date();
is($gb_release, undef, 'Release date');
my $gb_update = $genome->get_genebuild_last_geneset_update();
is($gb_update, undef, 'Update date');


#
# Test karyotype flag
#

my $empty_db = $multi->get_DBAdaptor("empty");
my $empty_genome = $empty_db->get_adaptor('GenomeContainer');

is(1, $genome->has_karyotype, "Human has some chromosomes");
is(0, $empty_genome->has_karyotype, "Empty db does not have chromosomes");


#
# Test high coverage flag
#

is(1, $genome->is_high_coverage, "Human genebuild is high coverage");
is(0, $empty_genome->is_high_coverage, "Empty db is not high coverage");


#
# Test if there is data
#

is(0, $genome->is_empty, "Human core database has genome data");
is(1, $empty_genome->is_empty, "Empty database has no genome data");

#
# Test polyploid genome support
#
# get a genome container for a non polyploid genome core db (human)
my $human = $multi->get_DBAdaptor("core");
my $hgdba = $human->get_adaptor('GenomeContainer');

ok($hgdba && $hgdba->isa('Bio::EnsEMBL::DBSQL::GenomeContainer'), 'GenomeContainer adaptor');
ok(!$hgdba->is_polyploid, "Human genome is not polyploid");
is_deeply($hgdba->get_genome_components(), [], "Human does not have genome components");

# get a genome container for a polyploid genome core db (bread wheat)
my $multi_polyploid = Bio::EnsEMBL::Test::MultiTestDB->new("polyploidy");
my $wheat = $multi_polyploid->get_DBAdaptor("core");
my $wgdba = $wheat->get_adaptor('GenomeContainer');

ok($wgdba && $wgdba->isa('Bio::EnsEMBL::DBSQL::GenomeContainer'), 'GenomeContainer adaptor');
ok($wgdba->is_polyploid, "Triticum aestivum genome is polyploid");
is_deeply($wgdba->get_genome_components(), ['A','B','D'], "Triticum aestivum genome components");

done_testing();
