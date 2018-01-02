-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2018] EMBL-European Bioinformatics Institute
-- 
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
-- 
--      http://www.apache.org/licenses/LICENSE-2.0
-- 
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

# patch_55_56_b.sql
#
# Title: Replace unnamed indexes with named ones.
#
# Description:
# Some indexes are not explicitly named.  Drop these and recreate them
# with non-default names.  Unbreaks health-checking.

-- The 'analysis' table has for seom reason two indexes on 'logic_name',
-- one unique and one not unique.  Drop both and recreate the unique
-- index.
ALTER TABLE analysis
  DROP INDEX logic_name,
  DROP INDEX logic_name_idx,
  ADD UNIQUE INDEX logic_name_idx (logic_name);

-- The 'translation' table has an unnamed index for its 'transcript_id'
-- field.
ALTER TABLE translation
  DROP INDEX transcript_id,
  ADD INDEX transcript_idx (transcript_id);

-- The 'assembly' table has two unnamed indexes.
ALTER TABLE assembly
  DROP INDEX cmp_seq_region_id,
  DROP INDEX asm_seq_region_id,
  ADD INDEX cmp_seq_region_idx (cmp_seq_region_id),
  ADD INDEX asm_seq_region_idx (asm_seq_region_id, asm_start);

-- The 'protein_feature' table...
ALTER TABLE protein_feature
  DROP INDEX translation_id,
  ADD INDEX translation_idx (translation_id);

-- The 'interpro' table...
ALTER TABLE interpro
  DROP INDEX interpro_ac,
  DROP INDEX id,
  ADD UNIQUE INDEX accession_idx (interpro_ac, id),
  ADD INDEX id_idx (id);

-- 'object_xref'...
ALTER TABLE object_xref
  DROP INDEX ensembl_object_type,
  ADD UNIQUE INDEX object_type_idx (ensembl_object_type, ensembl_id, xref_id);

-- 'go_xref'...
ALTER TABLE go_xref
  DROP INDEX source_xref_id,
  DROP INDEX object_xref_id,
  ADD INDEX source_idx (source_xref_id),
  ADD UNIQUE INDEX object_source_type_idx (object_xref_id, source_xref_id, linkage_type);

-- 'prediction_exon'...
ALTER TABLE prediction_exon
  DROP INDEX prediction_transcript_id,
  DROP INDEX seq_region_id,
  ADD INDEX transcript_idx (prediction_transcript_id),
  ADD INDEX seq_region_idx (seq_region_id, seq_region_start);

-- 'prediction_transcript'...
ALTER TABLE prediction_transcript
  DROP INDEX seq_region_id,
  ADD INDEX seq_region_idx (seq_region_id, seq_region_start);

-- The tables 'attrib_type' and 'misc_set' has an index (each) called
-- 'c', drop this and recreate with a more useful name.
ALTER TABLE attrib_type
  DROP INDEX c,
  ADD UNIQUE INDEX code_idx (code);
ALTER TABLE misc_set
  DROP INDEX c,
  ADD UNIQUE INDEX code_idx (code);

-- 'qtl_feature'...
ALTER TABLE qtl_feature
  DROP INDEX qtl_id,
  ADD INDEX qtl_idx (qtl_id);

-- 'density_type'...
ALTER TABLE density_type
  DROP INDEX analysis_id,
  ADD UNIQUE INDEX analysis_idx (analysis_id, block_size, region_features);

-- 'ditag_feature'...
ALTER TABLE ditag_feature
  DROP INDEX ditag_id,
  DROP INDEX ditag_pair_id,
  ADD INDEX ditag_idx (ditag_id),
  ADD INDEX ditag_pair_idx (ditag_pair_id);

-- 'unconventional_transcript_association'...
ALTER TABLE unconventional_transcript_association
  DROP INDEX transcript_id,
  DROP INDEX gene_id,
  ADD INDEX transcript_idx (transcript_id),
  ADD INDEX gene_idx (gene_id);

-- 'seq_region_mapping'...
ALTER TABLE seq_region_mapping
  DROP INDEX mapping_set_id,
  ADD INDEX mapping_set_idx (mapping_set_id);

-- Optimize affected tables.
OPTIMIZE TABLE analysis;
OPTIMIZE TABLE translation;
OPTIMIZE TABLE assembly;
OPTIMIZE TABLE protein_feature;
OPTIMIZE TABLE interpro;
OPTIMIZE TABLE object_xref;
OPTIMIZE TABLE go_xref;
OPTIMIZE TABLE prediction_exon;
OPTIMIZE TABLE prediction_transcript;
OPTIMIZE TABLE attrib_type;
OPTIMIZE TABLE misc_set;
OPTIMIZE TABLE qtl_feature;
OPTIMIZE TABLE density_type;
OPTIMIZE TABLE ditag_feature;
OPTIMIZE TABLE unconventional_transcript_association;
OPTIMIZE TABLE seq_region_mapping;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
VALUES (NULL, 'patch', 'patch_55_56_b.sql|add_index_names');
