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

# patch_43_44_a
#
# title: Key column types
#
# description:
# Change the types of some columns to make them more appropriate to the values they store.

ALTER TABLE analysis CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL AUTO_INCREMENT;
ALTER TABLE oligo_feature CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE analysis_description CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE simple_feature CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE protein_align_feature CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE dna_align_feature CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE repeat_feature CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE gene CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE transcript CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE protein_feature CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE identity_xref CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE prediction_transcript CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE marker_feature CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE qtl_feature CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE density_type CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE regulatory_feature CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE regulatory_search_region CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE unmapped_object CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE ditag_feature CHANGE COLUMN analysis_id analysis_id SMALLINT UNSIGNED NOT NULL;

ALTER TABLE external_db CHANGE COLUMN external_db_id external_db_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE xref CHANGE COLUMN external_db_id external_db_id SMALLINT UNSIGNED NOT NULL;
ALTER TABLE unmapped_object CHANGE COLUMN external_db_id external_db_id SMALLINT UNSIGNED NOT NULL;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_43_44_a.sql|rationalise_key_columns');

