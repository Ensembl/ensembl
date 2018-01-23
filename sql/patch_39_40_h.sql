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

# patch_39_40_h
#
# title: oligo_feature analysis id column type change
#
# description:
# Change oligo_feature.analysis_id to be int(10) unsigned. Should have been part of patch_39_40_a.sql

ALTER TABLE oligo_feature CHANGE COLUMN analysis_id analysis_id INT(10) UNSIGNED NOT NULL;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_39_40_h.sql|oligo_feature_analysis_id_type');

