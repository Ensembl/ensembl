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

# patch_40_41_g
#
# title: Genebuild version format change
#
# description: Change format of genebuild.version entries in meta table.

UPDATE meta set meta_value = concat('20',substring(meta_value,1,2),'-', substring(meta_value,3,2),'-',substring(meta_value,5)) where meta_key = 'genebuild.version' and meta_value not rlike '^[0-9][0-9][0-9][0-9]-[0-9][0-9]-';

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_41_42_g.sql|genebuild_version_format_change');
