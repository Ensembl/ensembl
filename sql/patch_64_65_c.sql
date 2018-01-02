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

# patch_62_63_c.sql
#
# Title: New table for file location storage
#
# Description:
# We now record the location of an external file to be used by the website
# for display. This is primarily for BAM data but we will extend this
# for other formats.

CREATE TABLE data_file (
	data_file_id int(11) unsigned NOT NULL AUTO_INCREMENT,
	coord_system_id int(11) NOT NULL,
	analysis_id int(11) NOT NULL,
	name varchar(100) NOT NULL,
	version_lock tinyint(1) DEFAULT 0 NOT NULL,
	absolute tinyint(1) DEFAULT 0 NOT NULL,
	url text,
	file_type enum('BAM','BIGBED','BIGWIG','VCF'),
	PRIMARY KEY (data_file_id),
  UNIQUE KEY df_unq_idx(coord_system_id, analysis_id, name, file_type),
  INDEX df_name_idx(name),
  INDEX df_analysis_idx(analysis_id)
) ENGINE=MyISAM;

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_64_65_c.sql|add_data_file');
