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


# patch_49_50_e.sql
#
# title: Create tables for new mapping_seq_region tables
#
# Description: Creates two new tables, seq_region_mapping and mapping_seq that will allow to upload user data and be able to map it
# even there is a change of seq_region in a previous release

################################################################################
#
# Table structure for seq_region mapping between releases
#
# Stores how the core seq_region_id have changed from release to release

CREATE TABLE seq_region_mapping (

	external_seq_region_id	INT(10) UNSIGNED NOT NULL,
	internal_seq_region_id	INT(10) UNSIGNED NOT NULL,
	mapping_set_id		INT(10) UNSIGNED NOT NULL,

	KEY (mapping_set_id)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;

################################################################################
#
# Table structure for seq_region mapping between releases
#
# Stores how which mapping group the seq_region are for a particular schema

CREATE TABLE mapping_set (

	mapping_set_id	INT(10)	UNSIGNED NOT NULL,
	schema_build	VARCHAR(20) NOT NULL,

	PRIMARY KEY(schema_build)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;

# Patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_49_50_e.sql|mapping_seq_region');



