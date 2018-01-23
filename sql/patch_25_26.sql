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

# 3 new tables to deal with affy probes

CREATE TABLE affy_feature (
       affy_feature_id INT NOT NULL auto_increment,
       seq_region_id INT UNSIGNED NOT NULL,
       seq_region_start INT NOT NULL,
       seq_region_end INT NOT NULL,
       seq_region_strand TINYINT NOT NULL,
       
       mismatches TINYINT,
       affy_probe_id INT NOT NULL,
       analysis_id INT NOT NULL,

       PRIMARY KEY (affy_feature_id),
       KEY seq_region_idx( seq_region_id, seq_region_start ),
       KEY probe_idx( affy_probe_id )
);


CREATE TABLE affy_probe (
       affy_probe_id INT NOT NULL auto_increment,
       affy_array_id INT NOT NULL,
       probeset VARCHAR(20),
       name VARCHAR(20),

       PRIMARY KEY ( affy_probe_id, affy_array_id ),
       KEY probeset_idx( probeset ),
       KEY array_idx( affy_array_id )
);


CREATE TABLE affy_array (
       affy_array_id INT NOT NULL auto_increment,
       parent_array_id INT,
       probe_setsize TINYINT NOT NULL,
       name VARCHAR(40) NOT NULL,

       PRIMARY KEY( affy_array_id )
);

# Shorten some db_name entries in the external_db table 
UPDATE external_db SET db_name='Genoscope_pred_transcript' WHERE db_name='Genoscope_predicted_transcript';
UPDATE external_db SET db_name='Genoscope_pred_gene' WHERE db_name='Genoscope_predicted_gene';
ALTER TABLE external_db MODIFY db_name VARCHAR(27);

# Adds an index no the seq_region_id column (only) of density_feature
# This speeds up some queries very significantly, particularly those 
# used in some healthchecks.

ALTER TABLE density_feature ADD INDEX seq_region_id_idx (seq_region_id);
