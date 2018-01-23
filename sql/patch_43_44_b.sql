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

# patch_43_44_b
#
# title: optimising ditag tables
#
# description:
# some small memory & speed improvements for the ditag tables

ALTER TABLE ditag CHANGE COLUMN type type varchar(30) NOT NULL, \
 CHANGE COLUMN tag_count tag_count smallint(6) unsigned NOT NULL default 1, \
 CHANGE COLUMN sequence sequence TINYTEXT NOT NULL, \
 CHANGE COLUMN name name varchar(30) NOT NULL;

ALTER TABLE ditag_feature CHANGE COLUMN cigar_line cigar_line TINYTEXT NOT NULL, \
 CHANGE COLUMN ditag_side ditag_side ENUM('F', 'L', 'R') NOT NULL;
CREATE INDEX seq_region_idx ON ditag_feature (seq_region_id, seq_region_start, seq_region_end);

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_43_44_b.sql|optimise_ditag_tables');

