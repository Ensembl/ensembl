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

# patch_39_40_g
#
# title: add ditag tables
#
# description:
# This patch adds two new tables, ditag and ditag_feature

CREATE TABLE ditag (

       ditag_id INT NOT NULL auto_increment,
       name VARCHAR(30),
       type VARCHAR(30),
       tag_count smallint(6) default 1,
       sequence TEXT,

       PRIMARY KEY ( ditag_id )

) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE ditag_feature (

       ditag_feature_id int(10) unsigned NOT NULL auto_increment,
       ditag_id int(10) unsigned NOT NULL default '0',
       ditag_pair_id int(10) unsigned NOT NULL default '0',
       seq_region_id int(10) unsigned NOT NULL default '0',
       seq_region_start int(10) unsigned NOT NULL default '0',
       seq_region_end int(10) unsigned NOT NULL default '0',
       seq_region_strand tinyint(1) NOT NULL default '0',
       analysis_id int(10) unsigned NOT NULL default '0',
       hit_start int(10) unsigned NOT NULL default '0',
       hit_end int(10) unsigned NOT NULL default '0',
       hit_strand tinyint(1) NOT NULL default '0',
       cigar_line text default '',
       ditag_side char default '',

       PRIMARY KEY  (ditag_feature_id),
       KEY (ditag_id),
       KEY (ditag_pair_id)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_39_40_g.sql|add_ditag_tables');

