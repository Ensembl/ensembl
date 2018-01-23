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


CREATE TABLE translation_attrib (
  translation_id              int(10) unsigned NOT NULL default '0',
  attrib_type_id              smallint(5) unsigned NOT NULL default '0',
  value                       varchar(255) NOT NULL default '',

  KEY type_val_idx( attrib_type_id, value ),
  KEY translation_idx( translation_id )
) TYPE=MyISAM;


CREATE TABLE transcript_attrib (
  transcript_id               int(10) unsigned NOT NULL default '0',
  attrib_type_id              smallint(5) unsigned NOT NULL default '0',
  value                       varchar(255) NOT NULL default '',

  KEY type_val_idx( attrib_type_id, value ),
  KEY transcript_idx( transcript_id )
) TYPE=MyISAM;

alter table misc_set modify code varchar(25) NOT NULL default ''

