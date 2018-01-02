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

# patch_67_68_b.sql
#
# Title: Xref unique constraint enforcement
#
# Description:
#   Remove null values from xref and object_xref tables enforcing the unique
#   index. Fields are unique on:
#
#   XREF - primary accession, version, DB, info type and info text
#   OBJECT XREF - xref id, ensembl id, ensembl object type and analysis id
#
#   This means we now force info type to be NONE, info text to be '' and 
#   analysis id to be 0 IF THEY WOULD HAVE BEEN NULL. 
#
#   See also DBEntryAdaptor thread safety changes

## NB A lot of the functions here use IFNULL(). We could have used the NULL safe comparison operator <=> but were unaware at the time of its existence

# Remove duplicate nulls in xref table

# Need to find the duplicates first and select the lowest xref_id as our "canonical" xref_id

create temporary table xref_dups
select `dbprimary_acc`,`version`,`external_db_id`,IFNULL(`info_type`, 'NONE') as info_type, IFNULL(`info_text`, '') as info_text, min(xref_id) as xref_id, count(*) as c
from xref
group by `dbprimary_acc`,`version`,`external_db_id`,IFNULL(`info_type`, 'NONE'),IFNULL(`info_text`, '')
having c > 1;

# Mark all other duplicate xrefs and flag their new canonical ID

create temporary table xref_MFD
select x.xref_id, xd.xref_id AS canonical_xref_id
from xref x join `xref_dups` xd on ( 
  x.`dbprimary_acc` = xd.`dbprimary_acc`
  and x.`version` = xd.`version`
  and x.`external_db_id` = xd.`external_db_id` 
  and IFNULL(x.`info_type`, 'NONE') = xd.`info_type`
  and IFNULL(x.`info_text`, '') = xd.info_text
  and xd.`xref_id` <> x.`xref_id`
);

# Remove the unique constraint
ALTER TABLE xref DROP KEY id_index;

# Delete the duplicates

DELETE FROM xref USING xref JOIN xref_MFD WHERE xref.xref_id = xref_MFD.xref_id;

# Update object_xref, dependent_xref, ontology_xref xref_ids to the canonical ones

UPDATE object_xref ox join xref_MFD xmfd using (xref_id) set ox.xref_id = xmfd.canonical_xref_id;

UPDATE dependent_xref dx join xref_MFD xmfd on (dx.master_xref_id = xmfd.xref_id) set dx.master_xref_id = xmfd.canonical_xref_id;
UPDATE dependent_xref dx join xref_MFD xmfd on (dx.dependent_xref_id = xmfd.xref_id) set dx.dependent_xref_id = xmfd.canonical_xref_id;

UPDATE ontology_xref ox join xref_MFD xmfd on (ox.source_xref_id = xmfd.xref_id) set ox.source_xref_id = xmfd.canonical_xref_id;

UPDATE gene g join xref_MFD xmfd on (g.display_xref_id = xmfd.xref_id) set g.display_xref_id = xmfd.canonical_xref_id;
UPDATE transcript t join xref_MFD xmfd on (t.display_xref_id = xmfd.xref_id) set t.display_xref_id = xmfd.canonical_xref_id;


# Apply info_text update and set not null constraint
UPDATE xref SET info_text='' WHERE info_text is NULL;

ALTER TABLE xref MODIFY info_text varchar(255) DEFAULT '' NOT NULL;

# Apply none type, update NULL to NONE and then apply not null constraint

ALTER TABLE xref MODIFY info_type enum('NONE','PROJECTION','MISC','DEPENDENT','DIRECT','SEQUENCE_MATCH','INFERRED_PAIR','PROBE','UNMAPPED','COORDINATE_OVERLAP','CHECKSUM') DEFAULT 'NONE';

UPDATE xref SET info_type='NONE' WHERE info_type is NULL;

ALTER TABLE xref MODIFY info_type enum('NONE','PROJECTION','MISC','DEPENDENT','DIRECT','SEQUENCE_MATCH','INFERRED_PAIR','PROBE','UNMAPPED','COORDINATE_OVERLAP','CHECKSUM') DEFAULT 'NONE' NOT NULL;

# Add the constraint back
ALTER TABLE xref ADD UNIQUE KEY id_index (dbprimary_acc, external_db_id, info_type, info_text, version);  

# Remove duplicate nulls in object_xref table

create temporary table object_xref_dups
select `ensembl_id`, `ensembl_object_type`, `xref_id`, IFNULL(`analysis_id`, 0) as analysis_id, min(object_xref_id) as object_xref_id, count(*) as c
from object_xref
group by `ensembl_id`, `ensembl_object_type`, `xref_id`, IFNULL(`analysis_id`, 0)
having c > 1;

create temporary table object_xref_MFD
select ox.object_xref_id
from object_xref ox join `object_xref_dups` oxd on ( 
  ox.`ensembl_id` = oxd.`ensembl_id` 
  and ox.ensembl_object_type = oxd.ensembl_object_type 
  and ox.xref_id = oxd.xref_id
  and IFNULL(ox.analysis_id, 0) = oxd.analysis_id
  and oxd.`object_xref_id` <> ox.`object_xref_id`
);
ALTER TABLE object_xref_MFD ADD INDEX dribbling_simpleton(object_xref_id);

DELETE FROM object_xref USING object_xref JOIN object_xref_MFD WHERE object_xref.object_xref_id = object_xref_MFD.object_xref_id;

UPDATE object_xref SET analysis_id = 0 WHERE analysis_id is NULL;
ALTER TABLE object_xref MODIFY analysis_id smallint(5) unsigned DEFAULT 0 NOT NULL;

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_67_68_b.sql|xref_uniqueness');
