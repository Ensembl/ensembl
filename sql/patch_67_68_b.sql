# patch_67_68_b.sql
#
# Title: 
#
# Description:
#   Remove null values from xref and object_xref tables. See also DBEntryAdaptor thread safety changes


ALTER TABLE xref MODIFY info_type enum('NONE','PROJECTION','MISC','DEPENDENT','DIRECT','SEQUENCE_MATCH','INFERRED_PAIR','PROBE','UNMAPPED','COORDINATE_OVERLAP','CHECKSUM') DEFAULT 'NONE';

UPDATE xref SET info_type='NONE' WHERE info_type is NULL;

ALTER TABLE xref MODIFY info_type enum('NONE','PROJECTION','MISC','DEPENDENT','DIRECT','SEQUENCE_MATCH','INFERRED_PAIR','PROBE','UNMAPPED','COORDINATE_OVERLAP','CHECKSUM') DEFAULT 'NONE' NOT NULL;

UPDATE xref SET info_text='' WHERE info_text is NULL;
ALTER TABLE xref MODIFY info_text varchar(255) DEFAULT '' NOT NULL;

# Remove duplicate nulls in object_xref table

create temporary table object_xref_dups
select `ensembl_id`, `ensembl_object_type`, `xref_id`, `linkage_annotation`, `analysis_id`, min(object_xref_id) as object_xref_id, count(*) as c
from object_xref
group by `ensembl_id`, `ensembl_object_type`, `xref_id`, `linkage_annotation`, `analysis_id`
having c > 1;

create temporary table object_xref_MFD
select ox.object_xref_id
from object_xref ox join `object_xref_dups` oxd on ( 
	ox.`ensembl_id` = oxd.`ensembl_id` 
	and ox.ensembl_object_type = oxd.ensembl_object_type 
	and ox.xref_id = oxd.xref_id
	and (ox.linkage_annotation = oxd.linkage_annotation || (ox.`linkage_annotation` IS NULL and oxd.`linkage_annotation` IS NULL))
	and (ox.analysis_id = oxd.analysis_id || (ox.`analysis_id` IS NULL and oxd.`analysis_id` IS NULL))
	and oxd.`object_xref_id` <> ox.`object_xref_id`
);
ALTER TABLE object_xref_MFD ADD INDEX dribbling_simpleton(object_xref_id);

-- DELETE FROM object_xref WHERE object_xref_id = ANY (SELECT object_xref_id FROM object_xref_MFD);
DELETE FROM object_xref USING object_xref JOIN object_xref_MFD WHERE object_xref.object_xref_id = object_xref_MFD.object_xref_id;

UPDATE object_xref SET analysis_id = 0 WHERE analysis_id is NULL;
ALTER TABLE object_xref MODIFY analysis_id smallint(5) unsigned DEFAULT 0 NOT NULL;

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_67_68_b.sql|xref_uniqueness');
