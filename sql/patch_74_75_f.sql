# patch_74_75_f.sql
#
# title: Longer code column in attrib_type table
#
# description:
# Code column in attrib_type needs to be longer

ALTER TABLE attrib_type MODIFY COLUMN code VARCHAR(20) NOT NULL DEFAULT '';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_74_75_f.sql|longer_code');

 
