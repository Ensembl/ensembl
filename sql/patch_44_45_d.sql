# patch_44_45_d.sql
#
# title: go_xref_source_field
#
# description: 
# Add a source_external_db_id field to go_xref table

ALTER TABLE `go_xref` ADD COLUMN 
  `source_external_db_id` int(10) unsigned default NULL;

ALTER TABLE `go_xref` ADD KEY (source_external_db_id);

# patch identifier
INSERT INTO meta (meta_key, meta_value) 
VALUES ('patch', 'patch_44_45_d.sql|go_xref_source_field');
