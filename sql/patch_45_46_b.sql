# patch_45_46_b.sql
#
# title: go_xref.source_xref_id
#
# description:
#   Addition of a source_xref_id field to the go_xref table

ALTER TABLE go_xref ADD COLUMN source_xref_id int(10) unsigned default NULL;

ALTER TABLE go_xref DROP KEY object_xref_id_2;

ALTER TABLE go_xref ADD UNIQUE KEY object_xref_id_2 
  ( object_xref_id, linkage_type, source_xref_id );

ALTER TABLE go_xref ADD KEY source_xref_id( source_xref_id );


# patch identifier
INSERT INTO meta (meta_key, meta_value) 
  VALUES ('patch', 'patch_45_46_b.sql|go_xref.source_xref_id');



