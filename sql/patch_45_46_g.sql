# patch_45_46_g.sql
#
# title: object_xref linkage_annotation
#
# description:
# new column added to object_xref to keep info about teh object_xref

ALTER TABLE object_xref ADD COLUMN linkage_annotation VARCHAR(255) DEFAULT NULL;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_45_46_g.sql|object_xref_linkage_annotation');



