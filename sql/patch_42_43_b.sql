# patch_42_43_b
#
# title: Add probe2transcript to type enum in unmapped_object
#
# description:
# Add probe2transcript to type enum in unmapped_object

ALTER TABLE unmapped_object CHANGE COLUMN type type ENUM('xref','cDNA','Marker', 'probe2transcript') NOT NULL;


# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_42_43_b.sql|unmapped_object_probe2transcript');