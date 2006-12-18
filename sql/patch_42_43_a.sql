# patch_42_43_a
#
# title: unmapped_object.parent
#
# description:
# Add parent column to unmapped_object

ALTER TABLE unmapped_object ADD COLUMN parent VARCHAR(255) DEFAULT NULL;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_42_43_a.sql|unmapped_object.parent');


