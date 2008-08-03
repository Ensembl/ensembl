# patch_50_51_i.sql
#
# title: Meta value binary
#
# description:
# Add BINARY flag to meta.meta_value to make it case sensitive when doing UNIQUE comparisions (especially for Gallus/gallus)

# Drop UNIQUE index first
ALTER TABLE meta DROP INDEX species_key_value_idx;

ALTER TABLE meta CHANGE COLUMN meta_value meta_value VARCHAR(255) BINARY NOT NULL;

# Redo index
ALTER TABLE meta ADD UNIQUE INDEX species_key_value_idx (species_id, meta_key, meta_value);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_50_51_i.sql|meta_value_binary');


