# patch_61_62_d.sql
#
# Title: Remove field display_label_linkable from external_db table.
#
# Description:
# Field is obsolete.

ALTER TABLE external_db DROP display_label_linkable;

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_61_62_d.sql|remove_display_label_linkable');
