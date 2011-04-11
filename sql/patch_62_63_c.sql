# patch_62_63_c.sql
#
# Title: Remove field dbprimary_acc_linkable from external_db table.
#
# Description:
# Field is obsolete.

ALTER TABLE external_db DROP dbprimary_acc_linkable;

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_62_63_c.sql|remove_dbprimary_acc_linkable');
