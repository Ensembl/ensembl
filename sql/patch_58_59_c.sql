# patch_58_59_c.sql
#
# Title:
#   Replace splicing_event.type with splicing_event.attrib_type_id
#
# Description:
#   The 'type' enumeration in the splicing_event table is too terse.
#   Replace it with a reference to a proper attrib_type.

# Modify the splicing_event table.
ALTER TABLE splicing_event
  DROP COLUMN `type`,
  ADD COLUMN attrib_type_id SMALLINT(5) UNSIGNED NOT NULL DEFAULT 0
    AFTER seq_region_strand;

# Insert patch identifier.
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch',
    'patch_58_59_c.sql|splicing_event_attrib_type_id');
