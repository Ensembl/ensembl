# patch_48_49_d.sql
#
# Title: Add new info_type to xref table
#
# Description:
# Add the ENUM 'COORDINATE_OVERLAP' to the xref.info_type column.

ALTER TABLE xref CHANGE COLUMN info_type
  info_type ENUM(
    'PROJECTION', 'MISC', 'DEPENDENT', 'DIRECT', 'SEQUENCE_MATCH',
    'INFERRED_PAIR', 'PROBE', 'UNMAPPED', 'COORDINATE_OVERLAP'
  );

# Patch identifier
INSERT INTO meta (meta_key, meta_value)
VALUES ('patch', 'patch_48_49_d.sql|new_info_type_enum');
