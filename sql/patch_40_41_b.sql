# patch_40_41_b
#
# title: info_type enum
#
# description:
# Add more types to the enum in xref.info_type

ALTER TABLE xref CHANGE COLUMN info_type info_type ENUM('PROJECTION', 'MISC', 'DEPENDENT', 'DIRECT', 'SEQUENCE_MATCH', 'INFERRED_PAIR');

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_40_41_b.sql|info_type_enum');