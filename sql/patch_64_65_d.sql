# patch_62_63_d.sql
#
# Title: Add new Xref info_type
#
# Description:
# We can now do checksum mapping in the Xref pipeline so Xrefs need to reflect
# this

ALTER TABLE xref modify COLUMN info_type ENUM( 
                                    'PROJECTION', 'MISC', 'DEPENDENT',
                                    'DIRECT', 'SEQUENCE_MATCH',
                                    'INFERRED_PAIR', 'PROBE',
                                    'UNMAPPED', 'COORDINATE_OVERLAP',
                                    'CHECKSUM');

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_64_65_d.sql|add_checksum_info_type');
