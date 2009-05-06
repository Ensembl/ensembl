# patch_54_55_b.sql
#
# title: Add missing GO_XREF types
#
# description:
# Add missing GO_XREF enums that are missing

ALTER TABLE go_xref MODIFY COLUMN linkage_type enum('IC','IDA','IEA','IEP','IGI','IMP','IPI','ISS','NAS','ND','TAS','NR','RCA','EXP','ISO','ISA','ISM','IGC');

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_b.sql|add_go_xrefs_types');

 