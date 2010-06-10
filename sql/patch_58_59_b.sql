# patch_58_59_b.sql
#
# Title:
#   Modify the assembly_exception.exc_type enumeration.
#
# Description:
#   Add 'PATCH_FIX' and 'PATCH_NOVEL' to the assembly_exception.exc_type
#   enumeration.

# Modify the assembly_exception table.
ALTER TABLE assembly_exception MODIFY exc_type
  ENUM('HAP', 'PAR', 'PATCH_FIX', 'PATCH_NOVEL') NOT NULL;

# Insert patch identifier.
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch',
    'patch_58_59_b.sql|assembly_exception_exc_type_enum');
