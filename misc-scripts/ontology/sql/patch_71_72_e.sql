-- patch_71_72_e.sql
--
-- Title: Add is_obsolete
--
-- Description:
--   Adds the is_obsolete flag to the term table

ALTER TABLE TERM
ADD COLUMN is_obsolete INT NOT NULL DEFAULT 0;

-- Patch identifier
INSERT INTO meta (meta_key, meta_value)
  VALUES ('patch', 'patch_71_72_e.sql|is_obsolete');


