# patch_66_67_b.sql
#
# Title: Drop stable ID views.
#
# Description:
# The stable ID views, introduced for release 65 as a way of providing a
# degree of backward compatibility, are dropped with this release.

-- Heavens!

DROP VIEW IF EXISTS exon_stable_id;
DROP VIEW IF EXISTS gene_stable_id;
DROP VIEW IF EXISTS operon_stable_id;
DROP VIEW IF EXISTS operon_transcript_stable_id;
DROP VIEW IF EXISTS translation_stable_id;
DROP VIEW IF EXISTS transcript_stable_id;

DROP TABLE IF EXISTS exon_stable_id;
DROP TABLE IF EXISTS gene_stable_id;
DROP TABLE IF EXISTS operon_stable_id;
DROP TABLE IF EXISTS operon_transcript_stable_id;
DROP TABLE IF EXISTS translation_stable_id;
DROP TABLE IF EXISTS transcript_stable_id;



# Patch identifier:
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_66_67_b.sql|drop_stable_id_views');
