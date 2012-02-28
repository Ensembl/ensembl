# patch_66_67_b.sql
#
# Title: Drop stable ID views.
#
# Description:
# The stable ID views, introduced for release 65 as a way of providing a
# degree of backward compatibility, are dropped with this release.

DROP VIEW exon_stable_id;
DROP VIEW gene_stable_id;
DROP VIEW operon_stable_id;
DROP VIEW operon_transcript_stable_id;
DROP VIEW translation_stable_id;
DROP VIEW transcript_stable_id;

# Patch identifier:
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_66_67_b.sql|drop_stable_id_views');
