# patch_66_67_e.sql
#
# Title: Adding an index to gene canonical transcript id
#
# Description: This is a lookup which is normally fast but when we have a lot
# of genes in a core schema (multispecies DBs) then these queries slow down 
# 

ALTER TABLE gene ADD KEY canonical_transcript_id_idx (canonical_transcript_id);

# Patch identifier:
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_66_67_e.sql|index_canonical_transcript_id');
