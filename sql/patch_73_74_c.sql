# patch_73_74_c.sql
#
# title: Unconventional_transcript_association removal
#
# description:
# Removal unconventional_transcript_association table which is not used any more

DROP TABLE unconventional_transcript_association;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_73_74_c.sql|remove_unconventional_transcript_association');

 
