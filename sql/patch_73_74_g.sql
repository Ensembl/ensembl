# patch_73_74_g.sql
#
# Title: Adding transcript index to transcript_intron_supporting_evidence
#
# Description:
#
# Adding an index on transcript id to transcript_intron_supporting_evidence to
# speed up retrieval of supporting features from a Transcript object 

CREATE INDEX transcript_idx ON transcript_intron_supporting_evidence(transcript_id);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_73_74_g.sql|add_transcript_idx_tise');

 
