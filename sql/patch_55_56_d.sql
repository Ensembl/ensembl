# patch_55_56_d.sql
#
# Title: Add an index to the splicing_event_feature table
#
# Description:
# With an index on transcript_id in splicing_event_feature, the
# generation of biomarts will be sped up.

ALTER TABLE splicing_event_feature
  ADD INDEX transcript_idx (transcript_id);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
VALUES (NULL, 'patch', 'patch_55_56_d.sql|add_index_to_splicing_event_feature');
