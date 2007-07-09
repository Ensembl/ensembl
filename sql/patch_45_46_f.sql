# patch_45_46_f.sql
#
# title stable_id_event.uni_idx
#
# description:
#   Drop old_version and new_version from uni_idx on stable_id_event

ALTER TABLE stable_id_event DROP KEY uni_idx;

ALTER TABLE stable_id_event ADD UNIQUE KEY uni_idx
  ( mapping_session_id, old_stable_id, new_stable_id, type );

# patch identifier
INSERT INTO meta (meta_key, meta_value) 
  VALUES ('patch', 'patch_45_46_f.sql|stable_id_event.uni_idx');



