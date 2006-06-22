# patch_39_40_f
#
# title: remove all_latest
#
# description:
# This patch removes the ALL/LATEST mapping session, and any associated stable_id_events.
# It is a data-only patch and does not involve a schema change.

DELETE s, m FROM stable_id_event s, mapping_session m WHERE m.mapping_session_id=s.mapping_session_id AND m.old_db_name='ALL' AND m.new_db_name='LATEST';

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_39_40_f.sql|remove_all_latest');

