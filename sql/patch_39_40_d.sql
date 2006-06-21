# patch_39_40_d
#
# title: stable id event score column
#
# description:
# Add stable_id_event.score column

ALTER TABLE stable_id_event ADD COLUMN score FLOAT NOT NULL DEFAULT 0;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_39_40_d.sql|stable_id_event_score_column');

