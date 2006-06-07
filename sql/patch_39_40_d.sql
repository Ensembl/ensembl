# patch_39_40_d
#
# title: stable id event score column
#
# description:
# Add stable_id_event.score column

ALTER TABLE stable_id_event ADD COLUMN score FLOAT NOT NULL DEFAULT 0;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_39_40_d.sql|stable_id_event_score_column');

























































# Similarly make all seq_region_start & seq_region_end columns INT(10) UNSIGNED
# Although these are not used as keys, they are used in many joins and having
# the same column type should make joins faster as MySQL will not have to do
# any type casting.








