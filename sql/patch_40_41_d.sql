# patch_40_41_d
#
# title: ditag primary key type
#
# description:
# Make ditag.ditag_id be INT(10)

ALTER TABLE ditag CHANGE COLUMN ditag_id ditag_id INT(10) UNSIGNED;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_40_41_d.sql|ditag_primary_key_type');


