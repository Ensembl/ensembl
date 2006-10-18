# patch_40_41_e
#
# title: ditag correction from patch_40_41_d
#
# description: put auto increment back in the table

alter table ditag change column ditag_id ditag_id int(10) unsigned  NOT NULL auto_increment;

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_41_42_e.sql|ditag_autoincrement');
