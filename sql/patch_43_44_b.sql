# patch_43_44_b
#
# title: optimising ditag tables
#
# description:
# some small memory & speed improvements for the ditag tables

ALTER TABLE ditag CHANGE COLUMN type type varchar(30) NOT NULL, \
 CHANGE COLUMN tag_count tag_count smallint(6) unsigned NOT NULL default 1, \
 CHANGE COLUMN sequence sequence TINYTEXT NOT NULL, \
 CHANGE COLUMN name name varchar(30) NOT NULL;

ALTER TABLE ditag_feature CHANGE COLUMN cigar_line cigar_line TINYTEXT NOT NULL, \
 CHANGE COLUMN ditag_side ditag_side ENUM('F', 'L', 'R') NOT NULL;
CREATE INDEX seq_region_idx ON ditag_feature (seq_region_id, seq_region_start, seq_region_end);

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_43_44_b.sql|optimise_ditag_tables');

