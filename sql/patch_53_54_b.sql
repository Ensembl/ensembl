# patch_53_54_b.sql
#
# title: Widen some text columns
#
# description:
# Change oligo_probe.name to 40 characters, external_db.db_name to 100 chars, analysis.logic_name to 128 characters

ALTER TABLE oligo_probe MODIFY name VARCHAR(40);

ALTER TABLE external_db MODIFY db_name VARCHAR(100) NOT NULL;

ALTER TABLE analysis MODIFY logic_name VARCHAR(128) NOT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_53_54_b.sql|widen_columns');


