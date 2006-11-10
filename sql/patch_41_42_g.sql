# patch_40_41_g
#
# title: Genebuild version format change
#
# description: Change format of genebuild.version entries in meta table.

UPDATE meta set meta_value = concat('20',substring(meta_value,1,2),'-', substring(meta_value,3,2),'-',substring(meta_value,5)) where meta_key = 'genebuild.version' and meta_value not rlike '^[0-9][0-9][0-9][0-9]-[0-9][0-9]-';

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_41_42_g.sql|genebuild_version_format_change');
