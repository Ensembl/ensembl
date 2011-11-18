# patch_62_63_c.sql
#
# Title: New table for file location storage
#
# Description:
# We now record the location of an external file to be used by the website
# for display. This is primarily for BAM data but we will extend this
# for other formats.

CREATE TABLE data_file (
	data_file_id int(11) unsigned NOT NULL AUTO_INCREMENT,
	coord_system_id int(11) NOT NULL,
	analysis_id int(11) NOT NULL,
	name varchar(100) NOT NULL,
	version_lock tinyint(1) DEFAULT 0 NOT NULL,
	absolute tinyint(1) DEFAULT 0 NOT NULL,
	url text,
	file_type enum('BAM','BIGBED','BIGWIG','VCF'),
	PRIMARY KEY (data_file_id),
  UNIQUE KEY df_unq_idx(coord_system_id, analysis_id, name, file_type),
  INDEX df_name_idx(name),
  INDEX df_analysis_idx(analysis_id)
) ENGINE=MyISAM;

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_64_65_c.sql|add_data_file');
