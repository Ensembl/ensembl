
# patch_49_50_e.sql
#
# title: Create tables for new mapping_seq_region tables
#
# Description: Creates two new tables, seq_region_mapping and mapping_seq that will allow to upload user data and be able to map it
# even there is a change of seq_region in a previous release

################################################################################
#
# Table structure for seq_region mapping between releases
#
# Stores how the core seq_region_id have changed from release to release

CREATE TABLE seq_region_mapping (

	external_seq_region_id	INT(10) UNSIGNED NOT NULL,
	internal_seq_region_id	INT(10) UNSIGNED NOT NULL,
	mapping_set_id		INT(10) UNSIGNED NOT NULL,

	KEY (mapping_set_id)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;

################################################################################
#
# Table structure for seq_region mapping between releases
#
# Stores how which mapping group the seq_region are for a particular schema

CREATE TABLE mapping_set (

	mapping_set_id	INT(10)	UNSIGNED NOT NULL,
	schema_build	VARCHAR(20) NOT NULL,

	PRIMARY KEY(schema_build)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;

# Patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_49_50_e.sql|mapping_seq_region');



