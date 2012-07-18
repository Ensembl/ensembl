# patch_67_68_e.sql
#
# Title: Xref unique index fix
#
# Description: Alters the table xref to turn on the unique index. The uniquness
#              required has already been applied by patch b. This is ensuring
#              users are synchronized with the 67 schema and we do so with
#              a stored procedure which can optionally add the index if it was 
#              missing
# 

DELIMITER $$
DROP PROCEDURE IF EXISTS `create_index_if_not_exists`$$

CREATE PROCEDURE `create_index_if_not_exists`(table_schema_vc varchar(64))
SQL SECURITY INVOKER
BEGIN

set @Index_version_count = (
select count(1)
FROM  INFORMATION_SCHEMA.STATISTICS
WHERE table_name = 'xref'
and index_name = 'id_index'
and TABLE_SCHEMA = table_schema_vc
and COLUMN_NAME = 'version'
);

IF @Index_version_count <> 0 THEN
  PREPARE idx_remove from 'alter table xref drop index id_index;';
  EXECUTE idx_remove;
  DEALLOCATE PREPARE idx_remove;
  
	PREPARE stmt FROM 'Alter table xref ADD UNIQUE INDEX id_index (dbprimary_acc,external_db_id,info_type,info_text,version);';
	EXECUTE stmt;
	DEALLOCATE PREPARE stmt;
END IF;


END$$
DELIMITER ;

call create_index_if_not_exists(database());

DROP PROCEDURE `create_index_if_not_exists`;

# Allows us to do the following without creating new indexes without good reason
# ALTER TABLE xref DROP KEY id_index;
# ALTER TABLE xref ADD UNIQUE KEY id_index (dbprimary_acc, version, external_db_id, info_type, info_text);

# Patch identifier:
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_67_68_e.sql|fix_67_68_e_xref_index');