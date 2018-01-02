-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2018] EMBL-European Bioinformatics Institute
-- 
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
-- 
--      http://www.apache.org/licenses/LICENSE-2.0
-- 
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

# patch_70_71_b.sql
#
# Title: Mapping_set unique index fix
#
# Description: Alters the table mapping_set to turn on the unique index.
#              This is ensuring users are synchronized with the 70 schema
#              and we do so with a stored procedure which can optionally
#              add the index if it was missing
# 

DELIMITER $$
DROP PROCEDURE IF EXISTS `create_index_if_not_exists`$$

CREATE PROCEDURE `create_index_if_not_exists`(table_schema_vc varchar(64))
SQL SECURITY INVOKER
BEGIN

set @Index_version_count = (
select count(1)
FROM  INFORMATION_SCHEMA.STATISTICS
WHERE table_name = 'mapping_set'
and index_name = 'PRIMARY'
and TABLE_SCHEMA = table_schema_vc
and COLUMN_NAME = 'mapping_set_id'
);

IF @Index_version_count <> 0 THEN
  PREPARE idx_remove from 'alter table mapping_set drop primary key;';
  EXECUTE idx_remove;
  DEALLOCATE PREPARE idx_remove;
  
END IF;

PREPARE stmt FROM 'Alter table mapping_set ADD PRIMARY KEY (mapping_set_id);';
EXECUTE stmt;
DEALLOCATE PREPARE stmt;

END$$
DELIMITER ;

call create_index_if_not_exists(database());

DROP PROCEDURE `create_index_if_not_exists`;

# Allows us to do the following without creating new indexes without good reason
# ALTER TABLE mapping_set DROP PRIMARY KEY;
# ALTER TABLE mapping_set ADD PRIMARY KEY (mapping_set_id);

# Patch identifier:
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_70_71_b.sql|mapping_set_index');
