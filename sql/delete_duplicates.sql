CREATE PROCEDURE delete_duplicates()
BEGIN
 
DECLARE v_identifier varChar(255);
DECLARE v_unmapped_object_id INT;
DECLARE existing_object_id INT DEFAULT 0;
DECLARE v_parent varChar(255);
DECLARE v_ensembl_id INT;
DECLARE v_unmapped_reason_id INT;
DECLARE v_ensembl_object_type varChar(15);
DECLARE v_external_db_id INT;
DECLARE done INT DEFAULT 0;
DECLARE row_count INT DEFAULT 0;

DECLARE cur1 CURSOR FOR SELECT unmapped_object_id, identifier, ensembl_id, parent, unmapped_reason_id, ensembl_object_type, external_db_id from unmapped_object where ensembl_id is not null and parent is not null and ensembl_object_type is not null and external_db_id is not null order by unmapped_object_id;
DECLARE CONTINUE HANDLER FOR NOT FOUND SET done = TRUE;

OPEN cur1;

main_loop: LOOP

	SET done = FALSE;
	FETCH cur1 INTO v_unmapped_object_id,v_identifier, v_ensembl_id, v_parent, v_unmapped_reason_id, v_ensembl_object_type, v_external_db_id;

	IF done THEN
   		LEAVE main_loop;
	END IF;

	SET row_count = 0;
	SELECT COUNT(1) INTO row_count FROM unmapped_object WHERE identifier = v_identifier and ensembl_id = v_ensembl_id and parent = v_parent and unmapped_reason_id = v_unmapped_reason_id and ensembl_object_type = v_ensembl_object_type and external_db_id = v_external_db_id;
	IF row_count > 1 THEN
		DELETE FROM unmapped_object where unmapped_object_id <> v_unmapped_object_id and identifier = v_identifier and ensembl_id = v_ensembl_id and parent = v_parent and unmapped_reason_id = v_unmapped_reason_id and ensembl_object_type = v_ensembl_object_type and external_db_id = v_external_db_id;

	END IF;
END LOOP;

CLOSE cur1;

END

