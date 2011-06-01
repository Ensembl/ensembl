CREATE PROCEDURE link_xrefs_to_genes()
BEGIN
 
DECLARE count_deleted, count_updated INT DEFAULT 0;
DECLARE v_object_xref_id INT;
DECLARE existing_object_xref INT DEFAULT 0;
DECLARE v_ensembl_id INT;
DECLARE v_ensembl_object_type varChar(15);
DECLARE v_xref_id INT;
DECLARE v_gene_id INT DEFAULT 0;
DECLARE done INT DEFAULT 0;

DECLARE cur1 CURSOR FOR SELECT object_xref_id, ensembl_id, ensembl_object_type, o.xref_id FROM object_xref o INNER JOIN xref x on o.xref_id = x.xref_id INNER JOIN external_db e on e.external_db_id =x.external_db_id WHERE db_name in ('UniGene', 'Uniprot_genename','DBASS3', 'DBASS5') AND ensembl_object_type <> 'Gene';
DECLARE CONTINUE HANDLER FOR NOT FOUND SET done = TRUE;

OPEN cur1;

main_loop: LOOP

	SET done = FALSE;
	FETCH cur1 INTO v_object_xref_id,v_ensembl_id, v_ensembl_object_type, v_xref_id;

	IF done THEN
   		LEAVE main_loop;
	END IF;
	-- find the corresponding gene_id
	SET v_gene_id = 0;
	CASE
		WHEN v_ensembl_object_type = 'Transcript' THEN
			SELECT COALESCE(gene_id,0) INTO v_gene_id FROM transcript where transcript_id = v_ensembl_id;
		WHEN v_ensembl_object_type = 'Translation' THEN
			SELECT COALESCE(gene_id,0) INTO v_gene_id FROM transcript tp inner join translation tl where tp.transcript_id = tl.transcript_id and translation_id = v_ensembl_id;
	END CASE;
	SET existing_object_xref = 0;
	IF v_gene_id > 0 THEN
    		SELECT COALESCE(object_xref_id,0) INTO existing_object_xref FROM object_xref where ensembl_id = v_gene_id AND xref_id = v_xref_id AND ensembl_object_type = 'Gene';
		IF existing_object_xref > 0 THEN
        		
			-- if gene_id already linked to this xref_id delete the object_xref
			DELETE FROM object_xref WHERE object_xref_id = v_object_xref_id;
			SET count_deleted = count_deleted + 1;
			DELETE FROM identity_xref WHERE object_xref_id = v_object_xref_id;
			DELETE FROM dependent_xref WHERE object_xref_id = v_object_xref_id;
    		ELSE
        		-- if not linked update the current object_xref
			UPDATE object_xref SET ensembl_id = v_gene_id, ensembl_object_type = 'Gene' WHERE object_xref_id = v_object_xref_id;
			SET count_updated = count_updated + 1;
    		END IF;
	END IF;
END LOOP;

CLOSE cur1;

-- populate the meta table with patch id and counts for updated rows
INSERT INTO meta (species_id, meta_key, meta_value) 
VALUES (NULL, 'data_changes', CONCAT('link_xrefs_to_genes.sql|updated ',count_updated,' xrefs, deleted ',count_deleted,' xrefs'));

END

