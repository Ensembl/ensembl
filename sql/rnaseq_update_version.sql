CREATE PROCEDURE rnaseq_version_update()
BEGIN
 
DECLARE v_count INT;
DECLARE v_stable_id varChar(128);
DECLARE v_version INT;
DECLARE v2_gene_id INT;
DECLARE v2_version INT;
DECLARE done INT DEFAULT 0;
DECLARE done2 INT DEFAULT 0;
DECLARE counter INT DEFAULT 0;

DECLARE cur1 CURSOR FOR select count(1) , stable_id, version from gene_stable_id group by stable_id, version having count(1) > 1 ;
DECLARE CONTINUE HANDLER FOR NOT FOUND SET done = TRUE;


OPEN cur1;


main_loop: LOOP

	SET done = FALSE;
	FETCH cur1 INTO v_count, v_stable_id, v_version;

	IF done THEN
   		LEAVE main_loop;
	END IF;

	DECLARE cur2 CURSOR FOR	select g.gene_id, gs.version  from seq_region s, gene_stable_id gs, gene g where g.gene_id=gs.gene_id and s.seq_region_id =g.seq_region_id and stable_id = v_stable_id;
	DECLARE CONTINUE HANDLER FOR NOT FOUND SET done = TRUE;
	OPEN cur2;
	counter = 0;
	loop2: LOOP2
	       SET done = FALSE;
	       FETCH cur2 INTO v2_gene_id, v2_version;

	       IF done2 THEN
   		  LEAVE loop2;
		  END IF;
		  SET counter = counter + 1;

	       IF (counter > 1) THEN
	       	  UPDATE gene_stable_id gs set version = (v2_version + counter - 1)  where gene_id = v2_gene_id;
	       END IF; 
		  
	END LOOP2;
	       CLOSE cur2;

END LOOP;

CLOSE cur1;

END

