CREATE TABLE affy_feature (
       affy_feature_id INT NOT NULL auto_increment,
       seq_region_id INT UNSIGNED NOT NULL,
       seq_region_start INT NOT NULL,
       seq_region_end INT NOT NULL,
       seq_region_strand TINYINT NOT NULL,
       
       mismatches TINYINT,
       affy_probe_id INT NOT NULL,
       analysis_id INT NOT NULL,

       PRIMARY KEY (affy_feature_id),
       KEY seq_region_idx( seq_region_id, seq_region_start ),
       KEY probe_idx( affy_probe_id )
);


CREATE TABLE affy_probe (
       affy_probe_id INT NOT NULL auto_increment,
       affy_array_id INT NOT NULL,
       probeset VARCHAR(20),
       name VARCHAR(20),

       PRIMARY KEY ( affy_probe_id, affy_array_id ),
       KEY probeset_idx( probeset ),
       KEY array_idx( affy_array_id )
);


CREATE TABLE affy_array (
       affy_array_id INT NOT NULL auto_increment,
       parent_array_id INT,
       probe_setsize TINYINT NOT NULL,
       name VARCHAR(40) NOT NULL,

       PRIMARY KEY( affy_array_id )
);

