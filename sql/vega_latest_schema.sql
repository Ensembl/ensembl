###################################################################
#FROM patch_23_24.sql	

# Adds a display_label column to the prediction_transcript table and
# populates it.
# This is done by creating a temporary table and renaming it rather than just
# using an alter table statement followed by an update statement.  This
# is because in mysql 3 it is not possible to use a join in an update 
# statement.


# create the temporary table

CREATE TABLE tmp_prediction_transcript (
    prediction_transcript_id int unsigned not null auto_increment,
    seq_region_id int unsigned not null,
    seq_region_start int unsigned not null,
    seq_region_end int unsigned not null,
    seq_region_strand tinyint not null,
    analysis_id int,
    display_label varchar(255),
    
    PRIMARY KEY( prediction_transcript_id ),
    KEY ( seq_region_id, seq_region_start ),
    KEY analysis_idx( analysis_id )
);


# populate the tmp table from the original table, and generate display_labels
# at the same time.  The display label generated looks like GENSCAN00000018082
# or SNAP00000018082, etc.

INSERT INTO tmp_prediction_transcript 
SELECT pt.prediction_transcript_id, pt.seq_region_id, pt.seq_region_start, 
       pt.seq_region_end, pt.seq_region_strand, pt.analysis_id,
       CONCAT(UPPER(a.logic_name), LPAD(pt.prediction_transcript_id, 11, '0'))
FROM   prediction_transcript pt, analysis a
WHERE  a.analysis_id = pt.analysis_id;


# drop the original table and rename the temp table

DROP TABLE prediction_transcript;
ALTER TABLE tmp_prediction_transcript RENAME prediction_transcript;
