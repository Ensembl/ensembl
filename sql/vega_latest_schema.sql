###################################################################
# FROM patch_23_24.sql	

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



###################################################################
# FROM patch_30_31.sql


# gene table changes

ALTER TABLE gene CHANGE type biotype VARCHAR(40) NOT NULL default 'protein_coding';
ALTER TABLE gene ADD source VARCHAR(20) NOT NULL default 'ensembl';
ALTER TABLE gene ADD confidence ENUM( 'KNOWN', 'NOVEL', 'PUTATIVE', 'PREDICTED' ) default 'NOVEL';
ALTER TABLE gene ADD description text;

UPDATE gene g, gene_description gd SET g.description = gd.description WHERE gd.gene_id = g.gene_id;

DROP TABLE gene_description;

# transcript related changes

ALTER TABLE transcript ADD biotype VARCHAR(40) NOT NULL DEFAULT 'protein_coding';
ALTER TABLE transcript ADD confidence ENUM( 'KNOWN', 'NOVEL', 'PUTATIVE', 'PREDICTED' ) default 'NOVEL';
ALTER TABLE transcript ADD description text;

# reasonable default for transcript description
# it might be questionable whether a separate transcript description is necessary
 
UPDATE transcript t, xref x SET t.description = x.description WHERE t.display_xref_id = x.xref_id;
UPDATE transcript SET description=NULL WHERE description="";

# usefull settings for the new tables

UPDATE gene SET source = 'vega';

#########################################################################################################
# may well want to update confidence for genes and transcripts based on source of the xref, but not yet #

#UPDATE gene g, xref x, external_db ed SET g.confidence='KNOWN' WHERE g.display_xref_id = x.xref_id and x.external_db_id = ed.external_db_id and g.display_xref_id != 0 and ed.status like 'KNOWN%';
#UPDATE transcript t, xref x, external_db ed SET t.confidence='KNOWN' WHERE t.display_xref_id = x.xref_id and x.external_db_id = ed.external_db_id and t.display_xref_id != 0 and ed.status like 'KNOWN%';

# HAWK -> biotype conversions

UPDATE gene SET biotype='processed_transcript' WHERE biotype = 'Transcript';
UPDATE gene SET biotype='pseudogene', confidence='NOVEL' WHERE biotype = 'Pseudogene';
UPDATE gene SET biotype='protein_coding', confidence='NOVEL' WHERE biotype = 'Novel_CDS';
UPDATE gene SET biotype='processed_transcript', confidence='NOVEL' WHERE biotype = 'Novel_Transcript';
UPDATE gene SET biotype='processed_transcript',confidence='PUTATIVE' WHERE biotype = 'Putative';
UPDATE gene SET biotype='protein_coding', confidence='KNOWN' WHERE biotype = 'Known';

UPDATE gene SET biotype='processed_pseudogene', confidence=NULL WHERE biotype = 'Processed_pseudogene';
UPDATE gene SET biotype='unprocessed_pseudogene', confidence=NULL WHERE biotype = 'Unprocessed_pseudogene';
UPDATE gene SET biotype='protein_coding',confidence='PREDICTED' WHERE biotype = 'Predicted_Gene';
UPDATE gene SET biotype='Ig_segment', confidence=NULL WHERE biotype = 'Ig_Segment';
UPDATE gene SET biotype='Ig_pseudogene_segment', confidence=NULL WHERE biotype = 'Ig_Pseudogene_Segment';

UPDATE gene SET biotype=replace( biotype, '-','_' );

# reasonable biotypes for the transcripts, take the one from the gene

UPDATE transcript t, gene g SET t.biotype = g.biotype, t.confidence = g.confidence WHERE g.gene_id = t.gene_id;


########################
# FROM patch_31_32.sql #
########################


# Add linkable & priority columns to external_db
ALTER TABLE external_db ADD COLUMN dbprimary_acc_linkable BOOLEAN DEFAULT 1 NOT NULL;
ALTER TABLE external_db ADD COLUMN display_label_linkable BOOLEAN DEFAULT 0 NOT NULL;
ALTER TABLE external_db ADD COLUMN priority INT NOT NULL;


########################
# FROM patch_32_33.sql #
########################

# Add db_display_name column to external_db
ALTER TABLE external_db ADD COLUMN db_display_name VARCHAR(255);

########################
# FROM patch_33_34.sql #
########################

ALTER table object_xref MODIFY ensembl_object_type ENUM( 'RawContig', 'Transcript', 'Gene', 'Translation', 'regulatory_factor', 'regulatory_feature' ) not null;

alter table transcript change confidence status  enum( 'KNOWN', 'NOVEL', 'PUTATIVE', 'PREDICTED' );
alter table gene change confidence status  enum( 'KNOWN', 'NOVEL', 'PUTATIVE', 'PREDICTED' );

########################
# FROM patch_36_37.sql #
########################

#manipulations of peptide and gene archive tables are in ensembl patch scripts
#but are not used for vega since tables are empty




###########################
#set latest schema version#
###########################
update meta set meta_value = '37' where meta_key = 'schema_version';
