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


################################################################################

# change the indexes on the align feature tables so that range queries
# can use the indexes properly

ALTER TABLE dna_align_feature DROP INDEX seq_region_idx;
ALTER TABLE protein_align_feature DROP INDEX seq_region_idx;
ALTER TABLE dna_align_feature ADD INDEX seq_region_idx( seq_region_id, analysis_id, seq_region_start, score);
ALTER TABLE protein_align_feature ADD INDEX seq_region_idx( seq_region_id, analysis_id, seq_region_start, score);


# reconstruct the meta coord table with the addition of a 
# max feature length column

DELETE FROM meta_coord;
ALTER TABLE meta_coord add column max_length int;

INSERT INTO meta_coord
SELECT 'assembly_exception', sr.coord_system_id, 
       MAX(IF(ae.seq_region_end - ae.seq_region_start > ae.exc_seq_region_end - ae.exc_seq_region_start, ae.seq_region_end - ae.seq_region_start + 1, ae.exc_seq_region_end - ae.exc_seq_region_start + 1))
FROM   assembly_exception ae, seq_region sr
WHERE  sr.seq_region_id = ae.seq_region_id
GROUP BY sr.coord_system_id;

INSERT INTO meta_coord
SELECT 'density_feature', sr.coord_system_id,
       MAX(df.seq_region_end - df.seq_region_start + 1)
FROM   density_feature df, seq_region sr
WHERE  sr.seq_region_id = df.seq_region_id
GROUP BY sr.coord_system_id;


INSERT INTO meta_coord
SELECT 'dna_align_feature', sr.coord_system_id,
       MAX(daf.seq_region_end - daf.seq_region_start + 1)
FROM   dna_align_feature daf, seq_region sr
WHERE  sr.seq_region_id = daf.seq_region_id
GROUP BY sr.coord_system_id;


INSERT INTO meta_coord
SELECT 'exon', sr.coord_system_id,
       MAX(e.seq_region_end - e.seq_region_start + 1)
FROM   exon e, seq_region sr
WHERE  sr.seq_region_id = e.seq_region_id
GROUP BY sr.coord_system_id;


INSERT INTO meta_coord
SELECT 'gene', sr.coord_system_id,
       MAX(g.seq_region_end - g.seq_region_start + 1)
FROM   gene g, seq_region sr
WHERE  sr.seq_region_id = g.seq_region_id
GROUP BY sr.coord_system_id;


INSERT INTO meta_coord
SELECT 'karyotype', sr.coord_system_id,
       MAX(k.seq_region_end - k.seq_region_start + 1)
FROM   karyotype k, seq_region sr
WHERE  sr.seq_region_id = k.seq_region_id
GROUP BY sr.coord_system_id;


INSERT INTO meta_coord
SELECT 'marker_feature', sr.coord_system_id,
       MAX(mf.seq_region_end - mf.seq_region_start + 1)
FROM   marker_feature mf, seq_region sr
WHERE  sr.seq_region_id = mf.seq_region_id
GROUP BY sr.coord_system_id;


INSERT INTO meta_coord
SELECT 'misc_feature', sr.coord_system_id,
       MAX(mf.seq_region_end - mf.seq_region_start + 1)
FROM   misc_feature mf, seq_region sr
WHERE  sr.seq_region_id = mf.seq_region_id
GROUP BY sr.coord_system_id;


INSERT INTO meta_coord
SELECT 'prediction_exon', sr.coord_system_id,
       MAX(pe.seq_region_end - pe.seq_region_start + 1)
FROM   prediction_exon pe, seq_region sr
WHERE  sr.seq_region_id = pe.seq_region_id
GROUP BY sr.coord_system_id;


INSERT INTO meta_coord
SELECT 'prediction_transcript', sr.coord_system_id,
       MAX(pt.seq_region_end - pt.seq_region_start + 1)
FROM   prediction_transcript pt, seq_region sr
WHERE  sr.seq_region_id = pt.seq_region_id
GROUP BY sr.coord_system_id;


INSERT INTO meta_coord
SELECT 'protein_align_feature', sr.coord_system_id,
       MAX(paf.seq_region_end - paf.seq_region_start + 1)
FROM   protein_align_feature paf, seq_region sr
WHERE  sr.seq_region_id = paf.seq_region_id
GROUP BY sr.coord_system_id;

INSERT INTO meta_coord
SELECT 'qtl_feature', sr.coord_system_id,
       MAX(qf.seq_region_end - qf.seq_region_start + 1)
FROM   qtl_feature qf, seq_region sr
WHERE  sr.seq_region_id = qf.seq_region_id
GROUP BY sr.coord_system_id;


INSERT INTO meta_coord
SELECT 'repeat_feature', sr.coord_system_id,
       MAX(rf.seq_region_end - rf.seq_region_start + 1)
FROM   repeat_feature rf, seq_region sr
WHERE  sr.seq_region_id = rf.seq_region_id
GROUP BY sr.coord_system_id;


INSERT INTO meta_coord
SELECT 'simple_feature', sr.coord_system_id,
       MAX(sf.seq_region_end - sf.seq_region_start + 1)
FROM   simple_feature sf, seq_region sr
WHERE  sr.seq_region_id = sf.seq_region_id
GROUP BY sr.coord_system_id;


INSERT INTO meta_coord
SELECT 'transcript', sr.coord_system_id,
       MAX(t.seq_region_end - t.seq_region_start + 1)
FROM   transcript t, seq_region sr
WHERE  sr.seq_region_id = t.seq_region_id
GROUP BY sr.coord_system_id;

################################################################################

# Add another index on the align feature tables that works without analysis id
ALTER TABLE dna_align_feature ADD INDEX seq_region_idx_2( seq_region_id, seq_region_start);
ALTER TABLE protein_align_feature ADD INDEX seq_region_idx_2( seq_region_id, seq_region_start);
