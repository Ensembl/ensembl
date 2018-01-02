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

# patch_41_42_b
#
# title: unconventional transcripts
#
# description:
# Add new table to support unconventional transcripts
# Also allow transcripts to not have genes (remove NOT NULL constraint on transcript.gene_id)

ALTER TABLE transcript CHANGE COLUMN gene_id gene_id INT(10) UNSIGNED;

CREATE TABLE unconventional_transcript_association (

  transcript_id    INT(10) UNSIGNED NOT NULL,
  gene_id          INT(10) UNSIGNED NOT NULL,
  interaction_type ENUM("antisense","sense_intronic","sense_overlaping_exonic","chimeric_sense_exonic"),

  KEY (transcript_id),
  KEY (gene_id)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;


# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_41_42_b.sql|unconventional_transcripts');
