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

# patch_62_63_b.sql
#
# Title: Indexing changes for core database.
#
# Description:

#change stable Id version to not null, default 1

ALTER TABLE exon_stable_id MODIFY version INT(10) NOT NULL DEFAULT 1; 

ALTER TABLE gene_stable_id MODIFY version INT(10) NOT NULL DEFAULT 1; 

ALTER TABLE transcript_stable_id MODIFY version INT(10) NOT NULL DEFAULT 1;  

ALTER TABLE translation_stable_id MODIFY version INT(10) NOT NULL DEFAULT 1;  

ALTER TABLE gene_archive MODIFY gene_version SMALLINT NOT NULL DEFAULT 1, MODIFY transcript_version SMALLINT NOT NULL DEFAULT 1, MODIFY translation_version SMALLINT NOT NULL DEFAULT 1;

DROP INDEX  gene_idx ON gene_archive;

CREATE INDEX  gene_idx ON gene_archive(gene_stable_id, gene_version);

DROP INDEX  transcript_idx ON gene_archive;

CREATE INDEX  transcript_idx ON gene_archive(transcript_stable_id, transcript_version);

DROP INDEX  translation_idx ON gene_archive;
 
CREATE INDEX  translation_idx ON gene_archive(translation_stable_id, translation_version);


#umapped_object new unique index

DROP INDEX anal_idx ON unmapped_object;

CREATE UNIQUE INDEX unique_unmapped_obj_idx ON unmapped_object(identifier, ensembl_id, parent, unmapped_reason_id, ensembl_object_type, external_db_id);

DROP INDEX id_idx ON unmapped_object;

#reduce identifier index to 50 characters – most ids are under 30, index will be faster

CREATE INDEX id_idx ON unmapped_object(identifier(50));

#add index for queries using external_db_id in the where clause

CREATE INDEX ext_db_identifier_idx ON unmapped_object(external_db_id, identifier);


# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_62_63_b.sql|indexing_changes');




