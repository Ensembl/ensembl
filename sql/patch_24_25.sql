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

# Add and populate created_date and modified_date to stable_id tables
# Schema 24-25

set @today = concat( curdate(), " 12:00:00" );

ALTER TABLE gene_stable_id ADD COLUMN created_date DATETIME NOT NULL;
ALTER TABLE gene_stable_id ADD COLUMN modified_date DATETIME NOT NULL;
UPDATE gene_stable_id SET created_date=@today;
UPDATE gene_stable_id SET modified_date=@today;

ALTER TABLE exon_stable_id ADD COLUMN created_date DATETIME NOT NULL;
ALTER TABLE exon_stable_id ADD COLUMN modified_date DATETIME NOT NULL;
UPDATE exon_stable_id SET created_date=@today;
UPDATE exon_stable_id SET modified_date=@today;

ALTER TABLE transcript_stable_id ADD COLUMN created_date DATETIME NOT NULL;
ALTER TABLE transcript_stable_id ADD COLUMN modified_date DATETIME NOT NULL;
UPDATE transcript_stable_id SET created_date=@today;
UPDATE transcript_stable_id SET modified_date=@today;

ALTER TABLE translation_stable_id ADD COLUMN created_date DATETIME NOT NULL;
ALTER TABLE translation_stable_id ADD COLUMN modified_date DATETIME NOT NULL;
UPDATE translation_stable_id SET created_date=@today;
UPDATE translation_stable_id SET modified_date=@today;
