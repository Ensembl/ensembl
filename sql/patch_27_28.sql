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

# Patch SQL required to convert Ensembl version 27 schemas to version 28

# Now have 2 separate RefSeq external_db entries
UPDATE external_db SET db_name='RefSeq_dna' WHERE db_name='RefSeq';

# Add Flybase external_db entries
INSERT INTO external_db VALUES (803, 'flybase_polypeptide_id', 1, 'KNOWNXREF');
INSERT INTO external_db VALUES (805, 'flybase_annotation_id',  1, 'KNOWNXREF');

# Modify column definition of xref.dbprimary_acc to remove binary property
ALTER TABLE xref MODIFY dbprimary_acc VARCHAR(40) NOT NULL;

# Add analysis_description table
CREATE TABLE analysis_description (
  analysis_id	               int(10) unsigned NOT NULL,
  description                  text,
  display_label                varchar(255),

  KEY analysis_idx( analysis_id )
);
