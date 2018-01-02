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

# changes to database structure from 26 to 27
# none of which should affect your current data
#  so you dont need to apply them


ALTER TABLE xref MODIFY dbprimary_acc VARCHAR(40) BINARY NOT NULL;
ALTER TABLE affy_probe MODIFY probeset VARCHAR(40);
ALTER TABLE interpro DROP INDEX interpro_ac;
ALTER TABLE interpro DROP INDEX id;

ALTER TABLE interpro ADD UNIQUE ( interpro_ac, id );
ALTER TABLE interpro ADD INDEX( id );



