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

# patch_54_55_d.sql
#
# title: Add table to store the dependent xrefs
#
# description:
# For a given object_xref store its master if it is dependent on another xref (master)

CREATE TABLE dependent_xref(
     object_xref_id         INT NOT NULL,
     master_xref_id         INT NOT NULL,
     dependent_xref_id      INT NOT NULL,

     PRIMARY KEY( object_xref_id ),
     KEY dependent ( dependent_xref_id ),
     KEY master_idx (master_xref_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_d.sql|add_dependent_xref_table');

