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

# patch_45_46_b.sql
#
# title: go_xref.source_xref_id
#
# description:
#   Addition of a source_xref_id field to the go_xref table

ALTER TABLE go_xref ADD COLUMN source_xref_id int(10) unsigned default NULL;

ALTER TABLE go_xref DROP KEY object_xref_id_2;

ALTER TABLE go_xref ADD UNIQUE KEY object_xref_id_2 
  ( object_xref_id, linkage_type, source_xref_id );

ALTER TABLE go_xref ADD KEY source_xref_id( source_xref_id );


# patch identifier
INSERT INTO meta (meta_key, meta_value) 
  VALUES ('patch', 'patch_45_46_b.sql|go_xref.source_xref_id');



