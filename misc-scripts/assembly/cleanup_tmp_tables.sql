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

#this is a SQL file intended to be run after mapping the assemblies
#once you are happy with the results, run it like this
#
#mysql -h your_host -u the_user -p database_updated < cleanup_tmp_table.sql
#will remove the *bak and *tmp* tables created during the mapping process

DROP TABLE assembly_bak;
DROP TABLE coord_system_bak;
DROP TABLE meta_bak;
DROP TABLE seq_region_bak;
DROP TABLE tmp_align;
