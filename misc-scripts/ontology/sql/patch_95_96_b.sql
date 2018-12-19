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

-- patch_95_96_b.sql
--
-- Title: Update schema version.
--
-- Description:
--  Added fields to Term and Ontology tables

-- Added columns
ALTER TABLE `ontology` ADD COLUMN `title` varchar(255) DEFAULT NULL;

ALTER TABLE `term` ADD COLUMN `iri` text NULL;

ALTER TABLE `meta` DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

ALTER TABLE `ontology` DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

ALTER TABLE `subset` DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

ALTER TABLE `term` DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

ALTER TABLE `synonym` DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

CREATE INDEX ix_closure_subparent_term_id ON closure (`subparent_term_id`);
CREATE INDEX ix_closure_ontology_id ON closure (`ontology_id`);

CREATE INDEX `ix_relation_parent_term_id` ON relation (`parent_term_id`);
CREATE INDEX `ix_relation_relation_type_id` ON relation (`relation_type_id`);
CREATE INDEX `ix_relation_ontology_id` ON relation (`ontology_id`);