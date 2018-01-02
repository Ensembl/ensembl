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

#
# some useful stable ID mapping related SQL
# just fragments for copy/paste, not intended to be run as an "SQL script"
#

# get counts from all stable ID related tables
SELECT COUNT(*) FROM gene_stable_id;
SELECT COUNT(*) FROM transcript_stable_id;
SELECT COUNT(*) FROM translation_stable_id;
SELECT COUNT(*) FROM exon_stable_id;
SELECT COUNT(*) FROM mapping_session;
SELECT COUNT(*) FROM stable_id_event;
SELECT COUNT(*) FROM gene_archive;
SELECT COUNT(*) FROM peptide_archive;

# backup all stable ID related tables
CREATE TABLE gene_stable_id_bak SELECT * FROM gene_stable_id;
CREATE TABLE transcript_stable_id_bak SELECT * FROM transcript_stable_id;
CREATE TABLE translation_stable_id_bak SELECT * FROM translation_stable_id;
CREATE TABLE exon_stable_id_bak SELECT * FROM exon_stable_id;
CREATE TABLE mapping_session_bak SELECT * FROM mapping_session;
CREATE TABLE stable_id_event_bak SELECT * FROM stable_id_event;
CREATE TABLE gene_archive_bak SELECT * FROM gene_archive;
CREATE TABLE peptide_archive_bak SELECT * FROM peptide_archive;

# now prune all of them
DELETE FROM gene_stable_id;
DELETE FROM transcript_stable_id;
DELETE FROM translation_stable_id;
DELETE FROM exon_stable_id;
DELETE FROM mapping_session;
DELETE FROM stable_id_event;
DELETE FROM gene_archive;
DELETE FROM peptide_archive;

# drop backup tables
DROP TABLE gene_stable_id_bak;
DROP TABLE transcript_stable_id_bak;
DROP TABLE translation_stable_id_bak;
DROP TABLE exon_stable_id_bak;
DROP TABLE mapping_session_bak;
DROP TABLE stable_id_event_bak;
DROP TABLE gene_archive_bak;
DROP TABLE peptide_archive_bak;

