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

# Patch to convert release 36 Ensembl schema to release 37

UPDATE meta set meta_value="37" where meta_key="schema_version";

# increase width of xref display_label column to allow for longer labels
ALTER TABLE xref CHANGE COLUMN display_label display_label VARCHAR(128) NOT NULL;

# update archives
ALTER TABLE gene_archive add COLUMN peptide_archive_id int NOT NULL AFTER translation_version;


Create TABLE new_peptide_archive(
  peptide_archive_id         INT NOT NULL AUTO_INCREMENT,
  md5_checksum             char(32),
  peptide_seq                mediumtext NOT NULL,

  translation_stable_id varchar(255),
  translation_version int,

  PRIMARY KEY( peptide_archive_id ),
  KEY checksum( md5_checksum )

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


INSERT INTO new_peptide_archive( md5_checksum, peptide_seq,
       translation_stable_id, translation_version)
SELECT UPPER(CAST( MD5( peptide_seq ) as CHAR CHARACTER SET LATIN1 )
       COLLATE latin1_swedish_ci),
       peptide_seq, translation_stable_id, translation_version
FROM   peptide_archive;

DROP TABLE peptide_archive;

UPDATE gene_archive ga, new_peptide_archive npa
SET ga.peptide_archive_id = npa.peptide_archive_id
WHERE ga.translation_stable_id = npa.translation_stable_id
AND ga.translation_version = npa.translation_version;

CREATE TABLE peptide_archive(
  peptide_archive_id         INT NOT NULL AUTO_INCREMENT,
  md5_checksum             char(32),
  peptide_seq                mediumtext NOT NULL,

  PRIMARY KEY( peptide_archive_id ),
  KEY checksum( md5_checksum )

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

INSERT INTO peptide_archive
SELECT peptide_archive_id, md5_checksum, peptide_seq
FROM new_peptide_archive;

DROP TABLE new_peptide_archive; 

# new gene_attrib table
CREATE TABLE gene_attrib (
  gene_id                     int(10) unsigned NOT NULL default '0',
  attrib_type_id              smallint(5) unsigned NOT NULL default '0',
  value                       varchar(255) NOT NULL default '',

  KEY type_val_idx( attrib_type_id, value ),
  KEY gene_idx( gene_id )
) COLLATE=latin1_swedish_ci TYPE=MyISAM;

