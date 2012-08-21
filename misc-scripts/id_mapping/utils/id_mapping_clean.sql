#
# some useful stable ID mapping related SQL
# just fragments for copy/paste, not intended to be run as an "SQL script"
#

# get counts from all stable ID related tables
SELECT COUNT(*) FROM gene where stable_id IS NOT NULL;
SELECT COUNT(*) FROM transcript where stable_id IS NOT NULL
SELECT COUNT(*) FROM translation where stable_id IS NOT NULL;
SELECT COUNT(*) FROM exon where stable_id IS NOT NULL;
SELECT COUNT(*) FROM mapping_session;
SELECT COUNT(*) FROM stable_id_event;
SELECT COUNT(*) FROM gene_archive;
SELECT COUNT(*) FROM peptide_archive;

# backup all stable ID related tables
CREATE TABLE gene_bak LIKE gene;
INSERT INTO  gene_bak SELECT * FROM gene;
CREATE TABLE transcript_bak LIKE transcript;
INSERT INTO  transcript_bak SELECT * FROM transcript;
CREATE TABLE translation_bak LIKE translation;
INSERT INTO  translation_bak SELECT * FROM translation;
CREATE TABLE exon_bak LIKE exon;
INSERT INTO  exon_bak SELECT * FROM exon;
CREATE TABLE mapping_session_bak SELECT * FROM mapping_session;
CREATE TABLE stable_id_event_bak SELECT * FROM stable_id_event;
CREATE TABLE gene_archive_bak SELECT * FROM gene_archive;
CREATE TABLE peptide_archive_bak SELECT * FROM peptide_archive;

# now prune all of them
UPDATE gene set stable_id = NULL;
UPDATE transcript set stable_id = NULL;
UPDATE translation set stable_id = NULL;
UPDATE exon set stable_id = NULL;

# drop backup tables
DROP TABLE gene_bak;
DROP TABLE transcript_bak;
DROP TABLE translation_bak;
DROP TABLE exon_bak;
DROP TABLE mapping_session_bak;
DROP TABLE stable_id_event_bak;
DROP TABLE gene_archive_bak;
DROP TABLE peptide_archive_bak;

