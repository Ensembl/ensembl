# Schema for internal-external database mappings (xrefs)


################################################################################
#
# General external annotation.

CREATE TABLE xref (

  xref_id                     int unsigned not null auto_increment,
  accession                   varchar(255) not null,
  label                       varchar(255),
  source_id                   int unsigned not null,
  species_id                  int unsigned not null,

  PRIMARY KEY (xref_id)

);

################################################################################

CREATE TABLE primary_xref (

  xref_id                     int unsigned not null,
  sequence                    mediumtext,
  sequence_type               enum('dna','peptide'),
  status                      enum('experimental','predicted'),
  source_id                   int unsigned not null,

  PRIMARY KEY (xref_id)

);

################################################################################

CREATE TABLE dependent_xref (

  master_xref_id              int unsigned not null,
  dependent_xref_id           int unsigned not null,
  linkage_annotation          varchar(255),
  source_id                   int unsigned not null,

  KEY master_idx(master_xref_id),
  KEY dependent_idx(dependent_xref_id)

);

################################################################################

CREATE TABLE synonym (

  xref_id                     int unsigned not null,
  synonym_xref_id             int unsigned not null,
  source_id                   int unsigned not null,

  KEY xref_idx(xref_id)

);

################################################################################

CREATE TABLE source (

  source_id                   int unsigned not null auto_increment,
  name                        varchar(255) not null,
  url                         varchar(255),
  checksum                    varchar(255),
  file_modified_date          datetime,
  upload_date                 datetime,
  release                     varchar(255),

  PRIMARY KEY (source_id),
  KEY name_idx(name) 

);

# Sources to fetch data from
# Note currently now UniProt/SwissProt data for fugu, anopheles, c.briggsae or chicken.
INSERT INTO source VALUES (1, 'UniProt_SwissProt_homo_sapiens', 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/9606.SPC', '', now(), now(), 1);
INSERT INTO source VALUES (2, 'UniProt_SwissProt_mus_musculus', 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/10090.SPC', '', now(), now(), 1);
INSERT INTO source VALUES (3, 'UniProt_SwissProt_rattus_norvegicus', 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/10116.SPC', '', now(), now(), 1);
INSERT INTO source VALUES (4, 'UniProt_SwissProt_drosophilla_melanogaster', 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/7227.SPC', '', now(), now(), 1);
INSERT INTO source VALUES (5, 'UniProt_SwissProt_caenorhabditis_elegans', 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/6239.SPC', '', now(), now(), 1);
INSERT INTO source VALUES (6, 'UniProt_SwissProt_gallus_gallus', 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/9031.SPC', '', now(), now(), 1);
INSERT INTO source VALUES (7, 'UniProt_SwissProt_pan_troglodytes', 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/9598.SPC', '', now(), now(), 1);

# Other sources - used to create dependent xrefs, but not to upload from
INSERT INTO source VALUES (100, 'EMBL', null, null, null, null, null);

################################################################################

CREATE TABLE species (

  species_id                  int unsigned not null auto_increment,
  taxonomy_id                 int unsigned not null,
  name                        varchar(255) not null,

  PRIMARY KEY(species_id),
  KEY taxonomy_idx(taxonomy_id),
  KEY name_idx(name)

);

INSERT INTO species (taxonomy_id, name) VALUES (9606,  'homo_sapiens');
INSERT INTO species (taxonomy_id, name) VALUES (10090, 'mus_musculus');
INSERT INTO species (taxonomy_id, name) VALUES (10116, 'rattus_norvegicus');
INSERT INTO species (taxonomy_id, name) VALUES (31033, 'fugu_rubripes');
INSERT INTO species (taxonomy_id, name) VALUES (7165,  'anopheles_gambiae');
INSERT INTO species (taxonomy_id, name) VALUES (7227,  'drosophila_melanogaster');
INSERT INTO species (taxonomy_id, name) VALUES (6239,  'caenorhabditis_elegans');
INSERT INTO species (taxonomy_id, name) VALUES (6238,  'caenorhabditis_briggsae');
INSERT INTO species (taxonomy_id, name) VALUES (7955,  'danio_rerio');
INSERT INTO species (taxonomy_id, name) VALUES (9598,  'pan_troglodytes');
INSERT INTO species (taxonomy_id, name) VALUES (9031,  'gallus_gallus');


################################################################################

CREATE TABLE direct_xref (

  general_xref_id             int unsigned not null,
  ensembl_stable_id           varchar(255),
  type                        enum('gene','transcript','translation'),
  linkage_xref                varchar(255),

  KEY primary_idx(general_xref_id),
  KEY ensembl_idx(ensembl_stable_id)

);

################################################################################

