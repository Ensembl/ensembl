# Schema for internal-external database mappings (xrefs)


################################################################################
#
# General external annotation.

CREATE TABLE xref (

  xref_id                     int unsigned not null auto_increment,
  accession                   varchar(255) not null,
  label                       varchar(255),
  description                 varchar(255),
  source_id                   int unsigned not null,
  species_id                  int unsigned not null,

  PRIMARY KEY (xref_id),
  UNIQUE acession_idx(accession,source_id)

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
  release                     varchar(255),
  download                    enum('Y', 'N') default 'Y',
   
  PRIMARY KEY (source_id),
  KEY name_idx(name) 

);

################################################################################

CREATE TABLE source_url (

  source_url_id               int unsigned not null auto_increment,
  source_id                   int unsigned not null,
  url                         varchar(255),
  checksum                    varchar(255),
  file_modified_date          datetime,
  upload_date                 datetime,
  parser                      varchar(255),

  PRIMARY KEY (source_url_id),
  KEY source_idx(source_id)

);

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

CREATE TABLE species (

  species_id                  int unsigned not null auto_increment,
  taxonomy_id                 int unsigned not null,
  name                        varchar(255) not null,

  PRIMARY KEY(species_id),
  KEY taxonomy_idx(taxonomy_id),
  KEY name_idx(name)

);

################################################################################

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
INSERT INTO species (taxonomy_id, name) VALUES (99883, 'tetraodon_nigroviridis');
INSERT INTO species (taxonomy_id, name) VALUES (170,   'apis_mellifera');

################################################################################

# "High level" sources that we will also download from (via source_url)

INSERT INTO source VALUES (1, "UniProtSwissProt", 1, 'Y');
INSERT INTO source VALUES (2, "RefSeq", 1, 'Y');

# Other sources - used to create dependent xrefs, but not to upload from

INSERT INTO source VALUES (1000, 'EMBL', 1, 'N');
INSERT INTO source VALUES (1001, 'PUBMED', 1, 'N');
INSERT INTO source VALUES (1002, 'MEDLINE', 1, 'N');

# Files to fetch data from

# --------------------------------------------------------------------------------
# UniProt/SwissProt
# Note currently no UniProt/SwissProt data for fugu, anopheles, c.briggsae or chicken.
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/9606.SPC', '', now(), now(), "SwissProtParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/10090.SPC', '', now(), now(), "SwissProtParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/10116.SPC', '', now(), now(), "SwissProtParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/7227.SPC', '', now(), now(), "SwissProtParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/6239.SPC', '', now(), now(), "SwissProtParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/9031.SPC', '', now(), now(), "SwissProtParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/9598.SPC', '', now(), now(), "SwissProtParser");

# --------------------------------------------------------------------------------
# RefSeq - release/ and cumulative/ directories, for protein and mRNA

# release/rna
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian1.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian2.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian3.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian4.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian5.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian6.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian7.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian8.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian9.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian10.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian11.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian12.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian13.rna.fna.gz', '', now(), now(), "RefSeqParser");

# release/protein
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian1.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian2.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian3.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian4.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian5.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian6.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian7.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian8.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian9.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian10.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian11.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian12.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian13.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");

# cumulative
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/cumulative/rscu.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/cumulative/rscu.fna.gz', '', now(), now(), "RefSeqParser");

################################################################################

