# Ensembl table definitions
# 
# Note that more information about each table can be found in
# ensembl/docs/schema_description/

# Conventions:
#  - use lower case and underscores
#  - internal ids are integers named tablename_id
#  - same name is given in foreign key relations


################################################################################
#
# Table structure for table 'oligo_feature'
#

CREATE TABLE oligo_feature (
  
  oligo_feature_id      INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id         INT(10) UNSIGNED NOT NULL,
  seq_region_start      INT(10) UNSIGNED  NOT NULL,
  seq_region_end        INT(10) UNSIGNED  NOT NULL,
  seq_region_strand     TINYINT NOT NULL,
  mismatches            TINYINT,
  oligo_probe_id        INT(10) UNSIGNED NOT NULL,
  analysis_id           SMALLINT UNSIGNED NOT NULL,

  PRIMARY KEY (oligo_feature_id),
  KEY seq_region_idx (seq_region_id, seq_region_start),
  KEY probe_idx (oligo_probe_id)
  
) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'oligo_probe'
#
# Note that the primary key contains both the probe ID and the array ID because
# it is often possible to get the same probe on different arrays 
# e.g. (older Affy arrays are often subsets of newer arrays).
# We give them the same oligo_probe_id so that we only have to store the 
# features once.

CREATE TABLE oligo_probe (
  
  oligo_probe_id      INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  oligo_array_id      INT(10) UNSIGNED NOT NULL,
  probeset            VARCHAR(40),
  name                VARCHAR(40),
  description         TEXT,
  length              SMALLINT NOT NULL,

  PRIMARY KEY (oligo_probe_id, oligo_array_id),
  KEY probeset_idx (probeset),
  KEY array_idx (oligo_array_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'oligo_array'
#

CREATE TABLE oligo_array (

  oligo_array_id      INT(10) UNSIGNED NOT NULL auto_increment,
  parent_array_id     INT(10) UNSIGNED,
  probe_setsize       TINYINT NOT NULL,
  name                VARCHAR(40) NOT NULL,
  type                ENUM( 'AFFY', 'OLIGO' ),

  PRIMARY KEY (oligo_array_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'alt_allele'
#

CREATE TABLE alt_allele (
  alt_allele_id         INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  gene_id               INT(10) UNSIGNED NOT NULL,

  UNIQUE gene_idx (gene_id),
  UNIQUE allele_idx (alt_allele_id, gene_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'analysis'
# 
# semantics:
#
# analysis_id - internal id
# created 
#   - date to distinguish newer and older versions off the same analysis. Not
#     well maintained so far.
# logic_name - string to identify the analysis. Used mainly inside pipeline.
# db, db_version, db_file
#  - db should be a database name, db version the version of that db
#    db_file the file system location of that database, 
#    probably wiser to generate from just db and configurations
# program, program_version,program_file
#  - The binary used to create a feature. Similar semantic to above
# module, module_version
#  - Perl module names (RunnableDBS usually) executing this analysis.
# parameters - a paramter string which is processed by the perl module
# gff_source, gff_feature 
#  - how to make a gff dump from features with this analysis


CREATE TABLE analysis (

  analysis_id                 SMALLINT UNSIGNED NOT NULL AUTO_INCREMENT,
  created                     datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  logic_name                  VARCHAR(128) NOT NULL,
  db                          VARCHAR(120),
  db_version                  VARCHAR(40),
  db_file                     VARCHAR(120),
  program                     VARCHAR(80),
  program_version             VARCHAR(40),
  program_file                VARCHAR(80),
  parameters                  TEXT,
  module                      VARCHAR(80),
  module_version              VARCHAR(40),
  gff_source                  VARCHAR(40),
  gff_feature                 VARCHAR(40),

  PRIMARY KEY (analysis_id),
  KEY logic_name_idx (logic_name),
  UNIQUE (logic_name)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'analysis_description'
#

CREATE TABLE analysis_description (

  analysis_id	               SMALLINT UNSIGNED NOT NULL,
  description                  TEXT,
  display_label                VARCHAR(255) NOT NULL,
  displayable                  BOOLEAN NOT NULL DEFAULT 1,
  web_data                     TEXT,

  UNIQUE KEY analysis_idx (analysis_id)
  
) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'dna'
#
# This table stores DNA sequence.

CREATE TABLE dna (

  seq_region_id       INT(10) UNSIGNED NOT NULL,
  sequence            LONGTEXT NOT NULL,

  PRIMARY KEY (seq_region_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM MAX_ROWS=750000 AVG_ROW_LENGTH=19000;


################################################################################
#
# Table structure for table 'dnac'
#
# Contains equivalent data to dna table, but 4 letters of DNA code are
# represented by a single binary character, based on 2 bit encoding
#
# do not need to worry about ambiguity of length, since this is stored in
# contig.length
#
# n_line column contains start-end pairs of coordinates in the string that are
# really Ns

CREATE TABLE dnac (

  seq_region_id     INT(10) UNSIGNED NOT NULL,
  sequence          MEDIUMBLOB NOT NULL,
  n_line            TEXT,  

  PRIMARY KEY (seq_region_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM MAX_ROWS=750000 AVG_ROW_LENGTH=19000;


################################################################################
#
# Table structure for table 'exon'
#
# Note seq_region_start always less that seq_region_end, i.e. when the exon is
# on the other strand the seq_region_start is specifying the 3prime end of the
# exon.

CREATE TABLE exon (
 
  exon_id       	      INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id     	      INT(10) UNSIGNED NOT NULL,
  seq_region_start  	      INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  seq_region_strand           TINYINT(2) NOT NULL,

  phase                       TINYINT(2) NOT NULL,
  end_phase                   TINYINT(2) NOT NULL,

  is_current                  BOOLEAN NOT NULL DEFAULT 1,
  is_constitutive             BOOLEAN NOT NULL DEFAULT 0,
  
  PRIMARY KEY (exon_id),
  KEY seq_region_idx (seq_region_id, seq_region_start)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'exon_stable_id'
#

CREATE TABLE exon_stable_id (
   
  exon_id   		      INT(10) UNSIGNED NOT NULL,
  stable_id                   VARCHAR(128) NOT NULL,
  version                     INT(10),
  created_date                DATETIME NOT NULL,
  modified_date               DATETIME NOT NULL,

  PRIMARY KEY (exon_id),
  KEY stable_id_idx (stable_id, version)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'exon_transcript'
#
# Note that the rank column indicates the 5prime to 3prime position of the 
# exon within the transcript ie rank of 1 means the exon is the 5prime most 
# within this transcript

CREATE TABLE exon_transcript (

  exon_id                     INT(10) UNSIGNED NOT NULL,
  transcript_id               INT(10) UNSIGNED NOT NULL,
  rank                        INT(10) NOT NULL,        

  PRIMARY KEY (exon_id,transcript_id,rank),
  KEY transcript (transcript_id),
  KEY exon (exon_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'simple_feature'
#

CREATE TABLE simple_feature (

  simple_feature_id 	      INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  seq_region_strand           TINYINT(1) NOT NULL,
  display_label               VARCHAR(40) NOT NULL,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  score                       DOUBLE,

  PRIMARY KEY (simple_feature_id),
  KEY seq_region_idx (seq_region_id, seq_region_start),
  KEY analysis_idx (analysis_id),
  KEY hit_idx (display_label)

) COLLATE=latin1_swedish_ci TYPE=MyISAM MAX_ROWS=100000000 AVG_ROW_LENGTH=80;


################################################################################
#
# Table structure for table 'protein_align_feature'
#

CREATE TABLE protein_align_feature (
 
  protein_align_feature_id    INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  seq_region_strand           TINYINT(1) DEFAULT '1' NOT NULL,
  hit_start                   INT(10) NOT NULL,
  hit_end                     INT(10) NOT NULL,
  hit_name                    VARCHAR(40) NOT NULL,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  score                       DOUBLE,
  evalue                      DOUBLE,
  perc_ident                  FLOAT,
  cigar_line                  TEXT,
  external_db_id              SMALLINT UNSIGNED,
  hcoverage                   DOUBLE,

  PRIMARY KEY (protein_align_feature_id),
  KEY seq_region_idx (seq_region_id, analysis_id, seq_region_start, score),
  KEY seq_region_idx_2 (seq_region_id, seq_region_start),
  KEY hit_idx (hit_name),
  KEY analysis_idx (analysis_id),
  KEY external_db_idx (external_db_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM MAX_ROWS=100000000 AVG_ROW_LENGTH=80;


################################################################################
#
# Table structure for table 'dna_align_feature'
#

CREATE TABLE dna_align_feature (

  dna_align_feature_id 	      INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  seq_region_strand           TINYINT(1) NOT NULL,
  hit_start                   INT NOT NULL,
  hit_end                     INT NOT NULL,
  hit_strand                  TINYINT(1) NOT NULL,
  hit_name                    VARCHAR(40) NOT NULL,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  score                       DOUBLE,
  evalue                      DOUBLE,
  perc_ident                  FLOAT,
  cigar_line                  TEXT,
  external_db_id              SMALLINT UNSIGNED,
  hcoverage                   DOUBLE,
  external_data               TEXT, 
  pair_dna_align_feature_id   INT(10) UNSIGNED,

  PRIMARY KEY (dna_align_feature_id),
  KEY seq_region_idx (seq_region_id, analysis_id, seq_region_start, score),
  KEY seq_region_idx_2 (seq_region_id, seq_region_start),
  KEY hit_idx (hit_name),
  KEY analysis_idx (analysis_id),
  KEY external_db_idx (external_db_id),
  KEY pair_idx (pair_dna_align_feature_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM MAX_ROWS=100000000 AVG_ROW_LENGTH=80;


################################################################################
#
# Table structure for table 'repeat_consensus'
#
# repeat_class examples: SINE, LINE, DNA Transposon, Retroviral LTR,
# Satellite, Tandem

CREATE TABLE repeat_consensus (

  repeat_consensus_id  	      INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  repeat_name                 VARCHAR(255) NOT NULL,
  repeat_class                VARCHAR(100) NOT NULL, 
  repeat_type                 VARCHAR(40) NOT NULL,
  repeat_consensus            TEXT,
  
  PRIMARY KEY (repeat_consensus_id),
  KEY name (repeat_name),
  KEY class (repeat_class),
  KEY consensus (repeat_consensus(10)),
  KEY type (repeat_type)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

################################################################################
#
# Table structure for table 'repeat_feature'
#

CREATE TABLE repeat_feature (

  repeat_feature_id 	      INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  seq_region_strand           TINYINT(1) DEFAULT '1' NOT NULL,
  repeat_start                INT(10) NOT NULL,
  repeat_end                  INT(10) NOT NULL,
  repeat_consensus_id         INT(10) UNSIGNED NOT NULL,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  score                       DOUBLE,
  
  PRIMARY KEY (repeat_feature_id),
  KEY seq_region_idx (seq_region_id, seq_region_start),
  KEY repeat_idx (repeat_consensus_id),
  KEY analysis_idx (analysis_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM MAX_ROWS=100000000 AVG_ROW_LENGTH=80;


################################################################################
#
# Table structure for table 'gene'
#

CREATE TABLE gene (

  gene_id                     INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  biotype                     VARCHAR(40) NOT NULL,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  seq_region_id               INT(10) UNSIGNED NOT NULL, 
  seq_region_start            INT(10) UNSIGNED NOT NULL, 
  seq_region_end              INT(10) UNSIGNED NOT NULL, 
  seq_region_strand           TINYINT(2) NOT NULL,       
  display_xref_id             INT(10) UNSIGNED,
  source                      VARCHAR(20) NOT NULL,
  status                      ENUM('KNOWN', 'NOVEL', 'PUTATIVE', 'PREDICTED', 'KNOWN_BY_PROJECTION', 'UNKNOWN'),
  description                 TEXT,
  is_current                  BOOLEAN NOT NULL DEFAULT 1,
  canonical_transcript_id     INT(10) UNSIGNED NOT NULL,
  canonical_annotation        VARCHAR(255) DEFAULT NULL,

  PRIMARY KEY (gene_id),
  KEY seq_region_idx (seq_region_id, seq_region_start),
  KEY xref_id_index (display_xref_id),
  KEY analysis_idx (analysis_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'gene_stable_id'
#

CREATE TABLE gene_stable_id (

  gene_id 		      INT UNSIGNED NOT NULL,
  stable_id                   VARCHAR(128) NOT NULL,
  version                     INT(10),
  created_date                DATETIME NOT NULL,
  modified_date               DATETIME NOT NULL,

  PRIMARY KEY (gene_id),
  KEY stable_id_idx (stable_id, version)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'supporting_feature'
#

CREATE TABLE supporting_feature (

  exon_id 		      INT(10) UNSIGNED DEFAULT '0' NOT NULL,
  feature_type                ENUM('dna_align_feature','protein_align_feature'),
  feature_id                  INT(10) UNSIGNED DEFAULT '0' NOT NULL,

  UNIQUE all_idx (exon_id,feature_type,feature_id),
  KEY feature_idx (feature_type,feature_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM MAX_ROWS=100000000 AVG_ROW_LENGTH=80;


################################################################################
#
# Table structure for table 'transcript_supporting_feature'
#

CREATE TABLE transcript_supporting_feature (

  transcript_id 	      INT(10) UNSIGNED DEFAULT '0' NOT NULL,
  feature_type                ENUM('dna_align_feature','protein_align_feature'),
  feature_id                  INT(10) UNSIGNED DEFAULT '0' NOT NULL,

  UNIQUE all_idx (transcript_id,feature_type,feature_id),
  KEY feature_idx (feature_type,feature_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM MAX_ROWS=100000000 AVG_ROW_LENGTH=80;


################################################################################
#
# Table structure for table 'transcript'
#

CREATE TABLE transcript (

  transcript_id               INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,  
  gene_id                     INT(10) UNSIGNED,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  seq_region_id               INT(10) UNSIGNED NOT NULL, 
  seq_region_start            INT(10) UNSIGNED NOT NULL, 
  seq_region_end              INT(10) UNSIGNED NOT NULL, 
  seq_region_strand           TINYINT(2) NOT NULL, 
  display_xref_id             INT(10) UNSIGNED,
  biotype                     VARCHAR(40) NOT NULL,
  status                      ENUM('KNOWN', 'NOVEL', 'PUTATIVE', 'PREDICTED', 'KNOWN_BY_PROJECTION', 'UNKNOWN'),
  description                 TEXT,
  is_current                  BOOLEAN NOT NULL DEFAULT 1,

  PRIMARY KEY (transcript_id),
  KEY seq_region_idx (seq_region_id, seq_region_start),
  KEY gene_index (gene_id),
  KEY xref_id_index (display_xref_id),
  KEY analysis_idx (analysis_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'transcript_stable_id'
#

CREATE TABLE transcript_stable_id (

  transcript_id               INT(10) UNSIGNED NOT NULL,
  stable_id                   VARCHAR(128) NOT NULL,
  version                     INT(10),
  created_date                DATETIME NOT NULL,
  modified_date               DATETIME NOT NULL,
  
  PRIMARY KEY (transcript_id),
  KEY stable_id_idx (stable_id, version)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'translation'
#
# The seq_start and seq_end are 1-based offsets into the *relative* coordinate 
# system of start_exon_id and end_exon_id. i.e. if the translation starts at the
# first base of the exon, seq_start would be 1

CREATE TABLE translation (

  translation_id              INT(10) UNSIGNED NOT NULL AUTO_INCREMENT, 
  transcript_id               INT(10) UNSIGNED NOT NULL, 
  seq_start                   INT(10) NOT NULL,       # relative to exon start
  start_exon_id               INT(10) UNSIGNED NOT NULL,
  seq_end                     INT(10) NOT NULL,       # relative to exon start
  end_exon_id                 INT(10) UNSIGNED NOT NULL,
  
  PRIMARY KEY (translation_id),
  KEY (transcript_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'translation_stable_id'
#
CREATE TABLE translation_stable_id (

  translation_id 	      INT(10) UNSIGNED NOT NULL,
  stable_id                   VARCHAR(128) NOT NULL,
  version                     INT(10),
  created_date                DATETIME NOT NULL,
  modified_date               DATETIME NOT NULL,

  PRIMARY KEY (translation_id),
  KEY stable_id_idx (stable_id, version)	
 
) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'assembly'
#
# This is a denormalised golden path. 
#
# The data in this table defines the "static golden path", i.e. the best effort
# draft full genome sequence as determined by the UCSC or NCBI (depending which
# assembly you are using).
#
# Each row represents a component, e.g. a contig,  (comp_seq_region_id, FK from
# seq_region table) at least part of which is present in the golden path. 
#
# The part of the component that is in the path is delimited by fields
# cmp_start and cmp_end (start < end), and the absolute position within the
# golden path chromosome (or other appropriate assembled structure)
# (asm_seq_region_id) is given by asm_start and asm_end. 

CREATE TABLE assembly (

  asm_seq_region_id           INT(10) UNSIGNED NOT NULL,
  cmp_seq_region_id           INT(10) UNSIGNED NOT NULL, 
  asm_start                   INT(10) NOT NULL,
  asm_end                     INT(10) NOT NULL,
  cmp_start                   INT(10) NOT NULL,
  cmp_end                     INT(10) NOT NULL,
  ori                         TINYINT  NOT NULL, 
  
  KEY (cmp_seq_region_id),
  KEY (asm_seq_region_id, asm_start),
  UNIQUE KEY all_idx (asm_seq_region_id, cmp_seq_region_id, asm_start, asm_end, cmp_start, cmp_end, ori)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'protein_feature'
#

CREATE TABLE protein_feature (

  protein_feature_id          INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  translation_id              INT(10) UNSIGNED NOT NULL,	
  seq_start                   INT(10) NOT NULL,
  seq_end                     INT(10) NOT NULL,
  hit_start                   INT(10) NOT NULL,
  hit_end                     INT(10) NOT NULL,
  hit_name                    VARCHAR(40) NOT NULL,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  score                       DOUBLE,
  evalue                      DOUBLE,
  perc_ident                  FLOAT,
  external_data               TEXT,

  PRIMARY KEY (protein_feature_id),
  KEY (translation_id),
  KEY hitname_idx (hit_name),
  KEY analysis_idx (analysis_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'interpro'
#

CREATE TABLE interpro (

  interpro_ac	              VARCHAR(40) NOT NULL,
  id		              VARCHAR(40) NOT NULL,

  UNIQUE (interpro_ac, id),
  KEY (id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'karyotype'
#

CREATE TABLE karyotype (
  karyotype_id                INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10)     NOT NULL,
  seq_region_end              INT(10)     NOT NULL,
  band                        VARCHAR(40) NOT NULL,
  stain                       VARCHAR(40) NOT NULL,

  PRIMARY KEY (karyotype_id),
  KEY region_band_idx (seq_region_id,band)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'object_xref'
#

CREATE TABLE object_xref (

  object_xref_id              INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  ensembl_id                  INT(10) UNSIGNED NOT NULL, 
  ensembl_object_type         ENUM('RawContig', 'Transcript', 'Gene',
                                   'Translation')
                              NOT NULL,
  xref_id                     INT UNSIGNED NOT NULL,
  linkage_annotation          VARCHAR(255) DEFAULT NULL,
  analysis_id                 SMALLINT UNSIGNED DEFAULT NULL,

  UNIQUE (ensembl_object_type, ensembl_id, xref_id),
  KEY oxref_idx (object_xref_id, xref_id, ensembl_object_type, ensembl_id),
  KEY xref_idx (xref_id, ensembl_object_type),
  KEY analysis_idx (analysis_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'identity_xref'
#

CREATE TABLE identity_xref (

  object_xref_id          INT(10) UNSIGNED NOT NULL,
  xref_identity 	  INT(5),
  ensembl_identity        INT(5),

  xref_start              INT,
  xref_end                INT,
  ensembl_start           INT,
  ensembl_end             INT,
  cigar_line              TEXT,
  
  score                   DOUBLE,
  evalue                  DOUBLE,

  PRIMARY KEY (object_xref_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'go_xref'
#

CREATE TABLE go_xref (

  object_xref_id          INT(10) UNSIGNED DEFAULT '0' NOT NULL,
  linkage_type            ENUM('IC', 'IDA', 'IEA', 'IEP', 'IGI', 'IMP', 
		               'IPI', 'ISS', 'NAS', 'ND', 'TAS', 'NR', 'RCA',
			       'EXP','ISO','ISA','ISM','IGC'),
  source_xref_id          INT(10) UNSIGNED DEFAULT NULL,
  KEY (source_xref_id),
  UNIQUE (object_xref_id, source_xref_id, linkage_type)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'xref'
#

CREATE TABLE xref (

   xref_id 		      INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
   external_db_id             SMALLINT UNSIGNED NOT NULL,
   dbprimary_acc              VARCHAR(40) NOT NULL,
   display_label              VARCHAR(128) NOT NULL,
   version                    VARCHAR(10) DEFAULT '0' NOT NULL,
   description                TEXT,
   info_type                  ENUM( 'PROJECTION', 'MISC', 'DEPENDENT',
                                    'DIRECT', 'SEQUENCE_MATCH',
                                    'INFERRED_PAIR', 'PROBE',
                                    'UNMAPPED', 'COORDINATE_OVERLAP' ),
   info_text                  VARCHAR(255),

   PRIMARY KEY (xref_id),
   UNIQUE KEY id_index (dbprimary_acc, external_db_id, info_type, info_text),
   KEY display_index (display_label),
   KEY info_type_idx (info_type)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'dependent_xref'
#

CREATE TABLE dependent_xref(

  object_xref_id         INT NOT NULL,
  master_xref_id         INT NOT NULL,
  dependent_xref_id      INT NOT NULL,

  PRIMARY KEY( object_xref_id ),
  KEY dependent ( dependent_xref_id ),
  KEY master_idx (master_xref_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

################################################################################
#
# Table structure for table 'external_synonym'
#

CREATE TABLE external_synonym (

  xref_id                     INT(10) UNSIGNED NOT NULL,
  synonym                     VARCHAR(40) NOT NULL,
  
  PRIMARY KEY (xref_id, synonym),
  KEY name_index (synonym)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'external_db' 
#

CREATE TABLE external_db (

  external_db_id 	      SMALLINT UNSIGNED NOT NULL,
  db_name                     VARCHAR(100) NOT NULL,
  db_release                  VARCHAR(255),
  status                      ENUM('KNOWNXREF','KNOWN','XREF','PRED','ORTH',
                                   'PSEUDO')
                              NOT NULL,
  dbprimary_acc_linkable      BOOLEAN DEFAULT 1 NOT NULL,
  display_label_linkable      BOOLEAN DEFAULT 0 NOT NULL,
  priority                    INT NOT NULL,
  db_display_name             VARCHAR(255),
  type                        ENUM('ARRAY', 'ALT_TRANS', 'MISC', 'LIT', 'PRIMARY_DB_SYNONYM', 'ENSEMBL'),
  secondary_db_name           VARCHAR(255) DEFAULT NULL,
  secondary_db_table          VARCHAR(255) DEFAULT NULL,
  description                 TEXT,
	
  PRIMARY KEY (external_db_id) 

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'prediction_exon' 
#

CREATE TABLE prediction_exon (

  prediction_exon_id          INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  prediction_transcript_id    INT(10) UNSIGNED NOT NULL,
  exon_rank                   SMALLINT UNSIGNED NOT NULL,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  seq_region_strand           TINYINT NOT NULL,
  start_phase                 TINYINT NOT NULL,
  score                       DOUBLE,
  p_value                     DOUBLE,

  PRIMARY KEY (prediction_exon_id),
  KEY (prediction_transcript_id),
  KEY (seq_region_id, seq_region_start)
  
) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'prediction_transcript' 
#

CREATE TABLE prediction_transcript (

  prediction_transcript_id    INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  seq_region_strand           TINYINT NOT NULL,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  display_label               VARCHAR(255),
  
  PRIMARY KEY (prediction_transcript_id),
  KEY (seq_region_id, seq_region_start),
  KEY analysis_idx (analysis_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'meta' 
#

CREATE TABLE meta (

  meta_id                     INT NOT NULL AUTO_INCREMENT,
  species_id                  INT UNSIGNED DEFAULT 1,
  meta_key                    VARCHAR(40) NOT NULL,
  meta_value                  VARCHAR(255) BINARY NOT NULL,

  PRIMARY   KEY (meta_id),
  UNIQUE    KEY species_key_value_idx (species_id, meta_key, meta_value),
            KEY species_value_idx (species_id, meta_value)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


# Auto add schema version to database
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, "schema_version", "55");

# patches included in this schema file
# NOTE: at beginning of release cycle, remove patch entries from last release
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_a.sql|schema_version');
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_b.sql|add_go_xrefs_types');
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_c.sql|add_splicing_event_tables');
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_d.sql|add_dependent_xref_table');
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_e.sql|add_is_constitutive_column');
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_f.sql|coord_system.version_null');
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_g.sql|analysis_description.display_label_NOT_NULL');
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_h.sql|gene_archive.allow_for_NULLs');

################################################################################
#
# Table structure for table 'marker_synonym'

CREATE TABLE marker_synonym (

  marker_synonym_id           INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  marker_id                   INT(10) UNSIGNED NOT NULL,
  source                      VARCHAR(20),
  name                        VARCHAR(50),    

  PRIMARY KEY (marker_synonym_id),
  KEY marker_synonym_idx (marker_synonym_id, name),
  KEY marker_idx (marker_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'marker'

CREATE TABLE marker (

  marker_id                   INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  display_marker_synonym_id   INT(10) UNSIGNED,
  left_primer                 VARCHAR(100) NOT NULL,
  right_primer                VARCHAR(100) NOT NULL,
  min_primer_dist             INT(10) UNSIGNED NOT NULL,
  max_primer_dist             INT(10) UNSIGNED NOT NULL,
  priority                    INT,
  type                        ENUM('est', 'microsatellite'),
  
  PRIMARY KEY (marker_id),
  KEY marker_idx (marker_id, priority),
  KEY display_idx (display_marker_synonym_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'marker_feature'

CREATE TABLE marker_feature (

  marker_feature_id           INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  marker_id                   INT(10) UNSIGNED NOT NULL,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  map_weight                  INT(10) UNSIGNED,

  PRIMARY KEY (marker_feature_id),
  KEY seq_region_idx (seq_region_id, seq_region_start),
  KEY analysis_idx (analysis_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;
   

################################################################################
#
# Table structure for table 'marker_map_location'
 
CREATE TABLE marker_map_location (

  marker_id                   INT(10) UNSIGNED NOT NULL,
  map_id                      INT(10) UNSIGNED NOT NULL,
  chromosome_name             VARCHAR(15)  NOT NULL, 
  marker_synonym_id           INT(10) UNSIGNED NOT NULL,
  position                    VARCHAR(15) NOT NULL,
  lod_score                   DOUBLE,
  
  PRIMARY KEY (marker_id, map_id),
  KEY map_idx (map_id, chromosome_name, position) 

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'map'

CREATE TABLE map (

  map_id                      INT(10) UNSIGNED NOT NULL AUTO_INCREMENT, 
  map_name                    VARCHAR(30) NOT NULL,

  PRIMARY KEY (map_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

    
################################################################################
#
# Table structure for table 'misc_feature'
#

CREATE TABLE misc_feature (

  misc_feature_id 	      INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id               INT(10) UNSIGNED NOT NULL DEFAULT '0',
  seq_region_start            INT(10) UNSIGNED NOT NULL DEFAULT '0',
  seq_region_end              INT(10) UNSIGNED NOT NULL DEFAULT '0',
  seq_region_strand           TINYINT(4) NOT NULL DEFAULT '0',

  PRIMARY KEY (misc_feature_id),
  KEY seq_region_idx (seq_region_id, seq_region_start)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'misc_attrib'
#

CREATE TABLE misc_attrib (

  misc_feature_id             INT(10) UNSIGNED NOT NULL DEFAULT '0',
  attrib_type_id              SMALLINT(5) UNSIGNED NOT NULL DEFAULT '0',
  value                       TEXT NOT NULL DEFAULT '',

  KEY type_val_idx (attrib_type_id, value(40)),
  KEY val_only_idx (value(40)),
  KEY misc_feature_idx (misc_feature_id)
  
) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'translation_attrib'
#

CREATE TABLE translation_attrib (

  translation_id              INT(10) UNSIGNED NOT NULL DEFAULT '0',
  attrib_type_id              SMALLINT(5) UNSIGNED NOT NULL DEFAULT '0',
  value                       TEXT NOT NULL DEFAULT '',

  KEY type_val_idx (attrib_type_id, value(40)),
  KEY val_only_idx (value(40)),
  KEY translation_idx (translation_id)
  
) COLLATE=latin1_swedish_ci TYPE=MyISAM;

################################################################################
#
# Table structure for table 'transcript_attrib'
#

CREATE TABLE transcript_attrib (

  transcript_id               INT(10) UNSIGNED NOT NULL DEFAULT '0',
  attrib_type_id              SMALLINT(5) UNSIGNED NOT NULL DEFAULT '0',
  value                       TEXT NOT NULL DEFAULT '',

  KEY type_val_idx (attrib_type_id, value(40)),
  KEY val_only_idx (value(40)),
  KEY transcript_idx (transcript_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'gene_attrib'
#

CREATE TABLE gene_attrib (

  gene_id                     INT(10) UNSIGNED NOT NULL DEFAULT '0',
  attrib_type_id              SMALLINT(5) UNSIGNED NOT NULL DEFAULT '0',
  value                       TEXT NOT NULL DEFAULT '',

  KEY type_val_idx (attrib_type_id, value(40)),
  KEY val_only_idx (value(40)),
  KEY gene_idx (gene_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'seq_region_attrib'
#

CREATE TABLE seq_region_attrib (

  seq_region_id               INT(10) UNSIGNED NOT NULL DEFAULT '0',
  attrib_type_id              SMALLINT(5) UNSIGNED NOT NULL DEFAULT '0',
  value                       TEXT NOT NULL DEFAULT '',

  KEY type_val_idx (attrib_type_id, value(40)),
  KEY val_only_idx (value(40)),
  KEY seq_region_idx (seq_region_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'attrib_type'
#

CREATE TABLE attrib_type (

  attrib_type_id              SMALLINT(5) UNSIGNED NOT NULL AUTO_INCREMENT,
  code                        VARCHAR(15) NOT NULL DEFAULT '',
  name                        VARCHAR(255) NOT NULL DEFAULT '',
  description                 TEXT,

  PRIMARY KEY (attrib_type_id),
  UNIQUE KEY c (code)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'misc_set'
#

CREATE TABLE misc_set (

  misc_set_id                 SMALLINT(5) UNSIGNED NOT NULL AUTO_INCREMENT,
  code                        VARCHAR(25) NOT NULL DEFAULT '',
  name                        VARCHAR(255) NOT NULL DEFAULT '',
  description                 TEXT NOT NULL,
  max_length                  INT UNSIGNED NOT NULL,

  PRIMARY KEY (misc_set_id),
  UNIQUE KEY c (code)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'misc_feature_misc_set'
#

CREATE TABLE misc_feature_misc_set (

  misc_feature_id 	      INT(10) UNSIGNED NOT NULL DEFAULT '0',
  misc_set_id 		      SMALLINT(5) UNSIGNED NOT NULL DEFAULT '0',

  PRIMARY KEY (misc_feature_id, misc_set_id),
  KEY reverse_idx (misc_set_id, misc_feature_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'qtl'
#

CREATE TABLE qtl (

  qtl_id 		      INT(10) UNSIGNED AUTO_INCREMENT NOT NULL,
  trait                       VARCHAR(255) NOT NULL,
  lod_score                   FLOAT,
  flank_marker_id_1 	      INT(10) UNSIGNED,
  flank_marker_id_2           INT(10) UNSIGNED,
  peak_marker_id              INT(10) UNSIGNED,

  PRIMARY KEY (qtl_id),
  KEY trait_idx (trait)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'qtl_synonym'
#

CREATE TABLE qtl_synonym (

  qtl_synonym_id 	      INT(10) UNSIGNED AUTO_INCREMENT NOT NULL,
  qtl_id                      INT(10) UNSIGNED NOT NULL,
  source_database             ENUM("rat genome database", "ratmap") NOT NULL,
  source_primary_id           VARCHAR(255) NOT NULL,

  PRIMARY KEY (qtl_synonym_id),
  KEY qtl_idx (qtl_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'qtl_feature'
#

CREATE TABLE qtl_feature (

  seq_region_id 	INT(10)	UNSIGNED NOT NULL,
  seq_region_start      INT(10)	UNSIGNED NOT NULL,
  seq_region_end        INT(10)	UNSIGNED NOT NULL,
  qtl_id                INT(10)	UNSIGNED NOT NULL,
  analysis_id           SMALLINT UNSIGNED NOT NULL,

  KEY (qtl_id),
  KEY loc_idx (seq_region_id, seq_region_start),
  KEY analysis_idx (analysis_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'mapping_session'
#

CREATE TABLE mapping_session (

  mapping_session_id 	      INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  old_db_name                 VARCHAR(80) NOT NULL DEFAULT '',
  new_db_name                 VARCHAR(80) NOT NULL DEFAULT '',
  old_release                 VARCHAR(5) NOT NULL DEFAULT '',
  new_release                 VARCHAR(5) NOT NULL DEFAULT '',
  old_assembly                VARCHAR(20) NOT NULL DEFAULT '',
  new_assembly                VARCHAR(20) NOT NULL DEFAULT '',
  created                     DATETIME NOT NULL,

  PRIMARY KEY (mapping_session_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'stable_id_event'
#

CREATE TABLE stable_id_event (

  old_stable_id 	    VARCHAR(128),
  old_version               SMALLINT,
  new_stable_id             VARCHAR(128),
  new_version               SMALLINT,
  mapping_session_id        INT(10) UNSIGNED NOT NULL DEFAULT '0',
  type                      ENUM('gene', 'transcript', 'translation') NOT NULL,
  score			    FLOAT NOT NULL DEFAULT 0,

  UNIQUE KEY uni_idx (mapping_session_id, old_stable_id, new_stable_id, type),

  KEY new_idx (new_stable_id),
  KEY old_idx (old_stable_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'gene_archive'
#

CREATE TABLE gene_archive (

  gene_stable_id	      VARCHAR(128) NOT NULL,
  gene_version                SMALLINT NOT NULL,
  transcript_stable_id        VARCHAR(128) NOT NULL,
  transcript_version          SMALLINT NOT NULL,
  translation_stable_id       VARCHAR(128),
  translation_version         SMALLINT,
  peptide_archive_id          INT(10) UNSIGNED,
  mapping_session_id          INT(10) UNSIGNED NOT NULL,

  KEY gene_idx (gene_stable_id, gene_version),
  KEY transcript_idx (transcript_stable_id, transcript_version),
  KEY translation_idx (translation_stable_id, translation_version),
  KEY peptide_archive_id_idx (peptide_archive_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'peptide_archive'
#

CREATE TABLE peptide_archive (

  peptide_archive_id         INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  md5_checksum               VARCHAR(32),
  peptide_seq                MEDIUMTEXT NOT NULL,

  PRIMARY KEY (peptide_archive_id),
  KEY checksum (md5_checksum)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'seq_region'
#

CREATE TABLE seq_region (

  seq_region_id               INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  name                        VARCHAR(40) NOT NULL,
  coord_system_id             INT(10) UNSIGNED NOT NULL,
  length                      INT(10) NOT NULL,

  PRIMARY KEY (seq_region_id),
  UNIQUE KEY name_cs_idx (name, coord_system_id),
  KEY cs_idx (coord_system_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'assembly_exception'
#

CREATE TABLE assembly_exception (

  assembly_exception_id       INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL, 
  exc_type                    ENUM('HAP', 'PAR') NOT NULL,
  exc_seq_region_id           INT(10) UNSIGNED NOT NULL, 
  exc_seq_region_start        INT(10) UNSIGNED NOT NULL, 
  exc_seq_region_end          INT(10) UNSIGNED NOT NULL,
  ori                         INT NOT NULL,

  PRIMARY KEY (assembly_exception_id),
  KEY sr_idx (seq_region_id, seq_region_start),
  KEY ex_idx (exc_seq_region_id, exc_seq_region_start)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'coord_system'
#

CREATE TABLE coord_system (

  coord_system_id             INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  species_id                  INT(10) UNSIGNED NOT NULL DEFAULT 1,
  name                        VARCHAR(40) NOT NULL,
  version                     VARCHAR(255) DEFAULT NULL,
  rank                        INT NOT NULL,
  attrib                      SET('default_version', 'sequence_level'),

  PRIMARY   KEY (coord_system_id),
  UNIQUE    KEY rank_idx (rank, species_id),
  UNIQUE    KEY name_idx (name, version, species_id),
            KEY species_idx (species_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'meta_coord'
#

CREATE TABLE meta_coord (

  table_name                  VARCHAR(40) NOT NULL,
  coord_system_id             INT(10) UNSIGNED NOT NULL,
  max_length                  INT,

  UNIQUE KEY cs_table_name_idx (coord_system_id, table_name)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'density_feature'
#

CREATE TABLE density_feature (

  density_feature_id    INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  density_type_id       INT(10) UNSIGNED NOT NULL,
  seq_region_id         INT(10) UNSIGNED NOT NULL,
  seq_region_start      INT(10) UNSIGNED NOT NULL,
  seq_region_end        INT(10) UNSIGNED NOT NULL,
  density_value         FLOAT NOT NULL,

  PRIMARY KEY (density_feature_id),
  KEY seq_region_idx (density_type_id, seq_region_id, seq_region_start),
  KEY seq_region_id_idx (seq_region_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'density_type'
#

CREATE TABLE density_type (

  density_type_id       INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  analysis_id           SMALLINT UNSIGNED NOT NULL,
  block_size            INT NOT NULL,
  region_features       INT NOT NULL,
  value_type            ENUM('sum','ratio') NOT NULL,
  
  PRIMARY KEY (density_type_id),
  UNIQUE (analysis_id, block_size, region_features)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;



################################################################################
#
# Table structure for table 'unmapped_object'
#
# Describes why a particular external entity was not mapped to an ensembl one.

CREATE TABLE unmapped_object (

  unmapped_object_id    INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  type                  ENUM('xref', 'cDNA', 'Marker', 'probe2transcript') NOT NULL,
  analysis_id           SMALLINT UNSIGNED NOT NULL,
  external_db_id        SMALLINT UNSIGNED,
  identifier            VARCHAR(255) NOT NULL,
  unmapped_reason_id    SMALLINT(5) UNSIGNED NOT NULL,
  query_score           DOUBLE,
  target_score          DOUBLE,
  ensembl_id            INT(10) UNSIGNED DEFAULT '0',
  ensembl_object_type   ENUM('RawContig','Transcript','Gene','Translation')
                        DEFAULT 'RawContig',
  parent                VARCHAR(255) DEFAULT NULL,
  PRIMARY KEY (unmapped_object_id),
  KEY id_idx (identifier),
  KEY anal_idx (analysis_id),
  KEY anal_exdb_idx (analysis_id, external_db_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'unmapped_reason'
#
# Describes the reason why a mapping failed.

CREATE TABLE unmapped_reason (

  unmapped_reason_id     SMALLINT(5) UNSIGNED NOT NULL AUTO_INCREMENT,
  summary_description    VARCHAR(255),
  full_description       VARCHAR(255),

  PRIMARY KEY (unmapped_reason_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

################################################################################
#
# Table structure for table 'ditag'
#
# Describes ditags

CREATE TABLE ditag (

       ditag_id          INT(10) UNSIGNED NOT NULL auto_increment,
       name              VARCHAR(30) NOT NULL,
       type              VARCHAR(30) NOT NULL,
       tag_count         smallint(6) UNSIGNED NOT NULL default 1,
       sequence          TINYTEXT NOT NULL,

       PRIMARY KEY (ditag_id)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;

################################################################################
#
# Table structure for table 'ditag_feature'
#
# Describes where ditags hit on the genome

CREATE TABLE ditag_feature (

       ditag_feature_id   INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
       ditag_id           INT(10) UNSIGNED NOT NULL default '0',
       ditag_pair_id      INT(10) UNSIGNED NOT NULL default '0',
       seq_region_id      INT(10) UNSIGNED NOT NULL default '0',
       seq_region_start   INT(10) UNSIGNED NOT NULL default '0',
       seq_region_end     INT(10) UNSIGNED NOT NULL default '0',
       seq_region_strand  TINYINT(1) NOT NULL default '0',
       analysis_id        SMALLINT UNSIGNED NOT NULL default '0',
       hit_start          INT(10) UNSIGNED NOT NULL default '0',
       hit_end            INT(10) UNSIGNED NOT NULL default '0',
       hit_strand         TINYINT(1) NOT NULL default '0',
       cigar_line         TINYTEXT NOT NULL,
       ditag_side         ENUM('F', 'L', 'R') NOT NULL,

       PRIMARY KEY  (ditag_feature_id),
       KEY (ditag_id),
       KEY (ditag_pair_id),
       KEY seq_region_idx (seq_region_id, seq_region_start, seq_region_end)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;

################################################################################
#
# Table structure for table 'unconventional_transcript_association'
#
# Describes transcripts that do not link to a single gene in the normal way.

CREATE TABLE unconventional_transcript_association (

       transcript_id    INT(10) UNSIGNED NOT NULL,
       gene_id          INT(10) UNSIGNED NOT NULL,
       interaction_type ENUM("antisense","sense_intronic","sense_overlaping_exonic","chimeric_sense_exonic"),

       KEY (transcript_id),
       KEY (gene_id)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;

################################################################################
#
# Table structure for seq_region mapping between releases
#
# Stores how the core seq_region_id have changed from release to release

CREATE TABLE seq_region_mapping (

	external_seq_region_id	INT(10) UNSIGNED NOT NULL,
	internal_seq_region_id	INT(10) UNSIGNED NOT NULL,
	mapping_set_id		INT(10) UNSIGNED NOT NULL,

	KEY (mapping_set_id)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;

################################################################################
#
# Table structure for seq_region mapping between releases
#
# Stores how which mapping group the seq_region are for a particular schema

CREATE TABLE mapping_set (

	mapping_set_id	INT(10)	UNSIGNED NOT NULL,
	schema_build	VARCHAR(20) NOT NULL,

	PRIMARY KEY(schema_build)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;


CREATE TABLE splicing_event (

  splicing_event_id       INT(10)  UNSIGNED NOT NULL AUTO_INCREMENT,
  name                    VARCHAR(134),
  gene_id                 INT(10) UNSIGNED NOT NULL,
  seq_region_id           INT(10) UNSIGNED NOT NULL,
  seq_region_start        INT(10) UNSIGNED NOT NULL,
  seq_region_end          INT(10) UNSIGNED NOT NULL,
  seq_region_strand       TINYINT(2) NOT NULL,
  type	                  ENUM('CNE','CE','AFE','A5SS','A3SS','MXE','IR','II','EI', 'AT', 'ALE', 'AI'),
  PRIMARY KEY (splicing_event_id),
  KEY gene_idx (gene_id),
  KEY seq_region_idx (seq_region_id, seq_region_start)

)  COLLATE=latin1_swedish_ci TYPE=MyISAM;


CREATE TABLE splicing_event_feature (

  splicing_event_feature_id INT(10)  UNSIGNED NOT NULL,
  splicing_event_id         INT(10)  UNSIGNED NOT NULL,
  exon_id                   INT(10)  UNSIGNED NOT NULL,
  transcript_id             INT(10)  UNSIGNED NOT NULL,
  feature_order             INT(10)  UNSIGNED NOT NULL,
  transcript_association    INT(10)  UNSIGNED NOT NULL,
  type                      ENUM('constitutive_exon','exon','flanking_exon'),
  start                     INT(10)  UNSIGNED NOT NULL,
  end                       INT(10)  UNSIGNED NOT NULL,

  PRIMARY KEY (splicing_event_feature_id,exon_id,transcript_id),
  KEY se_idx (splicing_event_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;



CREATE TABLE splicing_transcript_pair (


  splicing_transcript_pair_id INT(10)  UNSIGNED NOT NULL,
  splicing_event_id           INT(10)  UNSIGNED NOT NULL, 
  transcript_id_1             INT(10)  UNSIGNED NOT NULL,
  transcript_id_2             INT(10)  UNSIGNED NOT NULL,

  PRIMARY KEY (splicing_transcript_pair_id),
  KEY se_idx (splicing_event_id)
  

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

