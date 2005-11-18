# revisited schema naming issues
# Author: Arne Stabenau
# Date: 12.11.2001
#
# Glenn Proctor July 2003 - adapted for new schema structure
# 
# Note that more information about each table can be found in
# ensembl/docs/schema_description/

# Conventions:
#  - use lower case and underscores
#  - internal ids are integers named tablename_id
#  - same name is given in foreign key relations

CREATE TABLE affy_feature (
       affy_feature_id INT NOT NULL auto_increment,
       seq_region_id INT UNSIGNED NOT NULL,
       seq_region_start INT NOT NULL,
       seq_region_end INT NOT NULL,
       seq_region_strand TINYINT NOT NULL,
       
       mismatches TINYINT,
       affy_probe_id INT NOT NULL,
       analysis_id INT NOT NULL,

       PRIMARY KEY (affy_feature_id),
       KEY seq_region_idx( seq_region_id, seq_region_start ),
       KEY probe_idx( affy_probe_id )
) COLLATE=latin1_swedish_ci;		

CREATE TABLE affy_probe (
       affy_probe_id INT NOT NULL auto_increment,
       affy_array_id INT NOT NULL,
       probeset VARCHAR(40),
       name VARCHAR(20),

       PRIMARY KEY ( affy_probe_id, affy_array_id ),
       KEY probeset_idx( probeset ),
       KEY array_idx( affy_array_id )
) COLLATE=latin1_swedish_ci;

CREATE TABLE affy_array (
       affy_array_id INT NOT NULL auto_increment,
       parent_array_id INT,
       probe_setsize TINYINT NOT NULL,
       name VARCHAR(40) NOT NULL,

       PRIMARY KEY( affy_array_id )
) COLLATE=latin1_swedish_ci;


CREATE TABLE alt_allele (
  alt_allele_id INT NOT NULL auto_increment,
  gene_id INT NOT NULL,

  UNIQUE gene_idx( gene_id ),
  UNIQUE allele_idx( alt_allele_id, gene_id )
) COLLATE=latin1_swedish_ci;
  


################################################################################
#
# Table structure for table 'analysis'
# 
# semantics:
# analysis_id - internal id
# created   - date to distinguish newer and older versions off the 
#             same analysis. Not well maintained so far.
# logic_name  string to identify the analysis. Used mainly inside pipeline.
# db, db_version, db_file
#  - db should be a database name, db version the version of that db
#    db_file the file system location of that database, 
#    probably wiser to generate from just db and configurations
# program, program_version,program_file
#  - The binary used to create a feature. Similar semantic to above
# module, module_version
#  - Perl module names (RunnableDBS usually) executing this analysis.
# parameters a paramter string which is processed by the perl module
# gff_source, gff_feature 
#  - how to make a gff dump from features with this analysis


CREATE TABLE analysis (

  analysis_id                 int(10) unsigned NOT NULL auto_increment,
  created                     datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  logic_name                  varchar(40) not null,
  db                          varchar(120),
  db_version                  varchar(40),
  db_file                     varchar(120),
  program                     varchar(80),
  program_version             varchar(40),
  program_file                varchar(80),
  parameters                  varchar(255),
  module                      varchar(80),
  module_version              varchar(40),
  gff_source                  varchar(40),
  gff_feature                 varchar(40),

  PRIMARY KEY (analysis_id),
  KEY logic_name_idx( logic_name ),
  UNIQUE(logic_name)

) COLLATE=latin1_swedish_ci;


CREATE TABLE analysis_description (
  analysis_id	               int(10) unsigned NOT NULL,
  description                  text,
  display_label                varchar(255),

  KEY analysis_idx( analysis_id )
) COLLATE=latin1_swedish_ci;

################################################################################
#
# Table structure for table 'dna'
#
# This table stores DNA sequence.

CREATE TABLE dna (
  seq_region_id       int unsigned NOT NULL,
  sequence            mediumtext NOT NULL,

  PRIMARY KEY (seq_region_id)

) MAX_ROWS = 750000 AVG_ROW_LENGTH = 19000;

################################################################################
#
# Table structure for table 'dnac'
#
# Contains equivalent data to dna table, but 4 letters of DNA code are represented
# by a single binary character, based on 2 bit encoding
# do not need to worry about ambiguity of length, since this is stored in contig.length
# n_line column contains start-end pairs of coordinates in the string that are really Ns

CREATE TABLE dnac (
  seq_region_id  int unsigned NOT NULL,
  sequence  mediumblob NOT NULL,
  n_line    text,  

  PRIMARY KEY (seq_region_id)
) MAX_ROWS = 750000 AVG_ROW_LENGTH = 19000;

################################################################################
#
# Table structure for table 'exon'
#
# Note seq_region_start always less that seq_region_end, i.e.
# when the exon is on the other strand the seq_region_start
# is specifying the 3prime end of the exon.

CREATE TABLE exon (
 
  exon_id       	      int unsigned NOT NULL auto_increment,
  seq_region_id     	      int(10) unsigned NOT NULL,            # foreign key, seq_region:seq_region_id
  seq_region_start  	      int(10) unsigned NOT NULL,            # start of exon within seq_region
  seq_region_end              int(10) unsigned NOT NULL,            # end of exon within specified seq_region
  seq_region_strand           tinyint(2) NOT NULL,                  # 1 or -1 depending on the strand of the exon

  phase                       tinyint(2) NOT NULL,
  end_phase                   tinyint(2) NOT NULL,
  
  PRIMARY KEY (exon_id),
  KEY seq_region_idx (seq_region_id, seq_region_start )

) COLLATE=latin1_swedish_ci;

################################################################################
#
# Table structure for table 'exon_stable_id'
#

CREATE TABLE exon_stable_id (
   
  exon_id   		      int unsigned not null,       # foreign key exon:exon_id
  stable_id                   VARCHAR(128) not null,
  version                     int(10),
  created_date                DATETIME NOT NULL,
  modified_date               DATETIME NOT NULL,

  PRIMARY KEY( exon_id ),
  UNIQUE( stable_id, version )

) COLLATE=latin1_swedish_ci;

################################################################################
#
# Table structure for table 'exon_transcript'
#
# Note that the rank column indicates the 5prime to 3prime position of the 
# exon within the transcript ie rank of 1 means the exon is the 5prime most 
# within this transcript

CREATE TABLE exon_transcript (

  exon_id                     INT unsigned NOT NULL, # foreign key exon:exon_id
  transcript_id               INT unsigned NOT NULL, # foregin key transcript:transcript_id
  rank                        int(10) NOT NULL,        

  PRIMARY KEY (exon_id,transcript_id,rank),
  KEY transcript (transcript_id),
  KEY exon ( exon_id )

) COLLATE=latin1_swedish_ci;

################################################################################
#
# Table structure for table 'simple_feature'
#

CREATE TABLE simple_feature (

  simple_feature_id 	      int unsigned not null auto_increment,
  seq_region_id               int(10) unsigned NOT NULL,
  seq_region_start            int(10) unsigned NOT NULL,
  seq_region_end              int(10) unsigned NOT NULL,
  seq_region_strand           tinyint(1) NOT NULL,
  display_label               varchar(40) NOT NULL,
  analysis_id                 int(10) unsigned NOT NULL,
  score                       double,

  PRIMARY KEY ( simple_feature_id ),
  KEY seq_region_idx (seq_region_id, seq_region_start ),
  KEY analysis_idx( analysis_id ),
  KEY hit_idx( display_label )

) MAX_ROWS=100000000 AVG_ROW_LENGTH=80;

################################################################################
#
# Table structure for table 'protein_align_feature'
#

CREATE TABLE protein_align_feature (
 
  protein_align_feature_id    int unsigned not null auto_increment,
  seq_region_id               int(10) unsigned NOT NULL,
  seq_region_start            int(10) unsigned NOT NULL,
  seq_region_end              int(10) unsigned NOT NULL,
  seq_region_strand           tinyint(1) DEFAULT '1' NOT NULL,
  hit_start                   int(10) NOT NULL,
  hit_end                     int(10) NOT NULL,
  hit_name                    varchar(40) NOT NULL,
  analysis_id                 int(10) unsigned NOT NULL,
  score                       double,
  evalue                      double,
  perc_ident                  float,
  cigar_line                  text,

  PRIMARY KEY ( protein_align_feature_id ),
  KEY seq_region_idx( seq_region_id, analysis_id, seq_region_start, score ),
  KEY seq_region_idx_2( seq_region_id, seq_region_start),
  KEY hit_idx( hit_name ),
  KEY analysis_idx( analysis_id )

) MAX_ROWS=100000000 AVG_ROW_LENGTH=80;

################################################################################
#
# Table structure for table 'dna_align_feature'
#

CREATE TABLE dna_align_feature (

  dna_align_feature_id 	      int unsigned not null auto_increment,
  seq_region_id               int(10) unsigned NOT NULL,
  seq_region_start            int(10) unsigned NOT NULL,
  seq_region_end              int(10) unsigned NOT NULL,
  seq_region_strand           tinyint(1) NOT NULL,
  hit_start                   int NOT NULL,
  hit_end                     int NOT NULL,
  hit_strand                  tinyint(1) NOT NULL,
  hit_name                    varchar(40) NOT NULL,
  analysis_id                 int(10) unsigned NOT NULL,
  score                       double,
  evalue                      double,
  perc_ident                  float,
  cigar_line                  text,

  PRIMARY KEY ( dna_align_feature_id ),
  KEY seq_region_idx( seq_region_id, analysis_id, seq_region_start, score ),
  KEY seq_region_idx_2( seq_region_id, seq_region_start),
  KEY hit_idx( hit_name ),
  KEY analysis_idx( analysis_id )

) MAX_ROWS=100000000 AVG_ROW_LENGTH=80;

################################################################################
#
# Table structure for table 'repeat_consensus'
#

CREATE TABLE repeat_consensus (

  repeat_consensus_id  	      int unsigned NOT NULL auto_increment,
  repeat_name                 varchar(255) NOT NULL,
  repeat_class                varchar(100) NOT NULL,   # eg:  SINE, LINE, DNA Transposon,
                                                      # Retroviral LTR, Satellite,Tandem
  repeat_type                 varchar(40) NOT NULL,
  repeat_consensus            text,
  
  PRIMARY KEY( repeat_consensus_id ),
  KEY name (repeat_name),
  KEY class (repeat_class),
  KEY consensus(repeat_consensus(10)),
  KEY type( repeat_type )

) COLLATE=latin1_swedish_ci;

################################################################################
#
# Table structure for table 'repeat_feature'
#

CREATE TABLE repeat_feature (

  repeat_feature_id 	      int unsigned NOT NULL auto_increment,
  seq_region_id               int(10) unsigned NOT NULL,
  seq_region_start            int(10) unsigned NOT NULL,
  seq_region_end              int(10) unsigned NOT NULL,
  seq_region_strand           tinyint(1) DEFAULT '1' NOT NULL,
  repeat_start                int(10) NOT NULL,
  repeat_end                  int(10) NOT NULL,
  repeat_consensus_id         int(10) unsigned NOT NULL,
  analysis_id                 int(10) unsigned NOT NULL,
  score                       double,
  
  PRIMARY KEY (	repeat_feature_id ),
  KEY seq_region_idx( seq_region_id, seq_region_start ),
  KEY repeat_idx( repeat_consensus_id ),
  KEY analysis_idx( analysis_id )

) MAX_ROWS=100000000 AVG_ROW_LENGTH=80;

################################################################################
#
# Table structure for table 'gene'
#

CREATE TABLE gene (

  gene_id                     int unsigned NOT NULL auto_increment,
  biotype                     VARCHAR(40) NOT NULL,
  analysis_id                 int,
  seq_region_id               int(10) unsigned NOT NULL, 
  seq_region_start            int(10) unsigned NOT NULL, 
  seq_region_end              int(10) unsigned NOT NULL, 
  seq_region_strand           tinyint(2) NOT NULL,       
  display_xref_id             int unsigned,
  source                      VARCHAR(20) NOT NULL,
  status                      enum( 'KNOWN', 'NOVEL', 'PUTATIVE', 'PREDICTED' ),
  description                 text,

  PRIMARY KEY (gene_id),
  KEY seq_region_idx( seq_region_id, seq_region_start ),
  KEY xref_id_index ( display_xref_id ),
  KEY analysis_idx( analysis_id )

) COLLATE=latin1_swedish_ci;

################################################################################
#
# Table structure for table 'gene_stable_id'
#
CREATE TABLE gene_stable_id (

  gene_id 		      int unsigned not null,  # foreign key gene:gene_id
  stable_id                   VARCHAR(128) not null,
  version                     int(10),
  created_date                DATETIME NOT NULL,
  modified_date               DATETIME NOT NULL,

  PRIMARY KEY( gene_id ),
  UNIQUE( stable_id, version )

) COLLATE=latin1_swedish_ci;

################################################################################
#
# Table structure for table 'supporting_feature'
#

CREATE TABLE supporting_feature (

  exon_id 		      int(11) DEFAULT '0' NOT NULL,
  feature_type                enum('dna_align_feature','protein_align_feature'),
  feature_id                  int(11) DEFAULT '0' NOT NULL,

  UNIQUE all_idx (exon_id,feature_type,feature_id),
  KEY feature_idx (feature_type,feature_id)

) MAX_ROWS=100000000 AVG_ROW_LENGTH=80;

################################################################################
#
# Table structure for table 'transcript_supporting_feature'
#

CREATE TABLE transcript_supporting_feature (

  transcript_id 	      int(11) DEFAULT '0' NOT NULL,
  feature_type                enum('dna_align_feature','protein_align_feature'),
  feature_id                  int(11) DEFAULT '0' NOT NULL,

  UNIQUE all_idx (transcript_id,feature_type,feature_id),
  KEY feature_idx (feature_type,feature_id)

) MAX_ROWS=100000000 AVG_ROW_LENGTH=80;

################################################################################
#
# Table structure for table 'transcript'
#

CREATE TABLE transcript (

  transcript_id               INT UNSIGNED NOT NULL auto_increment,  
  gene_id                     INT UNSIGNED NOT NULL,  # foreign key gene:gene_id
  seq_region_id               int(10) unsigned NOT NULL, 
  seq_region_start            int(10) unsigned NOT NULL, 
  seq_region_end              int(10) unsigned NOT NULL, 
  seq_region_strand           tinyint(2) NOT NULL, 
  display_xref_id             int unsigned,
  biotype                     VARCHAR(40) NOT NULL,
  status                      enum( 'KNOWN', 'NOVEL', 'PUTATIVE', 'PREDICTED' ),
  description                 text,

  PRIMARY KEY (transcript_id),
  KEY seq_region_idx( seq_region_id, seq_region_start ),
  KEY gene_index (gene_id),
  KEY xref_id_index ( display_xref_id )

) COLLATE=latin1_swedish_ci;

################################################################################
#
# Table structure for table 'transcript_stable_id'
#

CREATE TABLE transcript_stable_id (

  transcript_id               int unsigned not null,  # foreign key transcript:transcript_id
  stable_id                   VARCHAR(128) not null,
  version                     int(10),
  created_date                DATETIME NOT NULL,
  modified_date               DATETIME NOT NULL,
  
  PRIMARY KEY( transcript_id ),
  UNIQUE( stable_id, version )

) COLLATE=latin1_swedish_ci;

################################################################################
#
# Table structure for table 'translation'
#
# The seq_start and seq_end are 1-based offsets into the *relative* coordinate 
# system of start_exon_id and end_exon_id. i.e. if the translation starts at the
# first base of the exon, seq_start would be 1

CREATE TABLE translation (

  translation_id              INT UNSIGNED NOT NULL auto_increment, 
  transcript_id               INT UNSIGNED NOT NULL, 
  seq_start                   INT(10) NOT NULL,       # relative to exon start
  start_exon_id               INT UNSIGNED NOT NULL,  # foreign key exon:exon_id
  seq_end                     INT(10) NOT NULL,       # relative to exon start
  end_exon_id                 INT UNSIGNED NOT NULL,  # foreign key exon:exon_id
  
  PRIMARY KEY (translation_id),
  KEY (transcript_id)
) COLLATE=latin1_swedish_ci;


################################################################################
#
# Table structure for table 'translation_stable_id'
#
CREATE TABLE translation_stable_id (

  translation_id 	      INT unsigned NOT NULL, # foreign key translation:translation_id
  stable_id                   VARCHAR(128) NOT NULL,
  version                     INT(10),
  created_date                DATETIME NOT NULL,
  modified_date               DATETIME NOT NULL,

  PRIMARY KEY( translation_id ),
  UNIQUE( stable_id, version )

) COLLATE=latin1_swedish_ci;


################################################################################
#
# Table structure for table 'assembly'
#
# This is a denormalised golden path. 
# The data in this table defines the "static golden path", i.e. the
# best effort draft full genome sequence as determined by the UCSC or NCBI
# (depending which assembly you are using)
#
# Each row represents a component, e.g. a contig,  (comp_seq_region_id, 
# FK from seq_region table) at least part of which is present in the golden path. 
# The part of the component that is in the path is delimited by fields cmp_start
# and cmp_end (start < end), and the absolute position within the golden path 
# chromosome (or other appropriate assembled structure) (asm_seq_region_id) is 
# given by asm_start and asm_end. 
# 

CREATE TABLE assembly (

  asm_seq_region_id           int unsigned NOT NULL,
  cmp_seq_region_id           int(10) unsigned NOT NULL, 
  asm_start                   int(10) NOT NULL,
  asm_end                     int(10) NOT NULL,
  cmp_start                   int(10) NOT NULL,
  cmp_end                     int(10) NOT NULL,
  ori                         tinyint  NOT NULL, 
  
  KEY(cmp_seq_region_id),
  KEY(asm_seq_region_id, asm_start)

) COLLATE=latin1_swedish_ci;

################################################################################
#
# Table structure for table 'protein_feature'
#

CREATE TABLE protein_feature (

  protein_feature_id          int(10) unsigned NOT NULL auto_increment,
  translation_id              int NOT NULL,	
  seq_start                   int(10) NOT NULL,
  seq_end                     int(10) NOT NULL,
  hit_start                   int(10) NOT NULL,
  hit_end                     int(10) NOT NULL,
  hit_id                      varchar(40) NOT NULL,
  analysis_id                 int(10) unsigned NOT NULL,
  score                       double NOT NULL,
  evalue                      double,
  perc_ident                  float,

  PRIMARY KEY   (protein_feature_id),
  KEY (translation_id),
  KEY hid_index ( hit_id ),
  KEY analysis_idx( analysis_id )

) COLLATE=latin1_swedish_ci;

################################################################################
#
# Table structure for table 'interpro'
#

CREATE TABLE interpro (

  interpro_ac	              varchar(40) NOT NULL,
  id		              varchar(40) NOT NULL,

  UNIQUE (interpro_ac, id),
  KEY (id)

) COLLATE=latin1_swedish_ci;


################################################################################
#
# Table structure for table 'karyotype'
#

CREATE TABLE karyotype (
  karyotype_id                int unsigned NOT NULL auto_increment,
  seq_region_id               int unsigned NOT NULL,
  seq_region_start            int(10)     NOT NULL,
  seq_region_end              int(10)     NOT NULL,
  band                        varchar(40) NOT NULL,
  stain                       varchar(40) NOT NULL,

  PRIMARY KEY (karyotype_id),
  KEY region_band_idx (seq_region_id,band)

) COLLATE=latin1_swedish_ci;

################################################################################
#
# Table structure for table 'object_xref'
#

CREATE TABLE object_xref (

  object_xref_id              INT not null auto_increment,
  ensembl_id                  int unsigned not null, 
  ensembl_object_type         ENUM( 'RawContig', 'Transcript', 'Gene', 'Translation', 'regulatory_factor', 'regulatory_feature' ) not null,
  xref_id                     INT unsigned not null,

  UNIQUE ( ensembl_object_type, ensembl_id, xref_id ),
  KEY oxref_idx( object_xref_id, xref_id, ensembl_object_type, ensembl_id ),
  KEY xref_idx(xref_id, ensembl_object_type)

) COLLATE=latin1_swedish_ci;

################################################################################
#
# Table structure for table 'identity_xref'
#

CREATE TABLE identity_xref(
        object_xref_id INT unsigned not null ,
	query_identity 	int(5),
        target_identity int(5),

	hit_start int,
	hit_end int,
	translation_start int,
	translation_end int,
	cigar_line text,
	
	score double,
	evalue double,
	analysis_id int,

  PRIMARY KEY (object_xref_id),
  KEY analysis_idx( analysis_id )
) COLLATE=latin1_swedish_ci;

################################################################################
#
# Table structure for table 'go_xref'
#

CREATE TABLE go_xref (

  object_xref_id int(10) unsigned DEFAULT '0' NOT NULL,
  linkage_type enum('IC', 'IDA', 'IEA', 'IEP', 'IGI', 'IMP', 
		    'IPI', 'ISS', 'NAS', 'ND', 'TAS', 'NR') NOT NULL,
  KEY (object_xref_id),
  UNIQUE(object_xref_id, linkage_type)

) COLLATE=latin1_swedish_ci;

################################################################################
#
# Table structure for table 'xref'
#

CREATE TABLE xref (

   xref_id 		      INT unsigned not null auto_increment,
   external_db_id             int not null,
   dbprimary_acc              VARCHAR(40) not null,
   display_label              VARCHAR(40) not null,
   version                    VARCHAR(10) DEFAULT '' NOT NULL,
   description                VARCHAR(255),

   PRIMARY KEY( xref_id ),
   UNIQUE KEY id_index( dbprimary_acc, external_db_id ),
   KEY display_index ( display_label )

) COLLATE=latin1_swedish_ci;

################################################################################
#
# Table structure for table 'external_synonym'
#

CREATE TABLE external_synonym (

  xref_id                     INT unsigned not null,
  synonym                     VARCHAR(40) not null,
  PRIMARY KEY( xref_id, synonym ),
  KEY name_index( synonym )

) COLLATE=latin1_swedish_ci;

################################################################################
#
# Table structure for table 'external_db' 
#

CREATE TABLE external_db (

  external_db_id 	      INT not null,
  db_name                     VARCHAR(27) NOT NULL,
  release                     VARCHAR(40) NOT NULL,
  status                      ENUM ('KNOWNXREF','KNOWN','XREF','PRED','ORTH', 'PSEUDO') not null,

  dbprimary_acc_linkable      BOOLEAN DEFAULT 1 NOT NULL,
  display_label_linkable      BOOLEAN DEFAULT 0 NOT NULL,

  priority                    INT NOT NULL,

  db_display_name             VARCHAR(255),

  PRIMARY KEY( external_db_id ) 

) COLLATE=latin1_swedish_ci;



CREATE TABLE prediction_exon (
    prediction_exon_id int unsigned not null auto_increment,
    prediction_transcript_id int unsigned not null,
    exon_rank smallint unsigned not null,
    seq_region_id int unsigned not null,
    seq_region_start int unsigned not null,
    seq_region_end int unsigned not null,
    seq_region_strand tinyint not null,
    start_phase tinyint not null,
    score double,
    p_value double,

    PRIMARY KEY( prediction_exon_id ),
    KEY (prediction_transcript_id),
    KEY ( seq_region_id, seq_region_start )
) COLLATE=latin1_swedish_ci;


CREATE TABLE prediction_transcript (
    prediction_transcript_id int unsigned not null auto_increment,
    seq_region_id int unsigned not null,
    seq_region_start int unsigned not null,
    seq_region_end int unsigned not null,
    seq_region_strand tinyint not null,
    analysis_id int,
    display_label varchar(255),
    
    PRIMARY KEY( prediction_transcript_id ),
    KEY ( seq_region_id, seq_region_start ),
    KEY analysis_idx( analysis_id )
) COLLATE=latin1_swedish_ci;



################################################################################
#
# Table structure for table 'meta' 
#

CREATE TABLE meta (

  meta_id 		      INT not null auto_increment,
  meta_key                    varchar( 40 ) not null,
  meta_value                  varchar( 255 ) not null,

  PRIMARY KEY( meta_id ),
  KEY meta_key_index ( meta_key ),
  KEY meta_value_index ( meta_value )

) COLLATE=latin1_swedish_ci;

# Auto add schema version to database

INSERT INTO meta (meta_key, meta_value) VALUES ("schema_version", "36");

################################################################################
#
# Table structure for table 'marker_synonym'

CREATE TABLE marker_synonym (

  marker_synonym_id           int unsigned not null auto_increment,
  marker_id                   int unsigned not null,  # foreign key marker:marker_id
  source                      varchar(20),
  name                        varchar(30),    

  PRIMARY KEY (marker_synonym_id),
  KEY marker_synonym_idx (marker_synonym_id, name),
  KEY marker_idx (marker_id)

) COLLATE=latin1_swedish_ci;


################################################################################
#
# Table structure for table 'marker'

CREATE TABLE marker (

  marker_id                   int unsigned not null auto_increment,
  display_marker_synonym_id   int unsigned, #foreign key marker_synonym:marker_synonym_id
  left_primer                 varchar(100) not null,
  right_primer                varchar(100) not null,
  min_primer_dist             int(10) unsigned not null,
  max_primer_dist             int(10) unsigned not null,
  priority                    int,
  type                        enum('est', 'microsatellite'),
  
  PRIMARY KEY (marker_id),
  KEY marker_idx (marker_id, priority)

) COLLATE=latin1_swedish_ci;

################################################################################
#
# Table structure for table 'marker_feature'

CREATE TABLE marker_feature (

  marker_feature_id           int unsigned not null auto_increment,
  marker_id                   int unsigned not null,     #foreign key marker:marker_id
  seq_region_id               int(10) unsigned NOT NULL, #foreign key contig:seq_region_id
  seq_region_start            int(10) unsigned NOT NULL,
  seq_region_end              int(10) unsigned NOT NULL,
  analysis_id                 int(10) unsigned NOT NULL, #foreign key analysis:analysis_id
  map_weight                  int(10) unsigned,

  PRIMARY KEY (marker_feature_id),
  KEY seq_region_idx (seq_region_id, seq_region_start ),
  KEY analysis_idx( analysis_id )

) COLLATE=latin1_swedish_ci;
   
################################################################################
#
# Table structure for table 'marker_map_location'
 
CREATE TABLE marker_map_location (

  marker_id                   int unsigned not null, #foreign key marker:marker_id
  map_id                      int unsigned not null, #foreign key map:map_id
  chromosome_name             varchar(15)  not null, 
  marker_synonym_id           int unsigned not null, #foreign key marker_synonym:marker_synonym_id
  position                    varchar(15) not null,
  lod_score                   double,
  
  PRIMARY KEY (marker_id, map_id),
  KEY map_idx( map_id, chromosome_name, position) 

) COLLATE=latin1_swedish_ci;

################################################################################
#
# Table structure for table 'map'

CREATE TABLE map (

  map_id                      int unsigned not null auto_increment, 
  map_name                    varchar(30) not null,

  PRIMARY KEY (map_id)
) COLLATE=latin1_swedish_ci;
    
################################################################################
#
# Table structure for table 'misc_feature'
#

CREATE TABLE misc_feature (

  misc_feature_id 	      int(10) unsigned NOT NULL auto_increment,
  seq_region_id               int(10) unsigned NOT NULL default '0',
  seq_region_start            int(10) unsigned NOT NULL default '0',
  seq_region_end              int(10) unsigned NOT NULL default '0',
  seq_region_strand           tinyint(4) NOT NULL default '0',

  PRIMARY KEY (misc_feature_id),
  KEY seq_region_idx( seq_region_id, seq_region_start )

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'misc_attrib'
#

CREATE TABLE misc_attrib (
  misc_feature_id             int(10) unsigned NOT NULL default '0',
  attrib_type_id              smallint(5) unsigned NOT NULL default '0',
  value                       varchar(255) NOT NULL default '',

  KEY type_val_idx( attrib_type_id, value ),
  KEY misc_feature_idx( misc_feature_id )
) COLLATE=latin1_swedish_ci TYPE=MyISAM;

################################################################################
#
# Table structure for table 'translation_attrib'
#

CREATE TABLE translation_attrib (
  translation_id              int(10) unsigned NOT NULL default '0',
  attrib_type_id              smallint(5) unsigned NOT NULL default '0',
  value                       varchar(255) NOT NULL default '',

  KEY type_val_idx( attrib_type_id, value ),
  KEY translation_idx( translation_id )
) COLLATE=latin1_swedish_ci TYPE=MyISAM;

################################################################################
#
# Table structure for table 'transcript_attrib'
#

CREATE TABLE transcript_attrib (
  transcript_id               int(10) unsigned NOT NULL default '0',
  attrib_type_id              smallint(5) unsigned NOT NULL default '0',
  value                       varchar(255) NOT NULL default '',

  KEY type_val_idx( attrib_type_id, value ),
  KEY transcript_idx( transcript_id )
) COLLATE=latin1_swedish_ci TYPE=MyISAM;




################################################################################
#
# Table structure for table 'seq_region_attrib'
#

CREATE TABLE seq_region_attrib (
  seq_region_id               int(10) unsigned NOT NULL default '0',
  attrib_type_id              smallint(5) unsigned NOT NULL default '0',
  value                       varchar(255) NOT NULL default '',

  KEY type_val_idx( attrib_type_id, value),
  KEY seq_region_idx (seq_region_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'attrib_type'
#

CREATE TABLE attrib_type (

  attrib_type_id              smallint(5) unsigned NOT NULL auto_increment,
  code                        varchar(15) NOT NULL default '',
  name                        varchar(255) NOT NULL default '',
  description                 text,

  PRIMARY KEY ( attrib_type_id),
  UNIQUE KEY c(code)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

################################################################################
#
# Table structure for table 'misc_set'
#

CREATE TABLE misc_set (

  misc_set_id                 smallint(5) unsigned NOT NULL auto_increment,
  code                        varchar(25) NOT NULL default '',
  name                        varchar(255) NOT NULL default '',
  description                 text NOT NULL,
  max_length                  int unsigned not null,

  PRIMARY KEY (misc_set_id),
  UNIQUE KEY c(code)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

################################################################################
#
# Table structure for table 'misc_feature_misc_set'
#

CREATE TABLE misc_feature_misc_set (

  misc_feature_id 	      int(10) unsigned NOT NULL default '0',
  misc_set_id 		      smallint(5) unsigned NOT NULL default '0',

  PRIMARY KEY ( misc_feature_id, misc_set_id ),
  KEY reverse_idx( misc_set_id, misc_feature_id )

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

################################################################################
#
# Tables for QTLs
#
################################################################################

################################################################################
#
# Table structure for table 'qtl'
#

CREATE TABLE qtl (

  qtl_id 		      int unsigned auto_increment not null,
  trait                       varchar(255) not null,
  lod_score                   float,
  flank_marker_id_1 	      int,
  flank_marker_id_2           int,
  peak_marker_id              int,

  PRIMARY KEY ( qtl_id ),
  KEY trait_idx( trait )

) COLLATE=latin1_swedish_ci;

################################################################################
#
# Table structure for table 'qtl_synonym'
#

CREATE TABLE qtl_synonym (

  qtl_synonym_id 	      int unsigned auto_increment not null,
  qtl_id                      int unsigned not null,
  source_database             enum("rat genome database", "ratmap") not null,
  source_primary_id           varchar(255) not null,

  PRIMARY KEY (qtl_synonym_id),
  KEY qtl_idx(qtl_id)

) COLLATE=latin1_swedish_ci;

################################################################################
#
# Table structure for table 'qtl_feature'
#

CREATE TABLE qtl_feature (

  seq_region_id 	      int not null,
  seq_region_start      int not null,
  seq_region_end        int not null,
  qtl_id                int not null,
  analysis_id           int not null,

  KEY( qtl_id ),
  KEY loc_idx( seq_region_id, seq_region_start ),
  KEY analysis_idx( analysis_id )
) COLLATE=latin1_swedish_ci;

################################################################################
#
# Tables for stable ID mapping tracking
#
################################################################################

################################################################################
#
# Table structure for table 'mapping_session'
#

CREATE TABLE mapping_session (

  mapping_session_id 	      int(11) NOT NULL auto_increment,
  old_db_name                 varchar(80) NOT NULL default '',
  new_db_name                 varchar(80) NOT NULL default '',
  created                     timestamp(14) NOT NULL,

  PRIMARY KEY  (mapping_session_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

################################################################################
#
# Table structure for table 'stable_id_event'
#

CREATE TABLE stable_id_event (

  old_stable_id 	      varchar(128),
  old_version                 smallint,
  new_stable_id               varchar(128),
  new_version                 smallint,
  mapping_session_id          int(11) NOT NULL default '0',
  type                        ENUM('gene', 'transcript', 'translation') NOT NULL,

   UNIQUE KEY uni_idx (mapping_session_id, old_stable_id, old_version, new_stable_id, new_version, type),

  KEY new_idx (new_stable_id),
  KEY old_idx (old_stable_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

################################################################################
#
# Table structure for table 'gene_archive'
#

CREATE TABLE gene_archive (

  gene_stable_id	      VARCHAR(128) NOT NULL,
  gene_version                smallint NOT NULL,
  transcript_stable_id        VARCHAR(128) NOT NULL,
  transcript_version          smallint NOT NULL,
  translation_stable_id       VARCHAR(128) NOT NULL,
  translation_version         smallint NOT NULL,
  mapping_session_id          int NOT NULL,

  KEY gene_idx( gene_stable_id, gene_version ),
  KEY transcript_idx( transcript_stable_id, transcript_version ),
  KEY translation_idx( translation_stable_id, translation_version )

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

################################################################################
#
# Table structure for table 'peptide_archive'
#

CREATE TABLE peptide_archive (

  translation_stable_id       VARCHAR(128) NOT NULL,
  translation_version         smallint NOT NULL,
  peptide_seq                 mediumtext NOT NULL,

  PRIMARY KEY( translation_stable_id, translation_version )

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

################################################################################
#
# Table structure for table 'seq_region'
#

CREATE TABLE seq_region (

  seq_region_id               INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  name                        VARCHAR(40) NOT NULL,
  coord_system_id             INT(10) NOT NULL,
  length                      INT(10) NOT NULL,

  UNIQUE(coord_system_id, name),
  PRIMARY KEY (seq_region_id),
  KEY name_idx(name)  

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'assembly_exception'
#

CREATE TABLE assembly_exception (

  assembly_exception_id       INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id               INT NOT NULL,
  seq_region_start            INT NOT NULL,
  seq_region_end              INT NOT NULL, 
  exc_type                    ENUM('HAP', 'PAR') NOT NULL,
  exc_seq_region_id           INT NOT NULL, 
  exc_seq_region_start        INT NOT NULL, 
  exc_seq_region_end          INT NOT NULL,
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

  coord_system_id             INT NOT NULL auto_increment,
  name                        VARCHAR(40) NOT NULL,
  version                     VARCHAR(40),
  rank                        INT NOT NULL,
  attrib                      SET ('default_version', 'sequence_level'),

  UNIQUE(name, version),
  UNIQUE(rank),
  PRIMARY KEY (coord_system_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

################################################################################
#
# Table structure for table 'meta_coord'
#

CREATE TABLE meta_coord (

  table_name                  VARCHAR(40) NOT NULL,
  coord_system_id             INT NOT NULL,
  max_length                  INT,

  UNIQUE(table_name, coord_system_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;



CREATE TABLE density_feature (
  density_feature_id    INT NOT NULL auto_increment,
  density_type_id       INT NOT NULL, #FK refs density_type
  seq_region_id         INT NOT NULL, #FK refs seq_region
  seq_region_start      INT NOT NULL,
  seq_region_end        INT NOT NULL,
  density_value         FLOAT NOT NULL,

  PRIMARY KEY(density_feature_id),
  KEY seq_region_idx (density_type_id, seq_region_id, seq_region_start),
  KEY seq_region_id_idx (seq_region_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;



CREATE TABLE density_type (
  density_type_id       INT NOT NULL auto_increment,
  analysis_id           INT NOT NULL, #FK refs analysis
  block_size            INT NOT NULL,
  region_features       INT NOT NULL,
  value_type            ENUM('sum','ratio') NOT NULL,
  PRIMARY KEY(density_type_id),
  UNIQUE(analysis_id, block_size, region_features)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

################################################################################
#
# Table structure for table 'regulatory_feature'
#
# Describes instances of regulatory_factor binding to the genome.

CREATE TABLE regulatory_feature (

  regulatory_feature_id INT NOT NULL auto_increment,
  name                  VARCHAR(255) NOT NULL,
  seq_region_id         INT NOT NULL,                  # FK refs seq_region
  seq_region_start      INT NOT NULL,
  seq_region_end        INT NOT NULL,
  seq_region_strand     TINYINT NOT NULL,
  analysis_id           INT NOT NULL,                  # FK refs analysis
  regulatory_factor_id  INT,                           # FK refs regulatory_factor


  PRIMARY KEY(regulatory_feature_id),
  KEY seq_region_idx(seq_region_id, analysis_id, seq_region_start),
  KEY seq_region_idx_2(seq_region_id, seq_region_start)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

################################################################################
#
# Table structure for table 'regulatory_factor'
#

CREATE TABLE regulatory_factor (

  regulatory_factor_id   INT NOT NULL auto_increment,
  name                   VARCHAR(255) NOT NULL,
  type                   ENUM('miRNA_target', 'transcription_factor', 'transcription_factor_complex'),

  PRIMARY KEY(regulatory_factor_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

################################################################################
#
# Table structure for table 'regulatory_feature_object'
#
# Relates regulatory regions to the Ensembl objects they influence. Many-many.

CREATE TABLE regulatory_feature_object (

  regulatory_feature_id INT NOT NULL,               # FK to regulatory_feature
  ensembl_object_type   ENUM( 'Transcript', 'Translation', 'Gene') NOT NULL,
  ensembl_object_id     INT NOT NULL,               # FK to transcript,gene etc
  influence             ENUM('positive', 'negative', 'mixed', 'unknown'),
  evidence              VARCHAR(255),

  KEY regulatory_feature_idx (regulatory_feature_id),
  KEY ensembl_object_idx (ensembl_object_type, ensembl_object_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;


################################################################################
#
# Table structure for table 'regulatory_factor_coding'
#
# Describes which genes/transcripts code for particular regulatory factors.

CREATE TABLE regulatory_factor_coding (

  regulatory_factor_id  INT NOT NULL,      # FK to regulatory_factor
  transcript_id         INT,               # FK to transcript
  gene_id         	INT,               # FK to gene

  KEY transcript_idx (transcript_id),
  KEY gene_idx (gene_id),
  KEY regulatory_factor_idx (regulatory_factor_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

################################################################################
#
# Table structure for table 'regulatory_search_region'
#
# Describes regions which were searched for regulatory features.

CREATE TABLE regulatory_search_region (

  regulatory_search_region_id  INT NOT NULL auto_increment,
  name                         VARCHAR(255) NOT NULL,
  seq_region_id                INT NOT NULL,                 # FK refs seq_region
  seq_region_start             INT NOT NULL,
  seq_region_end               INT NOT NULL,
  seq_region_strand            TINYINT NOT NULL,
  ensembl_object_type          ENUM( 'Transcript', 'Translation', 'Gene') NOT NULL,
  ensembl_object_id            INT,           # FK to gene/transcript/translation
  analysis_id                  INT NOT NULL,  # FK to analysis

  PRIMARY KEY (regulatory_search_region_id),
  KEY rsr_idx (regulatory_search_region_id),
  KEY ensembl_object_idx (ensembl_object_type, ensembl_object_id),
  KEY seq_region_idx(seq_region_id, seq_region_start),
  KEY seq_region_idx_2(seq_region_id, seq_region_start)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

