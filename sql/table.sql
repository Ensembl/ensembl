# revisited schema naming issues
# Author: Arne Stabenau
# Date: 12.11.2001

# conventions
# use lower case and underscores
# internal ids are integers named tablename_id
# same name is given in foreign key relations


#
# Table structure for table 'analysis'
#
# 

CREATE TABLE analysis (
  analysis_id int(10) unsigned DEFAULT '0' NOT NULL auto_increment,
  created datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  logic_name varchar(40) not null,
  db varchar(120),
  db_version varchar(40),
  db_file varchar(120),
  program varchar(80),
  program_version varchar(40),
  program_file varchar(80),
  parameters varchar(80),
  module varchar(80),
  module_version varchar(40),
  gff_source varchar(40),
  gff_feature varchar(40),

  PRIMARY KEY (analysis_id)
);

# semantics
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


#
# Table structure for table 'chromosome'
#
CREATE TABLE chromosome (
  chromosome_id     tinyint unsigned NOT NULL auto_increment,
  name              varchar(40) NOT NULL,
  id                int(11) NOT NULL,
  known_genes       int(11) NULL,
  unknown_genes     int(11) NULL,
  snps              int(11) NULL,
  length            int(11) NULL,
  
  PRIMARY KEY (chromosome_id)
);
# Statistics for the web site


#
# Table structure for table 'clone'
#
CREATE TABLE clone (
  clone_id      int(10) unsigned NOT NULL auto_increment,
  name          varchar(40) NOT NULL,
  embl_acc      varchar(40) NOT NULL,
  version       int(10) NOT NULL,
  embl_version  int(10) NOT NULL,
  htg_phase     int(10) DEFAULT '-1' NOT NULL,
  created       datetime NOT NULL,
  modified      datetime NOT NULL,
  
  PRIMARY KEY (clone_id),
  KEY embl (embl_acc,embl_version),
  KEY id   (name, version)
);

# semantics
# id - string we give to this clone in ensembl
#      should be same as embl_id unless clone is not in embl
# embl_id - hows the clone submitted to embl
# htg_phase - finished/unfinished: draft is 123, finished is 4


#
# Table structure for table 'map_density'
#
CREATE TABLE map_density (
   chromosome_id    tinyint unsigned NOT NULL,
   chr_start	    int(10) NOT NULL,
   chr_end	    int(10) NOT NULL,
   type		    varchar(20) NOT NULL,
   value	    int(10) NOT NULL,
    
   PRIMARY KEY(type,chromosome_id,chr_start) 
);

#
# Table structure for table 'contig'
#
CREATE TABLE contig (
  contig_id         int(10) unsigned NOT NULL auto_increment,
  name              varchar(40) NOT NULL,
  clone_id          int(10) NOT NULL,
  length            int(10) unsigned NOT NULL,   # foreign key clone:internal_id
  offset            int(10) unsigned,
  corder            int(10) unsigned,
  dna_id            int(10) NOT NULL,            # foreign key dna:id
  international_name  varchar(40),
  
  PRIMARY KEY (contig_id),
  UNIQUE name (name),
  KEY clone (clone_id),
  KEY dna (dna_id)
);



#
# Table structure for table 'dna'
#

# This table holds the sequence of the contigs from the contig table.
# The sequence is that of the contig, not that of the golden
# path. I.e. to construct the golden path from the dna entries,
# the sequence of contigs with an orientation of -1 must be
# reversed and bases complemented. The assembly
# table has the contig orientation (raw_ori).
# Note the length of the dna.sequence field is always equal
# to the appropriate length field in the contig table
# (probably a violation of some form of normal form since contig.length
# is an attibute of the dna.sequence field)

CREATE TABLE dna (
  dna_id    int(10) unsigned NOT NULL auto_increment,
  sequence  mediumtext NOT NULL,
  created   datetime NOT NULL,
  
  PRIMARY KEY (dna_id)
) MAX_ROWS = 750000 AVG_ROW_LENGTH = 13000;

#
# Table structure for table 'exon'
#

# Note seq_start always less that seq_end, i.e.
# when the exon is on the other strand the seq_start
# is specifying the 3' end of the exon.

# The Sticky Rank differentiates between fragments of
# the same exon. I.e for exons that
# span multiple contigs, all the fragments
# are in this table with the same id,
# but different sticky_rank values

CREATE TABLE exon (
  exon_id       int unsigned NOT NULL auto_increment,
  contig_id     int(10) unsigned NOT NULL,            # foreign key, contig:internal_id
  contig_start  int(10) unsigned NOT NULL,                     # start of exon within contig
  contig_end    int(10) unsigned NOT NULL,                     # end of exon within specified contig
  contig_strand tinyint(2) NOT NULL,                  # 1 or -1 depending on the strand of the exon

  phase         tinyint(2) NOT NULL,
  end_phase     tinyint(2) NOT NULL,
  sticky_rank   tinyint DEFAULT '1' NOT NULL,         # see note above
  
  PRIMARY KEY ( exon_id, sticky_rank),
  KEY contig_idx (contig_id, contig_start )
);

CREATE TABLE exon_stable_id (
    exon_id   int unsigned not null,       # foreign key exon:exon_id
    stable_id VARCHAR(40) not null,
    version   int(10) DEFAULT '1' NOT NULL,
    created   datetime NOT NULL,
    modified  datetime NOT NULL,
    
    PRIMARY KEY( exon_id ),
    UNIQUE( stable_id, version )
);



#
# Table structure for table 'exon_transcript'
#
CREATE TABLE exon_transcript (
  exon_id          INT unsigned NOT NULL, # foreign key exon:exon_id
  transcript_id    INT unsigned NOT NULL, # foregin key transcript:transcript_id
  rank          int(10) NOT NULL,         # Indicates the 5' to 3' position of the exon
                                          # within the transcript ie rank of 1 means
                                          # the exon is the 5' most within this transcript
  
  PRIMARY KEY (exon_id,transcript_id,rank),
  KEY transcript (transcript_id)
);



CREATE TABLE simple_feature (
  simple_feature_id int unsigned not null auto_increment,
  contig_id int(10) unsigned NOT NULL,
  contig_start int(10) unsigned NOT NULL,
  contig_end int(10) unsigned NOT NULL,
  contig_strand tinyint(1) NOT NULL,
  display_label varchar(40) NOT NULL, # what to show, may link to other things, depends on analysis
  analysis_id int(10) unsigned NOT NULL,

# What scoring do we need ?

  score double,

  PRIMARY KEY ( simple_feature_id ),
  KEY contig_idx( contig_id, analysis_id, contig_start ),
  KEY analysis_idx( analysis_id, contig_id ),
  KEY hit_idx( display_label )
) max_rows=300000000 avg_row_length=80;


CREATE TABLE protein_align_feature (
  protein_align_feature_id int unsigned not null auto_increment,
  contig_id int(10) unsigned NOT NULL,
  contig_start int(10) unsigned NOT NULL,
  contig_end int(10) unsigned NOT NULL,
  analysis_id int(10) unsigned NOT NULL,
  contig_strand tinyint(1) DEFAULT '1' NOT NULL,
  hit_start int(10) NOT NULL,
  hit_end int(10) NOT NULL,
  hit_name varchar(40) NOT NULL,
  cigar_line varchar(255),

  # cigar might be too small
  #  What scoring do we need ?

  evalue double,
  perc_ident float,
  score double,

  PRIMARY KEY (	protein_align_feature_id ),
  KEY hit_idx( hit_name ),
  KEY ctg_idx( contig_id, contig_start, analysis_id )
) max_rows=300000000 avg_row_length=80;


CREATE TABLE dna_align_feature (
  dna_align_feature_id int unsigned not null auto_increment,
  contig_id int(10) unsigned NOT NULL,
  analysis_id int(10) unsigned NOT NULL,

  contig_start int(10) unsigned NOT NULL,
  contig_end int(10) unsigned NOT NULL,
  contig_strand tinyint(1) NOT NULL,

  hit_start int NOT NULL,
  hit_end int NOT NULL,
  hit_name varchar(40) NOT NULL,
  hit_strand tinyint(1) NOT NULL,

  phase tinyint(1), # values 0,1,2. Only used if translations are involved
  cigar_line varchar(255),
  
#  What scoring do we need ?

  evalue double,
  perc_ident float,
  score double,

  PRIMARY KEY ( dna_align_feature_id ),
  KEY hit_idx( hit_name ),
  KEY ctg_idx( contig_id, analysis_id, contig_start )
) max_rows=300000000 avg_row_length=80;


CREATE TABLE repeat_consensus (
    repeat_id           int unsigned NOT NULL auto_increment,
    repeat_name         varchar(255) NOT NULL,
    repeat_class        varchar(40) NOT NULL,   # eg:  SINE, LINE, DNA Transposon,
                                                # Retroviral LTR, Satellite,Tandem
    repeat_consensus    text,   # Or dna_id with entry in DNA table?
    
    PRIMARY KEY( repeat_id )
);


CREATE TABLE repeat_feature (
  repeat_feature_id int unsigned NOT NULL auto_increment,
  contig_id int(10) unsigned NOT NULL,
  contig_start int(10) unsigned NOT NULL,
  contig_end int(10) unsigned NOT NULL,
  contig_strand tinyint(1) DEFAULT '1' NOT NULL,
  analysis_id int(10) unsigned NOT NULL,
  repeat_start int(10) NOT NULL,
  repeat_end int(10) NOT NULL,
  repeat_id int(10) unsigned NOT NULL,


#  What scoring do we need ?

  score double,
  
  PRIMARY KEY (	repeat_feature_id ),
  KEY contig_idx( contig_id, contig_start, analysis_id ),
  KEY repeat_idx( repeat_id, contig_id, contig_start )
) max_rows=300000000 avg_row_length=80;

#
# Table structure for table 'gene'
#
CREATE TABLE gene (
  gene_id   int unsigned not null auto_increment,
  type VARCHAR(40) not null,
  analysis_id int,
  transcript_count int not null,

  PRIMARY KEY (gene_id)
);


#
# Experimental table to decouple locations from located objects
#  Genes will have locations here. Others might follow if performance is ok.
#
# A sister fragment_location table might be constructed if this experiment is 
#  successful. This would enable better support for unassembled genome.
#
CREATE TABLE assembly_locations (
  type enum( 'gene' ) not null,
  type_id int unsigned not null,

  assembly_type varchar(20) not null,
  chromosome_id tinyint unsigned not null,
  chr_start int not null,
  chr_end int not null,
  chr_strand tinyint not null,

  KEY loc_idx( assembly_type, chromosome_id, chr_start, chr_end, chr_strand, type, type_id ),
  KEY obj_idx( type, type_id, assembly_type, chromosome_id, chr_start, chr_end, chr_strand )
);  


CREATE TABLE gene_stable_id (
    gene_id int unsigned not null,    # foreign key gene:gene_id
    stable_id VARCHAR(40) not null,
    version   int(10) DEFAULT '1' NOT NULL,
    created   datetime NOT NULL,
    modified  datetime NOT NULL,
    
    PRIMARY KEY( gene_id ),
    UNIQUE( stable_id, version )
);

# what are we going to do with this ??????

#
# Table structure for table 'supporting_feature'
#  try use the exon_feature table

CREATE TABLE supporting_feature (
  supporting_feature_id            int(10) unsigned NOT NULL auto_increment,
  exon_id          int NOT NULL,             # foreign key exon:exon_id
  contig_id     int(10) unsigned NOT NULL,
  contig_start  int(10) NOT NULL,
  contig_end    int(10) NOT NULL,
  score         int(10) NOT NULL,
  strand        int(1) DEFAULT '1' NOT NULL,
  analysis_id   int(10) unsigned NOT NULL,
  hit_start     int(11) NOT NULL,
  hit_end       int(11) NOT NULL,
  hit_id        varchar(40) NOT NULL,
  evalue        double,
  perc_ident    float,
  phase         tinyint(1),
  end_phase     tinyint(1),
  hit_strand    tinyint(1),
  
  PRIMARY KEY (supporting_feature_id),
  KEY exon( exon_id ),                   # most used index here!
  KEY analysis (contig_id, analysis_id),
  KEY hid (hit_id)
);


# Avoid feature duplication by just linking exons with their
#  belonging feature. Do we need to find exons by feature??

CREATE TABLE exon_feature (
  exon_id int unsigned not null,
  feature_id int unsigned not null,
  feature_type enum( "dna_align", "simple", "protein_align" ) not null,

  key exon_idx( exon_id, feature_id, feature_type )
) max_rows=300000000 avg_row_length=80;
  
 


#
# Table structure for table 'transcript'
#
CREATE TABLE transcript (
  transcript_id    INT UNSIGNED NOT NULL auto_increment,  
  gene_id          INT UNSIGNED NOT NULL,          # foreign key gene:gene_id
  translation_id   INT UNSIGNED NOT NULL,          # foreign key translation:translation_id
  exon_count int not null,

  PRIMARY KEY (transcript_id),
  KEY gene_index (gene_id),
  KEY translation_index ( translation_id )		
);

CREATE TABLE transcript_stable_id (
    transcript_id int unsigned not null,  # foreign key transcript:transcript_id
    stable_id     VARCHAR(40) not null,
    version       int(10) DEFAULT '1' NOT NULL,
    
    PRIMARY KEY( transcript_id ),
    UNIQUE( stable_id, version )
);


#
# Table structure for table 'translation'
#

# The seq_start and seq_end are 1-based offsets into the
# *relative* coordinate system of start_exon_id and end_exon_id.
# ie, if the translation starts at the first base of the exon, seq_start 
# would be 1

CREATE TABLE translation (
  translation_id  INT UNSIGNED NOT NULL auto_increment, 
  seq_start       INT(10) NOT NULL, # relative to exon start
  start_exon_id   INT UNSIGNED NOT NULL,  # foreign key exon:exon_id
  seq_end         INT(10) NOT NULL, # relative to exon start
  end_exon_id     INT UNSIGNED NOT NULL,  # foreign key exon:exon_id
  
  PRIMARY KEY (translation_id)
);

CREATE TABLE translation_stable_id (
    translation_id INT unsigned NOT NULL, # foreign key translation:translation_id
    stable_id VARCHAR(40) NOT NULL,
    version   INT(10) DEFAULT '1' NOT NULL,
    
    PRIMARY KEY( translation_id ),
    UNIQUE( stable_id, version )
);

# this is a denormalised golden path

#
# The data in this table defines the "static golden path", i.e. the
# best effort draft full genome sequence as determined by the UCSC or NCBI
# (depending which assembly you are using)
#
# Each row represents a contig (raw_id, FK from contig table) at least part of
# which is present in the golden path. The part of the contig that is
# in the path is delimited by fields raw_start and raw_end (start < end), and
# the absolute position within the golden path chromosome (chromosome_id) is given
# by chr_start and chr_end. Each contig is in some "supercontig" such as a
# "fingerprint clone contig" or NT contig and the super contig is identified
# by field superctg_name and the position of the specified bit of the contig
# within its super contig is given by fields superctg_start and superctg_end.
# With the data set at time of this writing, field type is always "NCBI_xx".
# 
# NB, chr_start <= chr_end, raw_start <= raw_end, and superctg_start <= superctg_end.
# 
 

CREATE TABLE assembly (
    superctg_name    varchar(20) NOT NULL,
    chromosome_id  tinyint unsigned  NOT NULL,
    contig_id      int(10) unsigned NOT NULL, # foreign key contig:internal_id
    chr_start      int(10) NOT NULL,
    chr_end        int(10) NOT NULL,
    superctg_start   int(10) NOT NULL,
    superctg_end     int(10) NOT NULL,
    superctg_ori     tinyint(2) NOT NULL,
    contig_start   int(10) NOT NULL,
    contig_end     int(10) NOT NULL,
    contig_ori     tinyint  NOT NULL, 
    type           varchar(20) NOT NULL,
    
    PRIMARY KEY(contig_id,type),
    KEY(superctg_name, superctg_start),
    KEY(chromosome_id,chr_start) 
);


#
# Table structure for table 'protein_feature'
#

CREATE TABLE protein_feature (
  protein_feature_id  int(10) unsigned NOT NULL auto_increment,
  translation_id int NOT NULL,	
  seq_start     int(10) NOT NULL,
  seq_end       int(10) NOT NULL,
  analysis_id      int(10) unsigned NOT NULL,
  hit_start        int(10) NOT NULL,
  hit_end          int(10) NOT NULL,
  hit_id           varchar(40) NOT NULL,
  score         double NOT NULL,
  evalue        double,
  perc_ident    float,

  PRIMARY KEY   (protein_feature_id),
  KEY (translation_id),
  KEY hid_index ( hit_id )
);

#
#Table structure for table 'interpro'
#

CREATE TABLE interpro (
  interpro_ac	varchar(40) NOT NULL,
  id		varchar(40) NOT NULL,
  KEY (interpro_ac),
  KEY (id)
);


#
#Table structure for table gene_description
#

CREATE TABLE gene_description (
  gene_id     int unsigned NOT NULL,
  description varchar(255),
  PRIMARY KEY (gene_id)
);

CREATE TABLE karyotype (
   chromosome_id  tinyint unsigned NOT NULL,
   chr_start      int(10)     NOT NULL,
   chr_end        int(10)     NOT NULL,
   band           varchar(40) NOT NULL,
   stain          varchar(40) NOT NULL,
   PRIMARY KEY (chromosome_id,band)
);


#
#Table structure for table objectXref
#

CREATE TABLE object_xref(
       object_xref_id INT not null auto_increment,
       ensembl_id int unsigned not null, 
       ensembl_object_type ENUM( 'RawContig', 'Transcript', 'Gene', 'Translation' ) not null,
       xref_id INT unsigned not null,

       UNIQUE ( ensembl_object_type, ensembl_id, xref_id ),
       KEY xref_index( object_xref_id, xref_id, ensembl_object_type, ensembl_id )
);

#
#Table structure for identity_xref
#

CREATE TABLE identity_xref(
        object_xref_id INT unsigned not null ,
	query_identity 	int(5),
        target_identity int(5),
        PRIMARY KEY (object_xref_id)
);


#
#Table structure for table xref
#

CREATE TABLE xref (
         xref_id INT unsigned not null auto_increment,
         external_db_id int not null,
         dbprimary_acc VARCHAR(40) not null,
	 display_label VARCHAR(40) not null,
         version VARCHAR(10) DEFAULT '' NOT NULL,
	 description VARCHAR(255),

         PRIMARY KEY( xref_id ),
         UNIQUE KEY id_index( dbprimary_acc, external_db_id ),
         KEY display_index ( display_label )
);


#
#Table structure for table externalSynonym
#

CREATE TABLE external_synonym(
         xref_id INT unsigned not null,
         synonym VARCHAR(40) not null,
         PRIMARY KEY( xref_id, synonym ),
	 KEY name_index( synonym )
   	);

#
#Table structure for table externalDB 
#

CREATE TABLE external_db(
         external_db_id INT not null auto_increment,
         db_name VARCHAR(40) not null,
	 release VARCHAR(40) DEFAULT '' NOT NULL,
	 status  ENUM ('KNOWN','XREF','PRED') not null,
         PRIMARY KEY( external_db_id ) 
);


#
# Table structure for table 'landmarkMarker'
#
CREATE TABLE landmark_marker (
  marker char(40) DEFAULT '' NOT NULL,
  name char(40) DEFAULT '' NOT NULL,
  chr_start bigint(17) DEFAULT '0' NOT NULL,
  chr_end bigint(17) DEFAULT '0' NOT NULL,
  chr_strand bigint(1) DEFAULT '0' NOT NULL,
  chromosome_id tinyint unsigned NOT NULL,
  KEY chromosome_id (chromosome_id,chr_start)
);

CREATE TABLE meta (
    meta_id INT not null auto_increment,
    meta_key varchar( 40 ) not null,
    meta_value varchar( 255 ) not null,

    PRIMARY KEY( meta_id ),
    KEY meta_key_index ( meta_key ),
    KEY meta_value_index ( meta_value )
);

CREATE TABLE prediction_transcript (
    prediction_transcript_id int unsigned not null auto_increment,
    exon_rank smallint unsigned not null,
    contig_id int unsigned not null,
    contig_start int unsigned not null,
    contig_end int unsigned not null,
    contig_strand tinyint not null,
    start_phase tinyint not null,
    score double,
    p_value double,
    analysis_id int,
    exon_count smallint,

    PRIMARY KEY( prediction_transcript_id, exon_rank ),
    KEY contig_idx( contig_id, contig_start )
);

    

# Auto add schema version to database
insert into meta (meta_key, meta_value) values ("schema_version", "$Revision$");


