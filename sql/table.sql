# MySQL dump 5.13
#
# Host: obi-wan    Database: ens500
#--------------------------------------------------------
# Server version	3.22.32

#
# Table structure for table 'analysisprocess'
#
CREATE TABLE analysisprocess (
  analysisId int(10) unsigned DEFAULT '0' NOT NULL auto_increment,
  created datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  logic_name varchar(40) not null,
  db varchar(120),
  db_version varchar(40),
  db_file varchar(120),
  program varchar(80),
  program_version varchar(40),
  program_file varchar(40),
  parameters varchar(80),
  module varchar(80),
  module_version varchar(40),
  gff_source varchar(40),
  gff_feature varchar(40),

  PRIMARY KEY (analysisId),
  index(gff_source),
  index(gff_feature),
  index(db)
);

#
# Table structure for table 'analysis'
#
CREATE TABLE analysis (
  id                int(10) unsigned NOT NULL auto_increment,
  db                varchar(40),
  db_version        varchar(40),
  program           varchar(40) NOT NULL,
  program_version   varchar(40),
  gff_source        varchar(40),
  gff_feature       varchar(40),
  
  PRIMARY KEY (id)
);

#
# Table structure for table 'chromosome'
#
CREATE TABLE chromosome (
  chromosome_id     int(10) unsigned NOT NULL auto_increment,
  name              varchar(40) NOT NULL,
  species_id        int(11) unsigned NOT NULL,
  id                int(11) unsigned NOT NULL,
  known_genes       int(11) unsigned NULL,
  unknown_genes     int(11) unsigned NULL,
  snps              int(11) unsigned NULL,
  length            int(11) unsigned NULL,
  
  PRIMARY KEY (chromosome_id)
);

#
# Table structure for table 'clone'
#
CREATE TABLE clone (
  internal_id   int(10) unsigned NOT NULL auto_increment,
  id            varchar(40) NOT NULL,
  embl_id       varchar(40) NOT NULL,
  version       int(10) NOT NULL,
  embl_version  int(10) NOT NULL,
  htg_phase     int(10) DEFAULT '-1' NOT NULL,
  created       datetime NOT NULL,
  modified      datetime NOT NULL,
  stored        datetime NOT NULL,
  
  PRIMARY KEY (internal_id),
  KEY embl (embl_id,embl_version),
  KEY id   (id,embl_version)
);


#
# Table structure for table 'map_density'
#
CREATE TABLE map_density (
   chr_name   varchar(15) NOT NULL,
   chr_start	int(10) NOT NULL,
   chr_end	int(10) NOT NULL,
   type		varchar(20) NOT NULL,
   value	int(10) NOT NULL,
    
   PRIMARY KEY(type,chr_name,chr_start) 
);

#
# Table structure for table 'contig'
#
CREATE TABLE contig (
  internal_id       int(10) unsigned NOT NULL auto_increment,
  id                varchar(40) NOT NULL,
  clone             int(10) unsigned NOT NULL, # foreign key clone:internal_id
  length            int(10) unsigned NOT NULL, 
  offset            int(10) unsigned,
  corder            int(10) unsigned,
  dna               int(10) unsigned NOT NULL, # foreign key dna:id
  chromosomeId      int(10) unsigned NOT NULL,
  international_id  varchar(40),
  
  PRIMARY KEY (internal_id),
  UNIQUE id (id),
  KEY clone (clone),
  KEY dna (dna)
);



#
# Table structure for table 'dna'
#

# This table holds the sequence of the contigs from the contig table.
# The sequence is that of the contig, not that of the golden
# path. I.e. to construct the golden path from the dna entries,
# the sequence of contigs with an orientation of -1 must be
# reversed and bases complemented. The static_golden_path
# table has the contig orientation (raw_ori).
# Note the length of the dna.sequence field is always equal
# to the appropriate length field in the contig table
# (probably a violation of some form of normal form since contig.length
# is an attibute of the dna.sequence field)

CREATE TABLE dna (
  id        int(10) unsigned NOT NULL auto_increment,
  sequence  mediumtext NOT NULL,
  created   datetime NOT NULL,
  
  PRIMARY KEY (id)
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
  seq_start     int(10) NOT NULL,                     # start of exon within contig
  seq_end       int(10) NOT NULL,                     # end of exon within specified contig
  strand        tinyint(2) NOT NULL,                  # 1 or -1 depending on the strand of the exon

  phase         tinyint(2) NOT NULL,
  end_phase     tinyint(2) NOT NULL,
  sticky_rank   int(10) DEFAULT '1' NOT NULL,         # see note above
  
  PRIMARY KEY ( exon_id, sticky_rank),
  KEY contig (contig_id)
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

#
# Table structure for table 'feature'
#
CREATE TABLE feature (
  id            int(10) unsigned NOT NULL auto_increment,
  contig        int(10) unsigned NOT NULL,      # foreign key contig:internal_id  
  seq_start     int(10) NOT NULL,
  seq_end       int(10) NOT NULL,
  score         double(16,4) NOT NULL,
  strand        int(1) DEFAULT '1' NOT NULL,
  analysis      int(10) unsigned NOT NULL,      # foreign key analysisprocess:analysisId
  name          varchar(40),
  hstart        int(11) NOT NULL,
  hend          int(11) NOT NULL,
  hid           varchar(40) NOT NULL,
  evalue        varchar(20),
  perc_id       tinyint(10),
  phase         tinyint(1),
  end_phase     tinyint(1),
  
  PRIMARY KEY (id),
  KEY contig_ana_score (contig, analysis, score),
  KEY hid (hid),
  KEY ana (analysis)
) MAX_ROWS = 300000000 AVG_ROW_LENGTH = 80;

#
# Table structure for table 'fset'
#

# NOT USED
# The fset and fset_feature table are not actively used and will
# be dropped in the future.

CREATE TABLE fset (
  id        int(10) unsigned NOT NULL auto_increment,
  score     double(16,4) DEFAULT '0.0000' NOT NULL,
  
  PRIMARY KEY (id)
);

#
# Table structure for table 'fset_feature'
#
# NOT USED

CREATE TABLE fset_feature (
  feature   int(10) unsigned NOT NULL,
  fset      int(10) unsigned NOT NULL,
  rank      int(11) NOT NULL,
  
  PRIMARY KEY (feature,fset,rank),
  KEY fset (fset)
);

#
# Table structure for table 'gene'
#
CREATE TABLE gene (
  gene_id   int unsigned not null auto_increment,
  type VARCHAR(40) not null,
  analysisId int unsigned not null,
     
  PRIMARY KEY (gene_id)
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


#
# Table structure for table 'repeat_feature'
#
CREATE TABLE repeat_feature (
  id        int(10) unsigned NOT NULL auto_increment,
  contig    int(10) unsigned NOT NULL,         # foreign key contig:internal_id
  seq_start int(10) NOT NULL,
  seq_end   int(10) NOT NULL,
  score     double(16,4) NOT NULL,
  strand    tinyint(1) DEFAULT '1' NOT NULL,
  analysis  int(10) unsigned NOT NULL,         # foregin key analysisprocess:analysisId
  hstart    int(11) NOT NULL,
  hend      int(11) NOT NULL,
  hid       varchar(40) NOT NULL,
  
  PRIMARY KEY (id),
  KEY contig (contig),
  KEY hid (hid)
);


#
# Table structure for table 'supporting_feature'
#
CREATE TABLE supporting_feature (
  supporting_feature_id            int(10) unsigned NOT NULL auto_increment,
  exon_id          int unsigned NOT NULL,             # foreign key exon:exon_id
  contig_id     int(10) unsigned NOT NULL,
  seq_start     int(10) NOT NULL,
  seq_end       int(10) NOT NULL,
  score         int(10) NOT NULL,
  strand        int(1) DEFAULT '1' NOT NULL,
  analysis      int(10) unsigned NOT NULL,
  name          varchar(40) NOT NULL,
  hstart        int(11) NOT NULL,
  hend          int(11) NOT NULL,
  hid           varchar(40) NOT NULL,
  evalue        double(16,4),
  perc_id       int(10),
  phase         tinyint(1),
  end_phase     tinyint(1),
  hstrand       tinyint(1),
  
  PRIMARY KEY (supporting_feature_id),
  KEY exon (exon_id),
  KEY analysis (analysis),
  KEY hid (hid),
  KEY name (name)
);

#
# Table structure for table 'transcript'
#
CREATE TABLE transcript (
  transcript_id    INT UNSIGNED NOT NULL auto_increment,  
  gene_id          INT UNSIGNED NOT NULL,          # foreign key gene:gene_id
  translation_id   INT UNSIGNED NOT NULL,          # foreign key translation:translation_id
  
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
  seq_start       INT(10) NOT NULL,
  start_exon_id   INT UNSIGNED NOT NULL,  # foreign key exon:exon_id
  seq_end         INT(10) NOT NULL,
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
# the absolute position within the golden path chromosome (chr_name) is given
# by chr_start and chr_end. Each contig is in some "fingerprint clone contig"
# and the fpc contig is identified by field fpcctg_name and the position of
# the specified bit of the contig within its FPC contig is given by fields
# fpcctg_start and fpcctg_end.
# With the data set at time of this writing, field type is always "UCSC".
# 
# NB, chr_start <= chr_end, raw_start <= raw_end, and fpcctg_start <= fpcctg_end.
# 
 

CREATE TABLE static_golden_path (
    fpcctg_name    varchar(20) NOT NULL,
    chr_name       varchar(20)  NOT NULL,
    raw_id         int(10) unsigned NOT NULL, # foreign key contig:internal_id
    chr_start      int(10) NOT NULL,
    chr_end        int(10) NOT NULL,
    fpcctg_start   int(10) NOT NULL,
    fpcctg_end     int(10) NOT NULL,
    raw_start      int(10) NOT NULL,
    raw_end        int(10) NOT NULL,
    raw_ori        tinyint(2)  NOT NULL, 
    type           varchar(20) NOT NULL,
    
    PRIMARY KEY(raw_id,type),
    KEY(fpcctg_name, fpcctg_start),
    KEY(chr_name,chr_start) 
);



#
# Table structure for table 'protein_feature'
#

CREATE TABLE protein_feature (
  id            int(10) unsigned NOT NULL auto_increment,
  translation   varchar(40) NOT NULL,	
  seq_start     int(10) NOT NULL,
  seq_end       int(10) NOT NULL,
  analysis      int(10) unsigned NOT NULL,
  hstart        int(10) NOT NULL,
  hend          int(10) NOT NULL,
  hid           varchar(40) NOT NULL,
  score         double(16,4) NOT NULL,
  evalue        varchar(20),
  perc_id       int(10),

  PRIMARY KEY   (id),
  KEY (translation),
  KEY hid_index ( hid )
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

CREATE TABLE interpro_description (
  interpro_ac varchar(40) DEFAULT '' NOT NULL,
  description varchar(255),
  short_description varchar(255),
  PRIMARY KEY (interpro_ac)
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
   chr_name  varchar(40) NOT NULL,
   chr_start int(10)     NOT NULL,
   chr_end   int(10)     NOT NULL,
   band      varchar(40) NOT NULL,
   stain     varchar(40) NOT NULL,
   PRIMARY KEY (chr_name,band)
);


#
#Table structure for table objectXref
#

CREATE TABLE objectXref(
       objectxrefId INT unsigned not null auto_increment,
       ensembl_id VARCHAR(40) not null, 
       ensembl_object_type ENUM( 'RawContig', 'Transcript', 'Gene', 'Translation' ) not null,
       xrefId INT unsigned not null,

       UNIQUE ( ensembl_object_type, ensembl_id, xrefId ),
       PRIMARY KEY (objectxrefId)
   	);

#
#Table structure for identityXref
#
CREATE TABLE identityXref(
        objectxrefId INT unsigned not null ,
	query_identity 	int(5),
        target_identity int(5),
        PRIMARY KEY (objectxrefId)
        );


#
#Table structure for table Xref
#

CREATE TABLE Xref(
         xrefId INT unsigned not null auto_increment,
         externalDBId int unsigned not null,
         dbprimary_id VARCHAR(40) not null,
	 display_id VARCHAR(40) not null,
         version VARCHAR(10) DEFAULT '' NOT NULL,
	 description VARCHAR(255),

         PRIMARY KEY( xrefId ),
         UNIQUE KEY id_index( dbprimary_id, externalDBId ),
         KEY display_index ( display_id )

   	);


#
#Table structure for table externalSynonym
#

CREATE TABLE externalSynonym(
         xrefId INT unsigned not null,
         synonym VARCHAR(40) not null,
         PRIMARY KEY( xrefId, synonym ),
	 KEY name_index( synonym )
   	);

#
#Table structure for table externalDB 
#

CREATE TABLE externalDB(
         externalDBId INT unsigned not null auto_increment,
         db_name VARCHAR(40) not null,
	 release VARCHAR(40) DEFAULT '' NOT NULL,
         PRIMARY KEY( externalDBId ) 
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
  chr_name char(20)  NOT NULL,
  KEY chr_name (chr_name,chr_start)
);





CREATE TABLE meta (
    meta_id INT unsigned not null auto_increment,
    meta_key varchar( 40 ) not null,
    meta_value varchar( 255 ) not null,

    PRIMARY KEY( meta_id ),
    KEY meta_key_index ( meta_key ),
    KEY meta_value_index ( meta_value )
);

# Auto add schema version to database
insert into meta (meta_key, meta_value) values (
	"schema_version", "$Revision$");
