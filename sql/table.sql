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
  db varchar(40),
  db_version varchar(40),
  db_file varchar(80),
  program varchar(80),
  program_version varchar(40),
  program_file varchar(40),
  parameters varchar(80),
  module varchar(80),
  module_version varchar(40),
  gff_source varchar(40),
  gff_feature varchar(40),

  PRIMARY KEY (analysisId)
);

#
# Table structure for table 'analysis'
#
CREATE TABLE analysis (
  id                int(10) unsigned NOT NULL auto_increment,
  db                varchar(40),
  db_version        varchar(5),
  program           varchar(40) NOT NULL,
  program_version   varchar(5),
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
  species_id        int(11) NOT NULL,
  id                int(11) NOT NULL,
  known_genes       int(11) NULL,
  unknown_genes     int(11) NULL,
  snps              int(11) NULL,
  length            int(11) NULL,
  
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
   chromosome_id   int(10) NOT NULL,
   chr_start	int(10) NOT NULL,
   chr_end	int(10) NOT NULL,
   type		varchar(20) NOT NULL,
   value	int(10) NOT NULL,
    
   PRIMARY KEY(type,chromosome_id,chr_start) 
);

#
# Table structure for table 'contig'
#
CREATE TABLE contig (
  internal_id       int(10) unsigned NOT NULL auto_increment,
  id                varchar(40) NOT NULL,
  clone             int(10) NOT NULL,
  length            int(10) unsigned NOT NULL,
  offset            int(10) unsigned,
  corder            int(10) unsigned,
  dna               int(10) NOT NULL,
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
CREATE TABLE dna (
  id        int(10) unsigned NOT NULL auto_increment,
  sequence  mediumtext NOT NULL,
  created   datetime NOT NULL,
  
  PRIMARY KEY (id)
) MAX_ROWS = 500000 AVG_ROW_LENGTH=100000;

#
# Table structure for table 'exon'
#
CREATE TABLE exon (
  id            varchar(40) NOT NULL,
  contig        int(10) unsigned NOT NULL,
  version       int(10) DEFAULT '1' NOT NULL,
  created       datetime NOT NULL,
  modified      datetime NOT NULL,
  stored        datetime NOT NULL,
  seq_start     int(10) NOT NULL,
  seq_end       int(10) NOT NULL,
  strand        int(2) NOT NULL,
  phase         int(11) NOT NULL,
  end_phase     int(11) NOT NULL,
  sticky_rank   int(10) DEFAULT '1' NOT NULL,
  
  PRIMARY KEY (id,sticky_rank),
  KEY contig (contig)
);

#
# Table structure for table 'exon_transcript'
#
CREATE TABLE exon_transcript (
  exon          varchar(40) NOT NULL,
  transcript    varchar(40) NOT NULL,
  rank          int(10) NOT NULL,
  
  PRIMARY KEY (exon,transcript,rank),
  KEY transcript (transcript)
);

#
# Table structure for table 'feature'
#
CREATE TABLE feature (
  id            int(10) unsigned NOT NULL auto_increment,
  contig        int(10) unsigned NOT NULL,
  seq_start     int(10) NOT NULL,
  seq_end       int(10) NOT NULL,
  score         double(16,4) NOT NULL,
  strand        int(1) DEFAULT '1' NOT NULL,
  analysis      int(10) unsigned NOT NULL,
  name          varchar(40),
  hstart        int(11) NOT NULL,
  hend          int(11) NOT NULL,
  hid           varchar(40) NOT NULL,
  evalue        varchar(20),
  perc_id       tinyint(10),
  phase         tinyint(1),
  end_phase     tinyint(1),
  
  PRIMARY KEY (id),
  KEY contig (contig),
  KEY hid (hid)
) MAX_ROWS = 100000000 AVG_ROW_LENGTH=150;

#
# Table structure for table 'fset'
#
CREATE TABLE fset (
  id        int(10) unsigned NOT NULL auto_increment,
  score     double(16,4) DEFAULT '0.0000' NOT NULL,
  
  PRIMARY KEY (id)
);

#
# Table structure for table 'fset_feature'
#
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
  id        varchar(40) NOT NULL,
  version   int(10) DEFAULT '1' NOT NULL,
  created   datetime NOT NULL,
  modified  datetime NOT NULL,
  stored    datetime NOT NULL,
  analysisId int,
     
  PRIMARY KEY (id)
);



#
# Table structure for table 'repeat_feature'
#
CREATE TABLE repeat_feature (
  id        int(10) unsigned NOT NULL auto_increment,
  contig    int(10) unsigned NOT NULL,
  seq_start int(10) NOT NULL,
  seq_end   int(10) NOT NULL,
  score     double(16,4) NOT NULL,
  strand    tinyint(1) DEFAULT '1' NOT NULL,
  analysis  int(10) unsigned NOT NULL,
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
  id            int(10) unsigned NOT NULL auto_increment,
  exon          varchar(40) NOT NULL,
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
  
  PRIMARY KEY (id),
  KEY id_exon (id,exon),
  KEY exon (exon),
  KEY analysis (analysis),
  KEY hid (hid),
  KEY name (name)
);

#
# Table structure for table 'transcript'
#
CREATE TABLE transcript (
  id            varchar(40) NOT NULL,
  version       int(10) DEFAULT '1' NOT NULL,
  gene          varchar(40) NOT NULL,
  translation   varchar(40) NOT NULL,
  
  PRIMARY KEY (id),
  KEY gene_index (gene),
  KEY translation_index ( translation )		
);

#
# Table structure for table 'translation'
#
CREATE TABLE translation (
  id            varchar(40) NOT NULL,
  version       int(10) DEFAULT '1' NOT NULL,
  seq_start     int(10) NOT NULL,
  start_exon    varchar(40) NOT NULL,
  seq_end       int(10) NOT NULL,
  end_exon      varchar(40) NOT NULL,
  
  PRIMARY KEY (id)
);


CREATE TABLE genetype (
   gene_id      varchar(40) NOT NULL,
   type  varchar(40) NOT NULL,
      
   PRIMARY KEY(gene_id),
   KEY(gene_id),
   KEY(type)
);

# this is a denormalised golden path

CREATE TABLE static_golden_path (
    fpcctg_name    varchar(20) NOT NULL,
    chr_name       varchar(20)  NOT NULL,
    raw_id         int(10) NOT NULL,
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
  evalue        double(16,4),
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
  gene_id varchar(40) NOT NULL,
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
       ensembl_id VARCHAR(40) not null, 
       ensembl_object_type ENUM( 'RawContig', 'Transcript', 'Gene', 'Translation' ) not null,
       xrefId INT not null,

       PRIMARY KEY( ensembl_object_type, ensembl_id, xrefId ),
       KEY xref_index( xrefId, ensembl_object_type, ensembl_id )
   	);			

#
#Table structure for table Xref
#

CREATE TABLE Xref(
         xrefId INT not null auto_increment,
         externalDBId int not null,
         dbprimary_id VARCHAR(40) not null,
	 display_id VARCHAR(40) not null,
         version VARCHAR(10),
	 description VARCHAR(255),

         PRIMARY KEY( xrefId ),
         KEY id_index( dbprimary_id ),
         KEY display_index ( display_id )

   	);


#
#Table structure for table externalSynonym
#

CREATE TABLE externalSynonym(
         xrefId INT not null,
         synonym VARCHAR(40) not null,
         PRIMARY KEY( xrefId, synonym ),
	 KEY name_index( synonym )
   	);

#
#Table structure for table externalDB 
#

CREATE TABLE externalDB(
         externalDBId INT not null auto_increment,
         db_name VARCHAR(40) not null,
	 release VARCHAR(40),
	 url_pattern varchar(255),
         PRIMARY KEY( externalDBId ) 
   	);




#
#Table structure for table contig_landmarkMarker
#

CREATE TABLE contig_landmarkMarker (
  contig int(10) unsigned DEFAULT '0' NOT NULL,
  marker char(40) DEFAULT '' NOT NULL,
  name char(40) DEFAULT '' NOT NULL,
  start bigint(17) DEFAULT '0' NOT NULL,
  end bigint(17) DEFAULT '0' NOT NULL,
  strand bigint(1) DEFAULT '0' NOT NULL,
  chr_name char(20) DEFAULT '' NOT NULL,
  KEY chr (chr_name),
  KEY contig (contig)
);


CREATE TABLE meta (
        meta_id INT not null auto_increment,
        meta_key varchar( 40 ) not null,
        meta_value varchar( 255 ) not null,

        PRIMARY KEY( meta_id ),
        KEY meta_key_index ( meta_key ),
        KEY meta_value_index ( meta_value )
	);









