# MySQL dump 5.13
#
# Host: obi-wan    Database: ens500
#--------------------------------------------------------
# Server version	3.22.32

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
  
  PRIMARY KEY (chromosome_id)
);

#
# Table structure for table 'clone'
#
CREATE TABLE clone (
  id            varchar(40) NOT NULL,
  embl_id       varchar(40) NOT NULL,
  version       int(10) NOT NULL,
  embl_version  int(10) NOT NULL,
  htg_phase     int(10) DEFAULT '-1' NOT NULL,
  created       datetime NOT NULL,
  modified      datetime NOT NULL,
  stored        datetime NOT NULL,
  
  PRIMARY KEY (id),
  UNIQUE embl (embl_id,embl_version)
);

#
# Table structure for table 'contig'
#
CREATE TABLE contig (
  internal_id       int(10) unsigned NOT NULL auto_increment,
  id                varchar(40) NOT NULL,
  clone             varchar(40) NOT NULL,
  length            int(10) unsigned,
  offset            int(10) unsigned,
  corder            int(10) unsigned,
  dna               int(10) NOT NULL,
  chromosomeId      int(10) unsigned NOT NULL,
  international_id  varchar(40),
  
  PRIMARY KEY (internal_id),
  KEY id (id),
  KEY clone (clone),
  KEY dna (dna)
);

#
# Table structure for table 'contigoverlap'
#
CREATE TABLE contigoverlap (
  dna_a_id              int(10) unsigned NOT NULL,
  dna_b_id              int(10) unsigned  NOT NULL,
  type                  varchar(40) NOT NULL,
  contig_a_position     int(10) unsigned,
  contig_b_position     int(10) unsigned,
  overlap_size          int(10) unsigned,
  overlap_type          enum('right2left','left2right','left2left','right2right'),
  
  PRIMARY KEY (dna_a_id,dna_b_id,type),
  KEY (type)
);

#
# Table structure for table 'db_update'
#
CREATE TABLE db_update (
  id                    int(10) NOT NULL auto_increment,
  time_started          datetime NOT NULL,
  time_finished         datetime NOT NULL,
  clones                int(10) NOT NULL,
  contigs               int(10) NOT NULL,
  genes                 int(10) NOT NULL,
  exons                 int(10) NOT NULL,
  basepairs             int(10) NOT NULL,
  features              int(10) NOT NULL,
  transcripts           int(10) NOT NULL,
  repeats               int(10) NOT NULL,
  supporting_features   int(10) NOT NULL,
  fsets                 int(10) NOT NULL,
  status                varchar(40) DEFAULT 'STARTED' NOT NULL,
  modified_start        datetime NOT NULL,
  modified_end          datetime NOT NULL,
  
  PRIMARY KEY (id)
);

#
# Table structure for table 'dna'
#
CREATE TABLE dna (
  id        int(10) unsigned NOT NULL auto_increment,
  sequence  mediumtext NOT NULL,
  created   datetime NOT NULL,
  
  PRIMARY KEY (id)
);

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
  rank          int(10) DEFAULT '1' NOT NULL,
  
  PRIMARY KEY (id,rank),
  KEY id_contig (id,contig),
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
  evalue        double(16,4),
  perc_id       int(10),
  
  PRIMARY KEY (id),
  KEY overlap (id,contig,seq_start,seq_end,analysis),
  KEY contig (contig),
  KEY hid (hid)
);

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
  
  PRIMARY KEY (id)
);

#
# Table structure for table 'geneclone_neighbourhood'
#
CREATE TABLE geneclone_neighbourhood (
  clone     varchar(40) NOT NULL,
  gene      varchar(40) NOT NULL,
  
  PRIMARY KEY (clone,gene),
  KEY gene (gene)
);

#
# Table structure for table 'ghost'
#
CREATE TABLE ghost (
  id        varchar(40) NOT NULL,
  version   varchar(5) NOT NULL,
  obj_type  enum('transcript','protein','exon') DEFAULT 'exon' NOT NULL,
  deleted   datetime NOT NULL,
  stored    datetime NOT NULL,
  
  PRIMARY KEY (id,version,obj_type)
);

#
# Table structure for table 'meta'
#
CREATE TABLE meta (
  donor_database_locator    varchar(100) NOT NULL,
  offset_time               time DEFAULT '00:30:00' NOT NULL,
  schema_version            varchar(40) NOT NULL
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
  strand    int(1) DEFAULT '1' NOT NULL,
  analysis  int(10) unsigned NOT NULL,
  hstart    int(11) NOT NULL,
  hend      int(11) NOT NULL,
  hid       varchar(40) NOT NULL,
  
  PRIMARY KEY (id),
  KEY overlap (id,contig,seq_start,seq_end,analysis),
  KEY contig (contig),
  KEY hid (hid)
);

#
# Table structure for table 'species'
#
CREATE TABLE species (
  species_id    int(10) NOT NULL auto_increment,
  nickname      varchar(40) NOT NULL,
  taxonomy_id   int(10) NOT NULL,
  
  PRIMARY KEY (species_id)
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
  
  PRIMARY KEY (id),
  KEY overlap (id,seq_start,seq_end,analysis),
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
  KEY gene_index (gene)
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


CREATE TABLE genedblink (
   gene_id      varchar(40) NOT NULL,
   external_db  varchar(40) NOT NULL,
   external_id  varchar(40) NOT NULL,
   
   PRIMARY KEY(gene_id,external_db,external_id)
);


CREATE TABLE transcriptdblink (
   transcript_id    varchar(40) NOT NULL,
   external_db      varchar(40) NOT NULL,
   external_id      varchar(40) NOT NULL,
   
   PRIMARY KEY(transcript_id,external_db,external_id)
);

