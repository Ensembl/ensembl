# MySQL dump 6.4
#
# Host: localhost    Database: ensdev
#--------------------------------------------------------
# Server version	3.22.27

#
# Table structure for table 'analysis'
#
CREATE TABLE analysis (
  db varchar(40),
  db_version varchar(5),
  program varchar(40) DEFAULT '' NOT NULL,
  program_version varchar(5),
  gff_source varchar(40),
  gff_feature varchar(40),
  id int(11) DEFAULT '0' NOT NULL auto_increment,
  PRIMARY KEY (id)
);

#
# Table structure for table 'analysis_history'
#
CREATE TABLE analysis_history (
  contig varchar(40),
  analysis varchar(40),
  created date
);

#
# Table structure for table 'clone'
#
CREATE TABLE clone (
  id varchar(40) DEFAULT '' NOT NULL,
  embl_id varchar(40) DEFAULT '' NOT NULL,
  version int(10) DEFAULT '0' NOT NULL,
  embl_version int(10) DEFAULT '0' NOT NULL,
  htg_phase int(10) DEFAULT '-1' NOT NULL,
  created datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  modified datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  stored datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,	
  PRIMARY KEY (id)
);

#
# Table structure for table 'contig'
#
CREATE TABLE contig (
  id varchar(40) DEFAULT '' NOT NULL,
  clone varchar(40) DEFAULT '' NOT NULL,
  mapbin varchar(40) DEFAULT '' NOT NULL,
  length int(10) unsigned,
  offset int(10) unsigned,
  orientation int(1) DEFAULT '1' NOT NULL,
  corder int(10) unsigned,
  dna varchar(40),
  PRIMARY KEY (id)
);

#
# Table structure for table 'contig_equiv'
#
CREATE TABLE contig_equiv (
  id varchar(40) DEFAULT '' NOT NULL,
  from_contig_id varchar(40),
  from_contig_postion int(11),
  from_contig_strand char(1),
  to_contig_id varchar(40),
  to_contig_position int(11),
  to_contig_strand char(1),
  PRIMARY KEY (id)
);

#
# Table structure for table 'dna'
#
CREATE TABLE dna (
  contig varchar(40) DEFAULT '' NOT NULL,
  sequence mediumtext NOT NULL,
  created datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  id int(10) unsigned DEFAULT '0' NOT NULL auto_increment,
  PRIMARY KEY (contig),
  KEY id (id),
  KEY idx1 (id,contig)
);

#
# Table structure for table 'exon'
#
CREATE TABLE exon (
  id varchar(40) DEFAULT '' NOT NULL,
  contig varchar(40) DEFAULT '' NOT NULL,
  version int(10) DEFAULT '1' NOT NULL,
  created datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  modified datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  stored datetime DEFAULT '0000-00-00 00:00:00' NOT NULL, 
  seq_start int(10) DEFAULT '0' NOT NULL,
  seq_end int(10) DEFAULT '0' NOT NULL,
  strand int (2) DEFAULT '1' NOT NULL,
  phase int(11) DEFAULT '0' NOT NULL,
  end_phase int(11) DEFAULT '0' NOT NULL,
  KEY idx1 (id,contig),
  PRIMARY KEY (id)
);

#
# Table structure for table 'exon_feature'
#
CREATE TABLE exon_feature (
  feature varchar(40),
  exon varchar(40)
);

#
# Table structure for table 'exon_transcript'
#
CREATE TABLE exon_transcript (
  exon varchar(40) DEFAULT '' NOT NULL,
  transcript varchar(40) DEFAULT '' NOT NULL,
  rank int(10) DEFAULT '0' NOT NULL,
  PRIMARY KEY(exon,transcript,rank),	
  KEY idx1 (exon,transcript)
);

#
# Table structure for table 'feature'
#
CREATE TABLE feature (
  id int(10) unsigned DEFAULT '0' NOT NULL auto_increment,
  contig varchar(40) DEFAULT '' NOT NULL,
  seq_start int(10) DEFAULT '0' NOT NULL,
  seq_end int(10) DEFAULT '0' NOT NULL,
  score int(10) DEFAULT '0' NOT NULL,
  strand int(1) DEFAULT 1 NOT NULL,
  analysis varchar(40) DEFAULT '' NOT NULL,
  name varchar(40),
  hstart int(11) DEFAULT '0' NOT NULL,
  hend int(11) DEFAULT '0' NOT NULL,
  hid varchar(40) DEFAULT '' NOT NULL,
  KEY overlap (id,contig,seq_start,seq_end,analysis),
  PRIMARY KEY(id)
);

# Table structure for repeat features

CREATE TABLE repeat_feature (
  id int(10) unsigned DEFAULT '0' NOT NULL auto_increment,
  contig varchar(40) DEFAULT '' NOT NULL,
  seq_start int(10) DEFAULT '0' NOT NULL,
  seq_end int(10) DEFAULT '0' NOT NULL,
  score int(10) DEFAULT '0' NOT NULL,
  strand int(1) DEFAULT 1 NOT NULL,
  analysis varchar(40) DEFAULT '' NOT NULL,
  hstart int(11) DEFAULT '0' NOT NULL,
  hend int(11) DEFAULT '0' NOT NULL,
  hid varchar(40) DEFAULT '' NOT NULL,
  KEY overlap (id,contig,seq_start,seq_end,analysis),
  PRIMARY KEY(id)
);


#
# Table structure for table 'fset'
#
CREATE TABLE fset (
  id int(10) unsigned DEFAULT '0' NOT NULL auto_increment,
  score double(16,4) DEFAULT '0.0000' NOT NULL,
  PRIMARY KEY (id)
);

#
# Table structure for table 'fset_feature'
#
CREATE TABLE fset_feature (
  fset varchar(40) DEFAULT '' NOT NULL,
  feature varchar(40) DEFAULT '' NOT NULL,
  rank int(11) DEFAULT '' NOT NULL,
  PRIMARY KEY(fset,feature,rank)
);

#
# Table structure for table 'gene'
#
CREATE TABLE gene (
  id varchar(40) DEFAULT '' NOT NULL,
  version int(10) DEFAULT '1' NOT NULL,
  created datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  modified datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  stored datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  PRIMARY KEY (id)
);

#
# Table structure for table 'ghost'
#
CREATE TABLE ghost (
  id varchar(40) DEFAULT '' NOT NULL,
  version varchar(5) DEFAULT '' NOT NULL,
  obj_type set("transcript","protein","exon") DEFAULT '' NOT NULL, 
  deleted datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  stored datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  PRIMARY KEY (id,version,obj_type)
);


#
# Table structure for table 'mapbin'
#
CREATE TABLE mapbin (
  id varchar(40) DEFAULT '' NOT NULL,
  chromosome char(2) DEFAULT '' NOT NULL,
  PRIMARY KEY (id)
);

#
# Table structure for table 'meta'
#
CREATE TABLE meta (
  last_update datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  donor_database_locator varchar(100) DEFAULT '' NOT NULL,
  offset_time time DEFAULT '00:30:00' NOT NULL,
  schema_version varchar(40) DEFAULT '' NOT NULL
);

#
# Table structure for table 'supporting_feature'
#
CREATE TABLE supporting_feature (
  id int(10) unsigned DEFAULT '0' NOT NULL auto_increment,
  exon varchar(40) DEFAULT '' NOT NULL,
  seq_start int(10) DEFAULT '0' NOT NULL,
  seq_end int(10) DEFAULT '0' NOT NULL,
  score int(10) DEFAULT '0' NOT NULL,
  strand int(1) DEFAULT 1 NOT NULL,
  analysis varchar(40) DEFAULT '' NOT NULL,
  name varchar(40),
  hstart int(11) DEFAULT '0' NOT NULL,
  hend int(11) DEFAULT '0' NOT NULL,
  hid varchar(40) DEFAULT '' NOT NULL,
  KEY overlap (id,seq_start,seq_end,analysis),
  KEY exon_id (id,exon),
  PRIMARY KEY(id)
);


#
# Table structure for table 'transcript'
#
CREATE TABLE transcript (
  id varchar(40) DEFAULT '' NOT NULL,
  version int(10) DEFAULT '1' NOT NULL,
  gene varchar(40) DEFAULT '' NOT NULL,
  translation varchar(40) DEFAULT '' NOT NULL,
  PRIMARY KEY (id),
  KEY id_geneid (id)
);

CREATE TABLE translation (
  id varchar(40) DEFAULT '' NOT NULL,
  version int(10) DEFAULT '1' NOT NULL,
  seq_start int(10) DEFAULT '0' NOT NULL,
  start_exon varchar(40) DEFAULT '' NOT NULL,
  seq_end   int(10) DEFAULT '0' NOT NULL,
  end_exon varchar(40) DEFAULT '' NOT NULL,
  PRIMARY KEY(id)
); 

CREATE TABLE geneclone_neighbourhood ( 
       clone varchar(40) NOT NULL, 
       gene varchar(40) NOT NULL, 
       primary key (clone,gene)
);










