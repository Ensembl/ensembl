# MySQL dump 6.4
#
# Host: localhost    Database: archive
#--------------------------------------------------------
# Server version	3.22.27


#
# Table structure for table 'sequence'
#
CREATE TABLE sequence (
  id varchar(40) DEFAULT '0' NOT NULL,	
  version varchar(5) DEFAULT '' NOT NULL,
  seq_type set("transcript","protein","exon") DEFAULT '' NOT NULL,   
  gene_id varchar(40) DEFAULT '' NOT NULL,  
  gene_version varchar(5) DEFAULT '' NOT NULL,
  clone_id varchar(40) DEFAULT '' NOT NULL,
  clone_version varchar(5) DEFAULT '' NOT NULL,
  sequence mediumtext NOT NULL,
  PRIMARY KEY (id,version,seq_type)
);

CREATE TABLE deleted_id (
  old_id varchar(40) DEFAULT '' NOT NULL,
  old_version int(10) DEFAULT '' NOT NULL,
  new_id varchar(40), 
  id_type varchar(40) DEFAULT '' NOT NULL, 
  PRIMARY KEY (old_id,id_type)
);
