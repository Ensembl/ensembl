# MySQL dump 6.4
#
# Host: localhost    Database: eliatest1
#--------------------------------------------------------
# Server version	3.22.27


#
# Table structure for table 'ghost'
#
CREATE TABLE ghost (
  id varchar(40) DEFAULT '' NOT NULL,
  version varchar(5) DEFAULT '' NOT NULL,
  seq_type set("transcript","protein","exon") DEFAULT '' NOT NULL, 
  deleted_time datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  PRIMARY KEY (id,version,seq_type)
);
