# MySQL dump 5.13
#
# Host: localhost    Database: bugtrack
#--------------------------------------------------------
# Server version	3.22.22

#
# Table structure for table 'bug'
#
CREATE TABLE bug (
  id int(10) DEFAULT '0' NOT NULL auto_increment,
  title varchar(255) DEFAULT '' NOT NULL,
  type set('modules','web','other') DEFAULT 'modules' NOT NULL,
  PRIMARY KEY (id)
);

#
# Table structure for table 'worknote'
#
CREATE TABLE worknote (
  bug int(10) DEFAULT '0' NOT NULL, 
  author varchar(255) DEFAULT '' NOT NULL,
  date datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  note text NOT NULL
);

#
# Table structure for table 'bug_archive'
#
CREATE TABLE bug_archive (
  id int(10) DEFAULT '0' NOT NULL auto_increment,
  title varchar(255) DEFAULT '' NOT NULL,
  type set('modules','web','other') DEFAULT 'modules' NOT NULL,
  last_author varchar(255) DEFAULT '' NOT NULL,
  last_date datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  last_note text DEFAULT '' NOT NULL,
  primary key (id)
);
