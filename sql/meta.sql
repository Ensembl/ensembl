# MySQL dump 6.4
#
# Host: localhost    Database: eliatest1
#--------------------------------------------------------
# Server version	3.22.27


#
# Table structure for table 'meta'
#
CREATE TABLE meta (
  last_update datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
  donor_database_locator varchar(40) DEFAULT '' NOT NULL,
  offset_time time(5) DEFAULT '30' NOT NULL,
  schema_version varchar(40) DEFAULT '' NOT NULL
);
