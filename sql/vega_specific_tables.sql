
# Script to create Vega specific database tables - this is used in 
# conjunction with table.sql (creates Ensembl Schema) to create a 
# Vega new schema database.
#
# Conventions are the same as for Ensmebl:
#  - use lower case and underscores
#  - internal ids are integers named tablename_id
#  - same name is given in foreign key relations
#
# Steve Trevanion (st3@sanger.ac.uk) created 2/7/04


################################################################################
#
# Table structure for table 'gene_synonym'
#

CREATE TABLE gene_synonym (
  synonym_id int(10) unsigned NOT NULL auto_increment,
  name varchar(100) default NULL,
  gene_info_id int(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (synonym_id)
) TYPE=MyISAM;

################################################################################
#
# Table structure for table 'translation_stable_id_pool'
#

CREATE TABLE translation_stable_id_pool (
  translation_pool_id int(10) unsigned NOT NULL auto_increment,
  translation_stable_id varchar(40) default NULL,
  timestamp datetime NOT NULL default '0000-00-00 00:00:00',
  PRIMARY KEY  (translation_pool_id)
) TYPE=MyISAM;

################################################################################
#
# Table structure for table 'transcript_info'
#

CREATE TABLE transcript_info (
  transcript_info_id int(10) unsigned NOT NULL auto_increment,
  transcript_stable_id varchar(40) default NULL,
  name varchar(40) default NULL,
  transcript_class_id int(10) unsigned default NULL,
  cds_start_not_found enum('true','false') NOT NULL default 'true',
  cds_end_not_found enum('true','false') NOT NULL default 'true',
  mRNA_start_not_found enum('true','false') NOT NULL default 'true',
  mRNA_end_not_found enum('true','false') NOT NULL default 'true',
  author_id int(10) unsigned NOT NULL default '0',
  timestamp datetime NOT NULL default '0000-00-00 00:00:00',
  PRIMARY KEY  (transcript_info_id)
) TYPE=MyISAM;

################################################################################
#
# Table structure for table 'current_gene_info'
#

CREATE TABLE current_gene_info (
  gene_info_id int(10) unsigned NOT NULL default '0',
  gene_stable_id varchar(40) default NULL,
  PRIMARY KEY  (gene_info_id),
  UNIQUE KEY gene_stable_id (gene_stable_id)
) TYPE=MyISAM;

################################################################################
#
# Table structure for table 'current_transcript_info'
#

CREATE TABLE current_transcript_info (
  transcript_info_id int(10) unsigned NOT NULL default '0',
  transcript_stable_id varchar(40) default NULL,
  PRIMARY KEY  (transcript_info_id),
  UNIQUE KEY transcript_stable_id (transcript_stable_id)
) TYPE=MyISAM;

################################################################################
#
# Table structure for table 'author'
#

CREATE TABLE author (
  author_id int(10) unsigned NOT NULL auto_increment,
  author_email varchar(50) default NULL,
  author_name varchar(50) default NULL,
  PRIMARY KEY  (author_id)
) TYPE=MyISAM;

################################################################################
#
# Table structure for table 'gene_name'
#

CREATE TABLE gene_name (
  gene_name_id int(10) unsigned NOT NULL auto_increment,
  name varchar(100) default NULL,
  gene_info_id int(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (gene_name_id)
) TYPE=MyISAM;


################################################################################
#
# Table structure for table 'transcript_class'
#

CREATE TABLE transcript_class (
  transcript_class_id int(10) unsigned NOT NULL auto_increment,
  name varchar(40) default NULL,
  description varchar(255) default NULL,
  PRIMARY KEY  (transcript_class_id),
  UNIQUE KEY name (name)
) TYPE=MyISAM;

################################################################################
#
# Table structure for table 'gene_remark'
#

CREATE TABLE gene_remark (
  gene_remark_id int(10) unsigned NOT NULL auto_increment,
  remark text,
  gene_info_id int(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (gene_remark_id)
) TYPE=MyISAM;


################################################################################
#
# Table structure for table 'gene_info'
#

CREATE TABLE gene_info (
  gene_info_id int(10) unsigned NOT NULL auto_increment,
  gene_stable_id varchar(40) default NULL,
  author_id int(10) unsigned NOT NULL default '0',
  is_known enum('true','false') default 'false',
  timestamp datetime NOT NULL default '0000-00-00 00:00:00',
  PRIMARY KEY  (gene_info_id)
) TYPE=MyISAM;


################################################################################
#
# Table structure for table 'evidence'
#

CREATE TABLE evidence (
  evidence_id int(10) unsigned NOT NULL auto_increment,
  evidence_name varchar(40) default NULL,
  transcript_info_id int(10) unsigned default NULL,
  type enum('EST','cDNA','Protein','Genomic','UNKNOWN') default NULL,
  PRIMARY KEY  (evidence_id)
) TYPE=MyISAM;

################################################################################
#
# Table structure for table 'transcript_remark'
#

CREATE TABLE transcript_remark (
  transcript_remark_id int(10) unsigned NOT NULL auto_increment,
  remark text,
  transcript_info_id int(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (transcript_remark_id)
) TYPE=MyISAM;

################################################################################
#
# Table structure for table 'clone_info'


CREATE TABLE clone_info (
  clone_info_id int(10) unsigned NOT NULL auto_increment,
  clone_id int(10) unsigned NOT NULL default '0',
  author_id int(10) default NULL,
  timestamp datetime NOT NULL default '0000-00-00 00:00:00',
  PRIMARY KEY  (clone_info_id),
  UNIQUE clone_id_idx (clone_id )	
) TYPE=MyISAM;



################################################################################
#
# Table structure for table 'input_id_analysis'
#

CREATE TABLE input_id_analysis (
  input_id varchar(40) NOT NULL default '',
  analysis_id smallint(10) unsigned NOT NULL default '0',
  created datetime NOT NULL default '0000-00-00 00:00:00',
  result smallint(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (analysis_id,input_id),
  KEY input_created (input_id,created)
) TYPE=MyISAM;

################################################################################
#
# Table structure for table 'rule_conditions'
#

CREATE TABLE rule_conditions (
  rule_id smallint(10) unsigned NOT NULL default '0',
  condition varchar(40) default NULL
) TYPE=MyISAM;


################################################################################
#
# Table structure for table 'rule_goal'
#

CREATE TABLE rule_goal (
  rule_id smallint(5) unsigned NOT NULL auto_increment,
  goal varchar(40) default NULL,
  PRIMARY KEY  (rule_id)
) TYPE=MyISAM;

################################################################################
#
# Table structure for table 'clone_remark'


CREATE TABLE clone_remark (
  clone_remark_id int(10) unsigned NOT NULL auto_increment,
  remark text,
  clone_info_id int(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (clone_remark_id)
) TYPE=MyISAM;

################################################################################
#
# Table structure for table 'clone_info_keyword'


CREATE TABLE clone_info_keyword (
  keyword_id int(10) unsigned NOT NULL default '0',
  clone_info_id int(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (keyword_id,clone_info_id),
  KEY clone_info_id_idx (clone_info_id)
) TYPE=MyISAM;



## The following tables are  probably not used, however create them anyway to be on the safe side.

################################################################################
#
# Table structure for table 'clone_lock'

CREATE TABLE clone_lock (
  clone_lock_id int(10) unsigned NOT NULL auto_increment,
  clone_id varchar(40) NOT NULL default '',
  clone_version int(10) unsigned NOT NULL default '0',
  timestamp datetime NOT NULL default '0000-00-00 00:00:00',
  author_id int(10) NOT NULL default '0',
  PRIMARY KEY  (clone_lock_id),
  UNIQUE KEY clone_index (clone_id,clone_version)
) TYPE=MyISAM;

################################################################################
#
# Table structure for table 'keyword'

CREATE TABLE keyword (
  keyword_id int(10) unsigned NOT NULL auto_increment,
  keyword_name varchar(50) default NULL,
  PRIMARY KEY  (keyword_id)
) TYPE=MyISAM;

################################################################################
#
# Table structure for table 'current_clone_info'


CREATE TABLE current_clone_info (
  clone_id int(10) unsigned NOT NULL default '0',
  clone_info_id int(10) unsigned NOT NULL default '0',
  clone_version int(10) default NULL,
  PRIMARY KEY  (clone_id)
) TYPE=MyISAM;


################################################################################
#
# Table structure for table 'job'


CREATE TABLE job (
  job_id int(10) unsigned NOT NULL auto_increment,
  input_id varchar(40) NOT NULL default '',
  analysis_id smallint(5) unsigned NOT NULL default '0',
  submission_id mediumint(10) unsigned NOT NULL default '0',
  stdout_file varchar(100) NOT NULL default '',
  stderr_file varchar(100) NOT NULL default '',
  retry_count tinyint(2) unsigned default '0',
  PRIMARY KEY  (job_id),
  KEY input_id (input_id),
  KEY analysis_id (analysis_id)
) TYPE=MyISAM;

################################################################################
#
# Table structure for table 'job_status'

CREATE TABLE job_status (
  job_id int(10) unsigned NOT NULL default '0',
  status varchar(40) NOT NULL default 'CREATED',
  time datetime NOT NULL default '0000-00-00 00:00:00',
  is_current enum('n','y') default 'n',
  KEY job_id (job_id),
  KEY status (status),
  KEY is_current (is_current)
) TYPE=MyISAM;

 
