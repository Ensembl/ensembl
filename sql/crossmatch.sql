
# this is symmetric feature pairs for raw contigs


CREATE TABLE symmetric_contig_pair_hit (
  symchid       int(10) unsigned NOT NULL auto_increment,
  score		int(10),
  perc_subs            int(10),
  perc_ins          int(10),
  perc_del          int(10),
  PRIMARY KEY(symchid)
);

CREATE TABLE symmetric_contig_feature (
  symcfid       int(10) unsigned NOT NULL auto_increment,
  symchid       int(10) NOT NULL,
  rawcontigid   varchar(40) NOT NULL,
  rawversion    int(10) NOT NULL,
  clone         varchar(40) NOT NULL,
  seq_start     int(10),
  seq_end       int(10),
  strand        int(2),
  PRIMARY KEY(symcfid),
  KEY(clone,rawcontigid),
  KEY(rawcontigid,rawversion,symchid),
  KEY(symchid)
);

CREATE TABLE dblocation (
  olddatabase  varchar(40) NOT NULL,
  newdatabase  varchar(40) NOT NULL
  );

CREATE TABLE clonelist (
  clone varchar(40) NOT NULL,
  PRIMARY KEY(clone)
);
