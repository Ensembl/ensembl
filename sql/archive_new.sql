CREATE TABLE seq (
       	seq_id   int(10) unsigned NOT NULL auto_increment,
       	name varchar(40) NOT NULL,
	type enum('exon','transcript','translation','gene') NOT NULL,
	created datetime NOT NULL DEFAULT 'now()',

	PRIMARY KEY(seq_id),
     	UNIQUE (name),
	KEY (type)
);

CREATE TABLE versioned_seq (
	versioned_seq_id int(10) unsigned NOT NULL auto_increment,
	seq_id  int(10) unsigned NOT NULL,
	version int(10),
	sequence mediumtext,
	start_clone varchar(40),
	start_coord int(10),
	end_clone varchar(40),
	end_coord int(10),
	modified datetime NOT NULL DEFAULT 'now()',
	release_version int(10),

	PRIMARY KEY(versioned_seq_id),
	UNIQUE(seq_id,version),
	KEY(modified),
	KEY(release_version)
);

CREATE TABLE versioned_seq_history (
	old_versioned_seq_id int(10) unsigned NOT NULL,
	new_versioned_seq_id int(10) unsigned NOT NULL,

	KEY (old_versioned_seq_id),
	KEY (new_versioned_seq_id)
);

CREATE TABLE versioned_seq_relatives (
	master_versioned_seq_id int(10) unsigned NOT NULL,
	relative_versioned_seq_id int(10) unsigned NOT NULL,
	
	KEY (master_versioned_seq_id),
	KEY (relative_versioned_seq_id)
);

	
