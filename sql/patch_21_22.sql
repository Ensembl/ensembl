
CREATE TABLE translation_attrib (
  translation_id              int(10) unsigned NOT NULL default '0',
  attrib_type_id              smallint(5) unsigned NOT NULL default '0',
  value                       varchar(255) NOT NULL default '',

  KEY type_val_idx( attrib_type_id, value ),
  KEY translation_idx( translation_id )
) TYPE=MyISAM;


CREATE TABLE transcript_attrib (
  transcript_id               int(10) unsigned NOT NULL default '0',
  attrib_type_id              smallint(5) unsigned NOT NULL default '0',
  value                       varchar(255) NOT NULL default '',

  KEY type_val_idx( attrib_type_id, value ),
  KEY transcript_idx( transcript_id )
) TYPE=MyISAM;

alter table misc_set modify code varchar(25) NOT NULL default ''

