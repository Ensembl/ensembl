# changes to database structure from 26 to 27
# none of which should affect your current data
#  so you dont need to apply them


ALTER TABLE xref MODIFY dbprimary_acc VARCHAR(40) BINARY NOT NULL;
ALTER TABLE affy_probe MODIFY probeset VARCHAR(40);
ALTER TABLE interpro DROP INDEX interpro_ac;
ALTER TABLE interpro DROP INDEX id;

ALTER TABLE interpro ADD UNIQUE ( interpro_ac, id );
ALTER TABLE interpro ADD INDEX( id );



