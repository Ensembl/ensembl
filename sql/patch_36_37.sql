# increase width of xref display_label column to allow for longer labels
ALTER TABLE xref CHANGE COLUMN display_label display_label VARCHAR(255) NOT NULL;
