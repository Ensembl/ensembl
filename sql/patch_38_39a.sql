# Add a new status to the enumeration in gene & transcript
ALTER TABLE gene CHANGE COLUMN status status enum('KNOWN','NOVEL','PUTATIVE','PREDICTED','KNOWN_BY_PROJECTION');
ALTER TABLE transcript CHANGE COLUMN status status enum('KNOWN','NOVEL','PUTATIVE','PREDICTED','KNOWN_BY_PROJECTION');