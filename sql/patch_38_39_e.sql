# patch_38_39_e
#
# title: unknown status enum
#
# description:
# this patch adds a new status of UNKNOWN to the enumeration in gene & transcript

ALTER TABLE gene CHANGE COLUMN status status enum('KNOWN','NOVEL','PUTATIVE','PREDICTED','KNOWN_BY_PROJECTION', 'UNKNOWN');

ALTER TABLE transcript CHANGE COLUMN status status enum('KNOWN','NOVEL','PUTATIVE','PREDICTED','KNOWN_BY_PROJECTION', 'UNKNOWN');

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_38_39_e.sql|unknown_status_enum');

