# patch_38_39_a
#
# title: status enum
#
# description:
# this patch adds a new status to the enumeration in gene & transcript

# Add a new status to the enumeration in gene & transcript

ALTER TABLE gene CHANGE COLUMN status status enum('KNOWN','NOVEL','PUTATIVE','PREDICTED','KNOWN_BY_PROJECTION');

ALTER TABLE transcript CHANGE COLUMN status status enum('KNOWN','NOVEL','PUTATIVE','PREDICTED','KNOWN_BY_PROJECTION');

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_38_39_a.sql|status_enum');

