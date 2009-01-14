# patch_52_53_d.sql
#
# title: drop go_xref index
#
# description:
# Remove the separate index on go_xref.object_xref as this is indexed in the UNIQUE constraint.

ALTER TABLE go_xref DROP INDEX object_xref_id;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_52_53_d.sql|drop_go_xref_index');


