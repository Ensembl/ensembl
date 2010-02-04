# patch_56_57_b.sql
#
# title: unmapped_object.type enum tidy
#
# description:
# Change unmapped_object.type enum to remove probe2transcript

ALTER table unmapped_object modify type ENUM('xref', 'cDNA', 'Marker') NOT NULL;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_56_57_b.sql|unmapped_object.typ_enum_tidy');


