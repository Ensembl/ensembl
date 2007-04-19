# patch_44_45_b.sql
#
# title: Marker index
#
# description: 
# Add index to marker.display_marker_synonym_id

ALTER TABLE marker ADD INDEX display_idx (display_marker_synonym_id);

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_44_45_b.sql|marker_index');


