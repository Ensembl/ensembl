# patch_56_57_b.sql
#
# title: Tidy old affy analyses
#
# description:
# Remove AffyAlign and probe2transcript from analysis & analysis_description table
# Change unmapped_object.type enum to remove probe2transcript

-- Delete in two parts just in case we don't have an ad entry
DELETE ad from analysis a, analysis_description ad where a.logic_name='AffyAlign' and a.analysis_id=ad.analysis_id;
DELETE ad from analysis a, analysis_description ad where a.logic_name='probe2transcript' and a.analysis_id=ad.analysis_id;

DELETE from analysis where logic_name='AffyAlign';
DELETE from analysis where logic_name='probe2transcript';

ALTER table unmapped_object modify type ENUM('xref', 'cDNA', 'Marker') NOT NULL;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_56_57_b.sql|affy_analysis_tidy');


