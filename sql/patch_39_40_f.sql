-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2018] EMBL-European Bioinformatics Institute
-- 
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
-- 
--      http://www.apache.org/licenses/LICENSE-2.0
-- 
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

# patch_39_40_f
#
# title: remove all_latest
#
# description:
# This patch removes the ALL/LATEST mapping session, and any associated stable_id_events.
# It is a data-only patch and does not involve a schema change.

DELETE s, m FROM stable_id_event s, mapping_session m WHERE m.mapping_session_id=s.mapping_session_id AND m.old_db_name='ALL' AND m.new_db_name='LATEST';

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_39_40_f.sql|remove_all_latest');

