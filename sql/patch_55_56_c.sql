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

# patch_55_56_c.sql
#
# title: oligo tables and xrefs
#
# description:
# Drop oligo tables and remove association xref schema data


-- Delete all AFFY object_xrefs
-- Do this on external_db_ids and dbname, to be safer i.e. we do not want to corrupt and custom data
DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id>=3000 and edb.external_db_id<3200 and edb.db_name like 'AFFY\_%' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;


-- And xrefs
DELETE x from external_db edb, xref x where edb.external_db_id>=3000 and edb.external_db_id<3200 and edb.db_name like 'AFFY\_%' and edb.external_db_id=x.external_db_id;

-- And external_dbs

DELETE from external_db where external_db_id>=3000 and external_db_id<3200 and db_name like 'AFFY\_%';


-- Species specific xrefs, to do this safely have to do these on a per db basis due to specific names
-- aedes_aegypti_core_55_1d
-- oligo_array_id  parent_array_id probe_setsize   name    type	external_db_id
-- 3       NULL    1       ARRAY_ND_TIGRTC_9_6K_v1 OLIGO	3207
-- 2       NULL    1       ARRAY_LIV_AEGDETOX_0_25K_v1     OLIGO	3206
-- 4       NULL    1       ARRAY_UCR_GillMgMT_0_2K_v2      OLIGO	3208

DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id=3207 and edb.db_name='ARRAY_ND_TIGRTC_9_6K_v1' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;
DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id=3206 and edb.db_name='ARRAY_LIV_AEGDETOX_0_25K_v1' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;
DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id=3208 and edb.db_name='ARRAY_UCR_GillMgMT_0_2K_v2' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;



-- And xrefs
DELETE x from external_db edb, xref x where edb.external_db_id=3207 and edb.db_name='ARRAY_ND_TIGRTC_9_6K_v1' and edb.external_db_id=x.external_db_id ;
DELETE x from external_db edb, xref x where edb.external_db_id=3206 and edb.db_name='ARRAY_LIV_AEGDETOX_0_25K_v1' and edb.external_db_id=x.external_db_id ;
DELETE x from external_db edb, xref x where edb.external_db_id=3208 and edb.db_name='ARRAY_UCR_GillMgMT_0_2K_v2' and edb.external_db_id=x.external_db_id ;


-- And external_dbs

DELETE from external_db  where external_db_id=3207 and db_name='ARRAY_ND_TIGRTC_9_6K_v1';
DELETE from external_db  where external_db_id=3206 and db_name='ARRAY_LIV_AEGDETOX_0_25K_v1';
DELETE from external_db  where external_db_id=3209 and db_name='ARRAY_UCR_GillMgMT_0_2K_v2';


-- anopheles_gambiae_core_55_3k
-- oligo_array_id  parent_array_id probe_setsize   name    type	external_db_id
-- 3       NULL    1       LIV_GAMDETOX_0.25k_v1   OLIGO	3204
-- 4       NULL    1       LIV_GAMDETOX_0.25k_v2   OLIGO	3210
-- 5       NULL    1       EMBL_MMC1_20k_v1        OLIGO	3201
-- 6       NULL    1       EMBL_MMC2_12k_v1        OLIGO	3202
-- 7       NULL    1       JHSPH_AG_GAMBER_15k_v1  OLIGO	3209


DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id=3204 and edb.db_name='ARRAY_LIV_GAMDETOX_0.25k_v1' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;
DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id=3210 and edb.db_name='ARRAY_LIV_GAMDETOX_0.25k_v2' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;
DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id=3201 and edb.db_name='ARRAY_EMBL_MMC1_20k_v1' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;
DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id=3202 and edb.db_name='ARRAY_EMBL_MMC2_12k_v1' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;
DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id=3209 and edb.db_name='ARRAY_JHSPH_AG_GAMBER_15k_v1' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;

DELETE x from external_db edb, xref x where edb.external_db_id=3204 and edb.db_name='ARRAY_LIV_GAMDETOX_0.25k_v1' and edb.external_db_id=x.external_db_id ;
DELETE x from external_db edb, xref x where edb.external_db_id=3210 and edb.db_name='ARRAY_LIV_GAMDETOX_0.25k_v2' and edb.external_db_id=x.external_db_id ;
DELETE x from external_db edb, xref x where edb.external_db_id=3201 and edb.db_name='ARRAY_EMBL_MMC1_20k_v1' and edb.external_db_id=x.external_db_id ;
DELETE x from external_db edb, xref x where edb.external_db_id=3202 and edb.db_name='ARRAY_EMBL_MMC2_12k_v1' and edb.external_db_id=x.external_db_id ;
DELETE x from external_db edb, xref x where edb.external_db_id=3209 and edb.db_name='ARRAY_JHSPH_AG_GAMBER_15k_v1' and edb.external_db_id=x.external_db_id ;


DELETE from external_db  where external_db_id=3204 and db_name='ARRAY_LIV_GAMDETOX_0.25k_v1';
DELETE from external_db  where external_db_id=3210 and db_name='ARRAY_LIV_GAMDETOX_0.25k_v2';
DELETE from external_db  where external_db_id=3201 and db_name='ARRAY_EMBL_MMC1_20k_v1';
DELETE from external_db  where external_db_id=3202 and db_name='ARRAY_EMBL_MMC2_12k_v1';
DELETE from external_db  where external_db_id=3200 and db_name='ARRAY_JHSPH_AG_GAMBER_15k_v1';


-- danio
-- oligo_array_id  parent_array_id probe_setsize   name    type	external_db_id
-- 2       NULL    1       AGILENT_G2518A  OLIGO	3251
-- 3       NULL    1       AGILENT_G2519F  OLIGO	3252
-- 4       NULL    1       LEIDEN2 OLIGO	3260
-- 5       NULL    1       LEIDEN3 OLIGO	3261

DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id=3251 and edb.db_name='AGILENT_G2518A' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;
DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id=3252 and edb.db_name='AGILENT_G2519F' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;
DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id=3260 and edb.db_name='ARRAY_ZFISH_LEIDEN2' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;
DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id=3261 and edb.db_name='ARRAY_ZFISH_LEIDEN3' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;

DELETE x from external_db edb, xref x where edb.external_db_id=3251 and edb.db_name='AGILENT_G2518A' and edb.external_db_id=x.external_db_id ;
DELETE x from external_db edb, xref x where edb.external_db_id=3252 and edb.db_name='AGILENT_G2519F' and edb.external_db_id=x.external_db_id ;
DELETE x from external_db edb, xref x where edb.external_db_id=3260 and edb.db_name='ARRAY_ZFISH_LEIDEN2' and edb.external_db_id=x.external_db_id ;
DELETE x from external_db edb, xref x where edb.external_db_id=3261 and edb.db_name='ARRAY_ZFISH_LEIDEN3' and edb.external_db_id=x.external_db_id ;

DELETE from external_db  where external_db_id=3251 and db_name='AGILENT_G2518A';
DELETE from external_db  where external_db_id=3252 and db_name='AGILENT_G2519F';
DELETE from external_db  where external_db_id=3260 and db_name='ARRAY_ZFISH_LEIDEN2';
DELETE from external_db  where external_db_id=3261 and db_name='ARRAY_ZFISH_LEIDEN3';

-- platypus
-- platypus_exon? Is done as affy for oligos, but no for xrefs.
-- Platypus_exon 3299
DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id=3299 and edb.db_name='Platypus_exon' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;
DELETE x from external_db edb, xref x where edb.external_db_id=3299 and edb.db_name='Platypus_exon' and edb.external_db_id=x.external_db_id ;
DELETE from external_db  where external_db_id=3299 and db_name='Platypus_exon';


-- Also have:
-- VB_array 3200
-- ARRAY_JHU_GAM3_0_21k_v1 3203
-- ARRAY_JHU_AEG1_0_20k_v1 3205	
-- ARRAY_LIV_AEGDETOX_0_25k_v2 3211
-- Seeking confirmation from Dan Lawson wrt removing these

DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id=3203 and edb.db_name='ARRAY_JHU_GAM3_0_21k_v1' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;
DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id=3205 and edb.db_name='ARRAY_JHU_AEG1_0_20k_v1' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;
DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id=3211 and edb.db_name='ARRAY_LIV_AEGDETOX_0_25k_v2' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;

DELETE x from external_db edb, xref x where edb.external_db_id=3203 and edb.db_name='ARRAY_JHU_GAM3_0_21k_v1' and edb.external_db_id=x.external_db_id ;
DELETE x from external_db edb, xref x where edb.external_db_id=3205 and edb.db_name='ARRAY_JHU_AEG1_0_20k_v1' and edb.external_db_id=x.external_db_id ;
DELETE x from external_db edb, xref x where edb.external_db_id=3211 and edb.db_name='ARRAY_LIV_AEGDETOX_0_25k_v2' and edb.external_db_id=x.external_db_id ;


DELETE from external_db  where external_db_id=3203 and db_name='ARRAY_JHU_GAM3_0_21k_v1';
DELETE from external_db  where external_db_id=3205 and db_name='ARRAY_JHU_AEG1_0_20k_v1';
DELETE from external_db  where external_db_id=3211 and db_name='ARRAY_LIV_AEGDETOX_0_25k_v2';


-- Also need to remove xref only array data for db_name like AGILENT, Illumina, CODELINK

--  +----------------+-------------+
--  | external_db_id | db_name     |
--  +----------------+-------------+
--  |           3250 | Illumina    | 
--  |           3255 | Illumina_V1 | 
--  |           3256 | Illumina_V2 | 
--  |           4300 | AgilentProbe | 
--  |           4305 | AgilentCGH   | 
--  |           6000 | Codelink | 


DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id=3250 and edb.db_name='Illumina' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;
DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id=3255 and edb.db_name='Illumina_V1' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;
DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id=3256 and edb.db_name='Illumina_V2' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;
DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id=4300 and edb.db_name='AgilentProbe' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;
DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id=4306 and edb.db_name='AgilentCGH' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;
DELETE ox from external_db edb, xref x, object_xref ox where edb.external_db_id=6000 and edb.db_name='Codelink' and edb.external_db_id=x.external_db_id and x.xref_id=ox.xref_id;


DELETE x from external_db edb, xref x where edb.external_db_id=3250 and edb.db_name='Illumina' and edb.external_db_id=x.external_db_id ;
DELETE x from external_db edb, xref x where edb.external_db_id=3255 and edb.db_name='Illumina_V1' and edb.external_db_id=x.external_db_id ;
DELETE x from external_db edb, xref x where edb.external_db_id=3256 and edb.db_name='Illumina_V2' and edb.external_db_id=x.external_db_id ;
DELETE x from external_db edb, xref x where edb.external_db_id=4300 and edb.db_name='AgilentProbe' and edb.external_db_id=x.external_db_id ;
DELETE x from external_db edb, xref x where edb.external_db_id=4305 and edb.db_name='AgilentCGH' and edb.external_db_id=x.external_db_id ;
DELETE x from external_db edb, xref x where edb.external_db_id=6000 and edb.db_name='Codelink' and edb.external_db_id=x.external_db_id ;




DELETE from external_db  where external_db_id=3250 and db_name='Illumina';
DELETE from external_db  where external_db_id=3255 and db_name='Illumina_V1';
DELETE from external_db  where external_db_id=3256 and db_name='Illumina_V2';
DELETE from external_db  where external_db_id=4300 and db_name='AgilentProbe';
DELETE from external_db  where external_db_id=4305 and db_name='AgilentCGH';
DELETE from external_db  where external_db_id=6000 and db_name='Codelink';


-- Delete unmapped_reasons

DELETE ur from unmapped_reason ur, unmapped_object uo where uo.type='probe2transcript' and uo.unmapped_reason_id=ur.unmapped_reason_id;

-- Delete unmapped_objects, separately just in case of orhpans
DELETE from unmapped_object where type='probe2transcript';


-- Drop oligo tables last so we do not get errors, should we need to rerun

DROP table oligo_array;
DROP table oligo_probe;
DROP table oligo_feature;

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_55_56_c.sql|drop_oligo_tables_and_xrefs');

