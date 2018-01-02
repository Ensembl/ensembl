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

# patch_72_73_b.sql
#
# Title: Introduce alt allele types
#
# Description: Copy out existing alt allele data, create new tables and
#              integrate old data into the new tables

# Relocate existing data out of the way
CREATE TABLE aa_bak LIKE alt_allele;
INSERT INTO aa_bak SELECT * FROM alt_allele;

DROP TABLE alt_allele;

# Make new table structure
CREATE TABLE alt_allele (alt_allele_id INT UNSIGNED AUTO_INCREMENT, 
                         alt_allele_group_id INT UNSIGNED NOT NULL, 
                         gene_id INT UNSIGNED NOT NULL,
                         PRIMARY KEY (alt_allele_id),
                         KEY (gene_id,alt_allele_group_id)
                         ) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

CREATE TABLE alt_allele_attrib (alt_allele_id INT UNSIGNED,
                                attrib ENUM('IS_REPRESENTATIVE',
                                            'IS_MOST_COMMON_ALLELE',
                                            'IN_CORRECTED_ASSEMBLY',
                                            'HAS_CODING_POTENTIAL',
                                            'IN_ARTIFICIALLY_DUPLICATED_ASSEMBLY',
                                            'IN_SYNTENIC_REGION',
                                            'HAS_SAME_UNDERLYING_DNA_SEQUENCE',
                                            'IN_BROKEN_ASSEMBLY_REGION',
                                            'IS_VALID_ALTERNATE',
                                            'SAME_AS_REPRESENTATIVE',
                                            'SAME_AS_ANOTHER_ALLELE',
                                            'MANUALLY_ASSIGNED',
                                            'AUTOMATICALLY_ASSIGNED'),
                                KEY aa_idx (alt_allele_id,attrib)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;;

CREATE TABLE alt_allele_group (alt_allele_group_id INT UNSIGNED AUTO_INCREMENT,
                               PRIMARY KEY (alt_allele_group_id)
                               ) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

# Port data into new structure
INSERT INTO alt_allele_group (alt_allele_group_id) SELECT DISTINCT alt_allele_id FROM aa_bak;
INSERT INTO alt_allele (alt_allele_group_id,gene_id) SELECT alt_allele_id,gene_id FROM aa_bak;

INSERT INTO alt_allele_attrib (alt_allele_id,attrib) SELECT a.alt_allele_id,'IS_REPRESENTATIVE' FROM alt_allele a, aa_bak b 
    WHERE b.is_ref = 1 AND a.gene_id = b.gene_id AND a.alt_allele_group_id = b.alt_allele_id;

# Clean up remains
DROP TABLE aa_bak;


# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_72_73_b.sql|alt_allele_type');
