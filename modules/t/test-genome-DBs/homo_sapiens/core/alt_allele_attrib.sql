DROP TABLE IF EXISTS `alt_allele_attrib`;

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
);