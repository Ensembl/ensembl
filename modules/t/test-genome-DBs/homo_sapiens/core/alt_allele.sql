--
-- Table structure for table `alt_allele`
--

DROP TABLE IF EXISTS `alt_allele`;

CREATE TABLE alt_allele (alt_allele_id INT UNSIGNED AUTO_INCREMENT, 
                         alt_allele_group_id INT UNSIGNED NOT NULL, 
                         gene_id INT UNSIGNED NOT NULL,
                         PRIMARY KEY (alt_allele_id),
                         KEY (gene_id,alt_allele_group_id)
                         );