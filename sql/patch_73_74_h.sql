# patch_73_74_h.sql
#
# Title: Creating a unique index on alt_allele(gene_id)
#
# Description:
#
# Alt alleles only allow a gene to appear in one group. Enforce this 
# at the DB level

CREATE UNIQUE INDEX gene_idx ON alt_allele(gene_id);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_73_74_h.sql|alt_allele_unique_gene_idx');

 
