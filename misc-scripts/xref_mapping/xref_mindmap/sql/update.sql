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


LOCK TABLES `external_db_type` WRITE;
UPDATE external_db_type SET db_type_id = 1 WHERE db_name IN ('MIM','MIM_GENE','MIM_MORBID','DBASS3','DBASS5','Orphanet');
UPDATE external_db_type SET db_type_id = 2 WHERE db_name IN ('UniGene', 'HPA');
UPDATE external_db_type SET db_type_id = 3 WHERE db_name IN ('EntrezGene','Uniprot/SPTREMBL','Uniprot/SPTREMBL_predicted','Uniprot/SWISSPROT','Uniprot/SWISSPROT_predicted','Uniprot/Varsplic', 'UniProtKB_all','MEROPS','WikiGene');
UPDATE external_db_type SET db_type_id = 4 WHERE db_name IN ('HGNC','MGI','ZFIN_ID','HGNC_curated_gene','HGNC_automatic_gene','Clone_based_vega_gene','Clone_based_ensembl_gene','HGNC_curated_transcript','HGNC_automatic_transcript','Clone_based_vega_transcript','Clone_based_ensembl_transcript','MGI_curated_gene','MGI_automatic_gene','MGI_curated_transcript','MGI_automatic_transcript','RFAM_gene_name','miRBase_gene_name','miRBase_transcript_name','RFAM_transcript_name','HGNC_transcript_name','MGI_transcript_name','ZFIN_ID_transcript_name','Uniprot_genename','RefSeq_gene_name');
UPDATE external_db_type SET db_type_id = 4 WHERE db_name like '%gene\_name' or db_name like '%transcript\_name' or db_name like '%translation\_name' or db_name like '%name\_gene' or  db_name like '%name\_transcript' or db_name like '%name\_translation';  
UPDATE external_db_type SET db_type_id = 5 WHERE db_name IN ('miRBase','miRBase_predicted','RFAM');
UPDATE external_db_type SET db_type_id = 6 WHERE db_name IN ('EMBL','EMBL_predicted','protein_id','protein_id_predicted','RefSeq_dna','RefSeq_mRNA','RefSeq_ncRNA','RefSeq_dna_predicted','RefSeq_mRNA_predicted','RefSeq_ncRNA_predicted','RefSeq_peptide','RefSeq_peptide_predicted','RefSeq_rna','RefSeq_rna_predicted','RefSeq_genomic','Vega_gene','Vega_gene_like','Vega_transcript','Vega_transcript_like','Vega_translation','Ens_Hs_gene','Ens_Hs_transcript','Ens_Hs_translation','IPI', 'OTTG','CCDS','OTTT','OTTP','shares_CDS_with','shares_CDS_with_ENST','shares_CDS_with_OTTT','shares_CDS_and_UTR_with_OTTT','ENSG','ENST','ENST_ident','ENST_CDS','UCSC','EMBLBANK_GENE', 'EMBLBANK_TRANSCRIPT','LRG','ENS_LRG_gene','ENS_LRG_transcript');
UPDATE external_db_type SET db_type_id = 7 WHERE db_name IN ('Interpro','PFAM');
UPDATE external_db_type SET db_type_id = 8 WHERE db_name IN ('GO','goslim_goa','goslim_generic');
UPDATE external_db_type SET db_type_id = 9 WHERE db_name IN ('PDB');
UPDATE external_db_type SET db_type_id = 10 WHERE db_type_id IS NULL;
UNLOCK TABLES;
