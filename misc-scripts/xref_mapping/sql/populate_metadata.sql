# Populate the appropriate tables in an xref metadata database


################################################################################
# SPECIES

INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (9606,9606,  'homo_sapiens',            'human,hsapiens,homosapiens');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (10090,10090, 'mus_musculus',            'mouse,mmusculus,musmusculus');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (10116, 10116, 'rattus_norvegicus',       'rat,rnovegicus,rattusnorvegicus');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (31033,31033, 'fugu_rubripes',           'pufferfish,fugu,frubripes,fugurubripes');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (7227, 7227, 'drosophila_melanogaster', 'drosophila,dmelongaster,drosophilamelanogaster' );
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (6239, 6239, 'caenorhabditis_elegans',  'elegans,celegans,caenorhabditiselegans');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (6238, 6238, 'caenorhabditis_briggsae', 'briggsae,cbriggsae,caenorhabditisbriggsae');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (7955, 7955, 'danio_rerio',             'zebrafish,danio,drerio,daniorerio' );
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (9598, 9598, 'pan_troglodytes',         'chimp,chimpanzee,ptroglodytes,pantroglodytes');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (9031, 9031,  'gallus_gallus',           'chicken,chick,ggallus,gallusgallus' );
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (99883, 99883,'tetraodon_nigroviridis',  'tetraodon,tnigroviridis,tetraodonnigroviridis');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (9913, 9913,  'bos_taurus',             'cow,btaurus,bostaurus');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (9615, 9615,  'canis_familiaris',        'dog,doggy,cfamiliaris,canisfamiliaris');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (8364, 8364,  'xenopus_tropicalis',        'pipid,pipidfrog,xenopus,xtropicalis,xenopustropicalis');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (13616, 13616,  'monodelphis_domestica',        'opossum,monodelphis,mdomestica,monodelphisdomestica');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (4932, 4932,  'saccharomyces_cerevisiae',  'yeast,saccharomyces,scerevisiae,saccharomycescerevisiae');

INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (7460, 7460,  'apis_mellifera',  'apismellifera,honeybee,amellifera');

INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (7719, 7719,  'ciona_intestinalis', 'cionaintestinalis,cintestinalis,seasquirt'); 

# Serveral strains of anopheles have different taxonomy IDs
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (7165,7165, 'anopheles_gambiae', 'mosquito,anopheles,agambiae,anophelesgambiae');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (7165,180454, 'anopheles_gambiae_strain', '');


################################################################################
# SOURCES - types of data we can read

# "High level" sources that we will also download from (via source_url)

INSERT INTO source VALUES (1, "Uniprot/SWISSPROT", 1, 'Y',1);
INSERT INTO source VALUES (2, "Uniprot/SPTREMBL", 1, 'Y',1);
INSERT INTO source VALUES (3, "RefSeq_peptide", 1, 'Y',1);
INSERT INTO source VALUES (4, "RefSeq_dna", 1, 'Y',1);
INSERT INTO source VALUES (5, "IPI", 1, 'Y',2);
INSERT INTO source VALUES (6, "UniGene", 1, 'Y',2);
INSERT INTO source VALUES (10, "RefSeq_peptide_predicted", 1, 'Y',1);
INSERT INTO source VALUES (11, "RefSeq_dna_predicted", 1, 'Y',1);

# Agilent probes
INSERT INTO source VALUES (2700, 'AgilentProbe', 1, 'Y', 4);
INSERT INTO source VALUES (2701, 'AgilentCGH', 1, 'Y', 4);

# Other sources - used to create dependent xrefs, but not to download from

INSERT INTO source VALUES (1010, 'EMBL', 1, 'N', 2);
INSERT INTO source VALUES (1020, 'MIM', 1, 'N', 2);
INSERT INTO source VALUES (1030, 'PDB', 1, 'N', 2);
INSERT INTO source VALUES (1040, 'protein_id', 1, 'N', 2);
INSERT INTO source VALUES (1050, 'PUBMED', 1, 'N', 2);
INSERT INTO source VALUES (1060, 'MEDLINE', 1, 'N', 2);
INSERT INTO source VALUES (1070, 'GO', 1, 'Y',5);
INSERT INTO source VALUES (1080, 'MarkerSymbol', 1, 'Y',2);
INSERT INTO source VALUES (1090, 'HUGO', 1, 'Y',2);

#INSERT INTO source VALUES (1100, 'LocusLink', 1, 'N', 2);
INSERT INTO source VALUES (1110, 'EntrezGene', 1, 'N', 2);

INSERT INTO source VALUES (1200, 'RGD', 1, 'Y',2);
INSERT INTO source VALUES (1300, 'Interpro', 1, 'Y', 2);
INSERT INTO source VALUES (1400, 'ZFIN_ID', 1, 'Y', 2);
#INSERT INTO source VALUES (1500, 'OMIM', 1, 'Y', 3);

INSERT INTO source VALUES (2000, 'CCDS', 1, 'Y', 4);

#INSERT INTO source VALUES (2400, 'WormBase', 1, 'Y',4);
INSERT INTO source VALUES (2400, 'wormpep_id', 1, 'Y', 4);
INSERT INTO source VALUES (2410, 'wormbase_gene', 1, 'N',4);
INSERT INTO source VALUES (2420, 'wormbase_transcript', 1, 'N', 4);
INSERT INTO source VALUES (2440, 'wormbase_pseudogene', 1, 'N', 4);

# drosphila melanogster sources 
INSERT INTO source VALUES (2500, 'flybase_gff', 1, 'Y', 4);
INSERT INTO source VALUES (2510, 'flybase_gene_id', 1, 'N', 4);
INSERT INTO source VALUES (2520, 'flybase_transcript_id', 1, 'N', 4);
INSERT INTO source VALUES (2530, 'flybase_polypeptide_id', 1, 'N', 4);
INSERT INTO source VALUES (2540, 'flybase_annotation_id', 1, 'N', 4);
INSERT INTO source VALUES (2550, 'flybase_synonym', 1, 'N', 4);
INSERT INTO source VALUES (2560, 'flybase_name', 1, 'N', 4);
INSERT INTO source VALUES (2570, 'gadfly_gene_cgid', 1, 'N', 4);
INSERT INTO source VALUES (2571, 'gadfly_transcript_cgid', 1, 'N', 4);
INSERT INTO source VALUES (2572, 'gadfly_translation_cgid', 1, 'N', 4);

# ciona intestinalis source 
INSERT INTO source VALUES (2601, 'c_int_proteins_jgi_v1', 1, 'Y', 4);

# anopheles gambiae sources
# predicted versions of Uniprot/EMBL - created as dependent in appropriate parse
INSERT INTO source VALUES (2801, "Uniprot/SWISSPROT_predicted", 1, 'N',1);
INSERT INTO source VALUES (2802, "Uniprot/SPTREMBL_predicted", 1, 'N',1);
INSERT INTO source VALUES (2810, "EMBL_predicted", 1, 'N',1);
INSERT INTO source VALUES (2840, 'protein_id_predicted', 1, 'N', 2);
INSERT INTO source VALUES (2850, "Celera_Gene", 1, 'Y', 3);
INSERT INTO source VALUES (2852, "Celera_Pep", 1, 'Y', 3);
INSERT INTO source VALUES (2854, "Celera_Trans", 1, 'Y', 3);
INSERT INTO source VALUES (2856, "Anopheles_symbol", 1, 'Y', 3);

################################################################################
# Files to fetch data from


# --------------------------------------------------------------------------------
### DROSOPHILA TESTING GFF
# chr 2L	
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gff/dmel-2L-r4.1.gff.gz', 'N', now(), now(), "Flybase_dmel_GFFv3_Parser");

# chr 2R
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gff/dmel-2R-r4.1.gff.gz', 'N', now(), now(), "Flybase_dmel_GFFv3_Parser");

# chr 3L
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gff/dmel-3L-r4.1.gff.gz', 'N', now(), now(), "Flybase_dmel_GFFv3_Parser");

# chr 3R
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gff/dmel-3R-r4.1.gff.gz', 'N', now(), now(), "Flybase_dmel_GFFv3_Parser");

# chr 4
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gff/dmel-4-r4.1.gff.gz', 'N', now(), now(), "Flybase_dmel_GFFv3_Parser");

# chr X 
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gff/dmel-X-r4.1.gff.gz', 'N', now(), now(), "Flybase_dmel_GFFv3_Parser");

## HETERO-CHROMATIN

# 2h 
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current_hetchr/gff/dmel-2h-hetr32b2.gff.gz', 'N', now(), now(), "Flybase_dmel_GFFv3_Parser");

# 3h 
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current_hetchr/gff/dmel-3h-hetr32b2.gff.gz', 'N', now(), now(), "Flybase_dmel_GFFv3_Parser");

# 4h 
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current_hetchr/gff/dmel-4h-hetr32b2.gff.gz', 'N', now(), now(), "Flybase_dmel_GFFv3_Parser");

# U 
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current_hetchr/gff/dmel-U-hetr32b2.gff.gz', 'N', now(), now(), "Flybase_dmel_GFFv3_Parser");

# Xh 
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current_hetchr/gff/dmel-Xh-hetr32b2.gff.gz', 'N', now(), now(), "Flybase_dmel_GFFv3_Parser");

# Yh	 
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current_hetchr/gff/dmel-Yh-hetr32b2.gff.gz', 'N', now(), now(), "Flybase_dmel_GFFv3_Parser");


# Uniprot for drosophila
 
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1,7227,'ftp://ftp.ebi.ac.uk/pub/databases/integr8/uniprot/proteomes/17.D_melanogaster.dat.gz',now(),now(),'UniProtParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1,7227,'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz',now(),now(),'UniProtParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (2,7227,'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz',now(),now(),'UniProtParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,7227,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate1.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,7227,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate2.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,7227,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate3.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,7227,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate4.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,7227,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate5.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,7227,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate6.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,7227,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate7.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,7227,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate8.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,7227,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate9.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1070,7227,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz',now(),now(),'GOParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1300,7227,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz',now(),now(),'InterproParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (6,7227,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Drosophila_melanogaster/Dm.seq.uniq.gz ftp://ftp.ncbi.nih.gov/repository/UniGene/Drosophila_melanogaster/Dm.data.gz',now(),now(),'UniGeneParser');
















### HUMAN
##       uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 9606, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', '', now(), now(), "UniProtParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 9606, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', '', now(), now(), "UniProtParser");


##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 9606,'ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4, 9606,'ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.rna.fna.gz', '', now(), now(), "RefSeqParser");

##       GO
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1070, 9606,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/gene_association.goa_human.gz', '', now(), now(), "GOParser");

##       HUGO
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1090, 9606,'http://www.gene.ucl.ac.uk/cgi-bin/nomenclature/gdlw.pl?title=Genew+output+data&col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_prev_sym&col=md_prot_id&status=Approved&status=Approved+Non-Human&status_opt=3&=on&where=&order_by=gd_hgnc_id&limit=&format=text&submit=submit&.cgifields=&.cgifields=status&.cgifields=chr', '', now(), now(), "HUGOParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1300, 9606,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', '', now(), now(), "InterproParser");

##      OMIM 
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1500, 9606,'ftp://ftp.ncbi.nih.gov/repository/OMIM/morbidmap', '', now(), now(), "MIMParser");

##      IPI
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (5, 9606,'ftp://ftp.ebi.ac.uk/pub/databases/IPI/current/ipi.HUMAN.fasta.gz', '', now(), now(), "IPIParser");

##      CCDS
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2000, 9606,'/dummy/CCDS.txt', '', now(), now(), "CCDSParser");


##      Agilent
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2700, 9606,'./Agilent/HumanExpression.fasta', '', now(), now(), "AgilentParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2701, 9606,'./Agilent/HumanCGH.fasta', '', now(), now(), "AgilentParser");

##      UniGene
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (6, 9606,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Homo_sapiens/Hs.seq.uniq.gz ftp://ftp.ncbi.nih.gov/repository/UniGene/Homo_sapiens/Hs.data.gz', '', now(), now(), "UniGeneParser");

# --------------------------------------------------------------------------------

###MOUSE
##      uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 10090, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', '', now(), now(), "UniProtParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 10090, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', '', now(), now(), "UniProtParser");


##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 10090,'ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4, 10090,'ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.rna.fna.gz', '', now(), now(), "RefSeqParser");

##      mgd (MGI -- MarkerSymbol)
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1080, 10090,'ftp://ftp.informatics.jax.org/pub/reports/MRK_SwissProt_TrEMBL.rpt ftp://ftp.informatics.jax.org/pub/reports/MRK_Synonym.sql.rpt', '', now(), now(), "MGDParser");

##      GO 
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1070, 10090,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/MOUSE/gene_association.goa_mouse.gz', '', now(), now(), "GOParser");

##      IPI
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (5, 10090,'ftp://ftp.ebi.ac.uk/pub/databases/IPI/current/ipi.MOUSE.fasta.gz', '', now(), now(), "IPIParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1300, 10090,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', '', now(), now(), "InterproParser");

##      UniGene
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (6, 10090,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Mus_musculus/Mm.seq.uniq.gz ftp://ftp.ncbi.nih.gov/repository/UniGene/Mus_musculus/Mm.data.gz', '', now(), now(), "UniGeneParser");

##      Agilent
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2700, 10090,'./Agilent/MouseExpression.fasta', '', now(), now(), "AgilentParser");

# --------------------------------------------------------------------------------
### RAT
##      uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 10116, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', '', now(), now(), "UniProtParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 10116, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', '', now(), now(), "UniProtParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 10116,'ftp://ftp.ncbi.nih.gov/refseq/R_norvegicus/mRNA_Prot/rat.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4, 10116,'ftp://ftp.ncbi.nih.gov/refseq/R_norvegicus/mRNA_Prot/rat.rna.fna.gz', '', now(), now(), "RefSeqParser");


##      GO 
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1070, 10116,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/RAT/gene_association.goa_rat.gz', '', now(), now(), "GOParser");

##  RGD
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1200, 10116,'ftp://rgd.mcw.edu/pub/data_release/genbank_to_gene_ids.txt', '', now(), now(), "RGDParser");

##  IPI
##      IPI
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (5, 10116,'ftp://ftp.ebi.ac.uk/pub/databases/IPI/current/ipi.RAT.fasta.gz', '', now(), now(), "IPIParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1300, 10116,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', '', now(), now(), "InterproParser");

##      UniGene
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (6, 10116,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Rattus_norvegicus/Rn.seq.uniq.gz ftp://ftp.ncbi.nih.gov/repository/UniGene/Rattus_norvegicus/Rn.data.gz', '', now(), now(), "UniGeneParser");

##      Agilent
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2700, 10116,'./Agilent/RatExpression.fasta', '', now(), now(), "AgilentParser");

# --------------------------------------------------------------------------------

### Zebrafish
##      uniprot

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 7955, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', '', now(), now(), "UniProtParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 7955, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', '', now(), now(), "UniProtParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 7955,'ftp://ftp.ncbi.nih.gov/refseq/D_rerio/mRNA_Prot/zebrafish.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4, 7955,'ftp://ftp.ncbi.nih.gov/refseq/D_rerio/mRNA_Prot/zebrafish.rna.fna.gz', '', now(), now(), "RefSeqParser");


##      GO 
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1070, 7955,'ftp://ftp.geneontology.org/pub/go/gene-associations/gene_association.zfin.gz', '', now(), now(), "GOParser");

##      ZFIN
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1400, 7955,'http://zfin.org/data_transfer/Downloads/refseq.txt http://zfin.org/data_transfer/Downloads/swissprot.txt', '', now(), now(), "ZFINParser");

##      IPI
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (5, 7955,'ftp://ftp.ebi.ac.uk/pub/databases/IPI/current/ipi.BRARE.fasta.gz', '', now(), now(), "IPIParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1300, 7955,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', '', now(), now(), "InterproParser");

##      UniGene
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (6, 7955,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Danio_rerio/Dr.seq.uniq.gz ftp://ftp.ncbi.nih.gov/repository/UniGene/Danio_rerio/Dr.data.gz', '', now(), now(), "UniGeneParser");

# --------------------------------------------------------------------------------
### Chicken
##      uniprot

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 9031, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', '', now(), now(), "UniProtParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 9031, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', '', now(), now(), "UniProtParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 9031,'ftp://ftp.ncbi.nih.gov/genomes/Gallus_gallus/protein/protein.gbk.gz', '', now(), now(), "RefSeqGPFFParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4, 9031,'ftp://ftp.ncbi.nih.gov/genomes/Gallus_gallus/RNA/rna.gbk.gz', '', now(), now(), "RefSeqGPFFParser");

##       GO
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1070, 9031,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz', '', now(), now(), "GOParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1300, 9031,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', '', now(), now(), "InterproParser");

##      UniGene
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (6, 9031,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Gallus_gallus/Gga.seq.uniq.gz ftp://ftp.ncbi.nih.gov/repository/UniGene/Gallus_gallus/Gga.data.gz', '', now(), now(), "UniGeneParser");

# --------------------------------------------------------------------------------
### XENOPUS
# uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 8364, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', '', now(), now(), "UniProtParser");

#uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 8364, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', '', now(), now(), "UniProtParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 8364,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other1.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 8364,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other2.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 8364,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other3.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4, 8364,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other1.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4, 8364,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other2.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4, 8364,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other3.rna.fna.gz', '', now(), now(), "RefSeqParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1300, 8364,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', '', now(), now(), "InterproParser");

##      UniGene
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (6, 8364,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Xenopus_tropicalis/Str.seq.uniq.gz ftp://ftp.ncbi.nih.gov/repository/UniGene/Xenopus_tropicalis/Str.data.gz', '', now(), now(), "UniGeneParser");

# --------------------------------------------------------------------------------

### Dog

#        uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 9615, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', '', now(), now(), "UniProtParser");

#        uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 9615, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', '', now(), now(), "UniProtParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 9615,'ftp://ftp.ncbi.nih.gov/genomes/Canis_familiaris/protein/protein.gbk.gz', '', now(), now(), "RefSeqGPFFParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4, 9615,'ftp://ftp.ncbi.nih.gov/genomes/Canis_familiaris/RNA/rna.gbk.gz', '', now(), now(), "RefSeqGPFFParser");

##       GO
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1070, 9615,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz', '', now(), now(), "GOParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1300, 9615,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', '', now(), now(), "InterproParser");

##      UniGene
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (6, 9615,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Canis_familiaris/Cfa.seq.uniq.gz ftp://ftp.ncbi.nih.gov/repository/UniGene/Canis_familiaris/Cfa.data.gz', '', now(), now(), "UniGeneParser");

# --------------------------------------------------------------------------------

### C elegans

#        uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 6239, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', '', now(), now(), "UniProtParser");

#        uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 6239, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', '', now(), now(), "UniProtParser");

##       GO
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1070, 6239,'ftp://ftp.geneontology.org/pub/go/gene-associations/gene_association.wb.gz', '', now(), now(), "GOParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1300, 6239,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', '', now(), now(), "InterproParser");

##      UniGene
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (6, 6239,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Caenorhabditis_elegans/Cel.seq.uniq.gz ftp://ftp.ncbi.nih.gov/repository/UniGene/Caenorhabditis_elegans/Cel.data.gz', '', now(), now(), "UniGeneParser");

##       refseq
#INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4, 6239,'ftp://ftp.ncbi.nih.gov/genomes/Caenorhabditis_elegans/CHR_I/NC_003279.gbk', '', now(), now(), "RefSeqGPFFParser");
#INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4, 6239,'ftp://ftp.ncbi.nih.gov/genomes/Caenorhabditis_elegans/CHR_II/NC_003280.gbk', '', now(), now(), "RefSeqGPFFParser");
#INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4, 6239,'ftp://ftp.ncbi.nih.gov/genomes/Caenorhabditis_elegans/CHR_III/NC_003281.gbk', '', now(), now(), "RefSeqGPFFParser");
#INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4, 6239,'ftp://ftp.ncbi.nih.gov/genomes/Caenorhabditis_elegans/CHR_IV/NC_003282.gbk', '', now(), now(), "RefSeqGPFFParser");
#INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4, 6239,'ftp://ftp.ncbi.nih.gov/genomes/Caenorhabditis_elegans/CHR_V/NC_003283.gbk', '', now(), now(), "RefSeqGPFFParser");
#INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4, 6239,'ftp://ftp.ncbi.nih.gov/genomes/Caenorhabditis_elegans/CHR_X/NC_003284.gbk', '', now(), now(), "RefSeqGPFFParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 6239,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate1.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 6239,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate2.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 6239,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate3.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 6239,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate4.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 6239,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate5.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 6239,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate6.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 6239,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate7.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 6239,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate8.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 6239,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate9.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");

##   WormBase 

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2400, 6239, 'ftp://ftp.sanger.ac.uk/pub/databases/wormpep/wormpep140/wormpep.table140', '', now(), now(), "WormPepParser");

##   Stable ID xref transfer - note use of WormbaseDatabaseStableIDParser
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2400, 6239, 'mysql:ecs4:3350:glenn_elegans_140:ensro', '', now(), now(), "WormbaseDatabaseStableIDParser");

# --------------------------------------------------------------------------------

#### HoneyBee

## Uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 7460, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', '', now(), now(), "UniProtParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 7460, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', '', now(), now(), "UniProtParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 7460,'ftp://ftp.ncbi.nih.gov/genomes/Apis_mellifera/protein/protein.gbk.gz', '', now(), now(), "RefSeqGPFFParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4, 7460,'ftp://ftp.ncbi.nih.gov/genomes/Apis_mellifera/RNA/rna.gbk.gz', '', now(), now(), "RefSeqGPFFParser");

##       GO
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1070, 7460,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz', '', now(), now(), "GOParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1300, 7460,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', '', now(), now(), "InterproParser");

##      UniGene
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (6, 7460,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Apis_mellifera/Ame.seq.uniq.gz ftp://ftp.ncbi.nih.gov/repository/UniGene/Apis_mellifera/Ame.data.gz', '', now(), now(), "UniGeneParser");

# --------------------------------------------------------------------------------


#### SeaSquirt

## Uniprot
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (1, 7719, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (2, 7719, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', now(), now(), "UniProtParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 7719,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other1.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 7719,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other2.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 7719,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other3.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 7719,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other1.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 7719,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other2.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 7719,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other3.rna.fna.gz',  now(), now(), "RefSeqParser");

##       GO
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(1070, 7719,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz', now(), now(), "GOParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(1300, 7719,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', now(), now(), "InterproParser");


##      UniGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(6, 7719,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Ciona_intestinalis/Cin.seq.uniq.gz ftp://ftp.ncbi.nih.gov/repository/UniGene/Ciona_intestinalis/Cin.data.gz',  now(), now(), "UniGeneParser");

## xrefs to protein-annotation from JGI on OLDER assembly version v1 
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (2601, 7719,'ftp://ftp.jgi-psf.org/pub/JGI_data/Ciona/v1.0/ciona.prot.fasta.gz',now(),now(),'FastaProteinSeqParser');


# --------------------------------------------------------------------------------

###### Anopheles

# Uniprot
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1,7165,'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz',now(),now(),'UniProtParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (2,7165,'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz',now(),now(),'UniProtParser');

# Interpro
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1300,7165,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz',now(),now(),'InterproParser');

# Unigene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (6,7165,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Anopheles_gambiae/Aga.seq.uniq.gz  ftp://ftp.ncbi.nih.gov/repository/UniGene/Anopheles_gambiae/Aga.data.gz',now(),now(),'UniGeneParser');

# Celera
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (2852, 7165,'/Celera_Pep/consensus-proteins_xref.fasta', now(), now(), "CeleraProteinParser");

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (2854, 7165,'/Celera_trans/consensus-transcripts_xref.fasta', now(), now(), "CeleraTranscriptParser");

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (2856, 7165,'/Anopheles_symbol/GeneName_translation_UniqID.txt', now(), now(), "AnophelesSymbolParser");

# --------------------------------------------------------------------------------



#### Cow

## Uniprot
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (1, 9913, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (2, 9913, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', now(), now(), "UniProtParser");

##       refseq protein
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 9913,'ftp://ftp.ncbi.nih.gov/genomes/Bos_taurus/protein/protein.gbk.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4, 9913,'ftp://ftp.ncbi.nih.gov/genomes/Bos_taurus/RNA/rna.gbk.gz', '', now(), now(), "RefSeqGPFFParser");

##       GO
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(1070, 9913,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz', now(), now(), "GOParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(1300, 9913,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', now(), now(), "InterproParser");

##      UniGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(6, 9913,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Bt.seq.uniq.gz ftp://ftp.ncbi.nih.gov/repository/UniGene/Bt.data.gz',  now(), now(), "UniGeneParser");


################################################################################

