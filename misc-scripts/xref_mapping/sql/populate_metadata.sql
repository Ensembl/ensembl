# Populate the appropriate tables in an xref metadata database


################################################################################
# SPECIES

INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (9606,9606,  'homo_sapiens',            'human,hsapiens,homosapiens');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (10090,10090, 'mus_musculus',            'mouse,mmusculus,musmusculus');
INSERT INTO species (species_id,taxonomy_id, name, aliases) VALUES (10116, 10116, 'rattus_norvegicus',       'rat,rnovegicus,rattusnorvegicus');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (31033,31033, 'takifugu_rubripes',           'pufferfish,fugu,fugu_rubripes,frubripes,fugurubripes');
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
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (51511, 51511,'ciona_savignyi',     'ciona_savi,csav,savi'); 

# Serveral strains of anopheles have different taxonomy IDs
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (7165,7165, 'anopheles_gambiae', 'mosquito,anopheles,agambiae,anophelesgambiae');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (7165,180454, 'anopheles_gambiae_strain', '');

INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (13616,13616, 'monodelphis_domestica', 'mono,opossum,possum,monodelphis,monodelphisdomestica');

INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (9544, 9544,  'macaca_mulatta',        'rmacaque, rhesus, rhesus macaque, macaque');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (9785, 9785,  'loxodonta_africana',    'elephant,lafricana,loxodontaafricana');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (9361, 9361,  'dasypus_novemcinctus',  'armadillo,dnovemcinctus,dasypusnovemcinctus');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (9986, 9986,  'oryctolagus_cuniculus', 'rabbit,ocuniculus,oryctolaguscuniculus');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (9371, 9371,  'echinops_telfairi',     'tenrec,etelfairi,echinopstelfairi');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (9685, 9685,  'felis_catus',           'cat,fcatus,feliscatus,moggy,tiddles');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (7159, 7159,  'aedes_aegypti',        'aedes,aedesaegypti,aaegypti');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (69293,69293,  'gasterosteus_aculeatus',  'stickleback, gaculeatus, gasterosteusaculeatus');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (8090, 8090, 'oryzias_latipes', 'medaka, japanese ricefish, japanese medaka');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (59463, 59463, 'myotis_lucifugus','microbat, little brown bat');
INSERT INTO species (species_id, taxonomy_id, name, aliases) VALUES (9258,9258, 'ornithorhynchus_anatinus', 'platypus,duckbilledplatypus,ornithorhynchusanatinus duckbilled_platypus');
################################################################################
# SOURCES - types of data we can read

# "High level" sources that we will also download from (via source_url)

INSERT INTO source VALUES (1020, 'MIM', 1, 'Y', 10, 1, "");
INSERT INTO source VALUES (2000, 'CCDS', 1, 'Y', 10, 1, "");
INSERT INTO source VALUES (1110, 'EntrezGene', 1, 'Y', 10, 1, "");
INSERT INTO source VALUES (1250, 'SGD', 1, 'Y',10, 1, "");



INSERT INTO source VALUES (1,  "Uniprot/SWISSPROT", 1, 'Y',20,1, "" );
INSERT INTO source VALUES (2,  "Uniprot/SPTREMBL", 1, 'Y',20, 1, "");
INSERT INTO source VALUES (3,  "RefSeq_peptide", 1, 'Y',20, 1, "");

INSERT INTO source VALUES (4,  "RefSeq_dna", 1, 'Y',20, 2, "refseq");
INSERT INTO source VALUES (14, "RefSeq_dna", 1, 'Y',20, 1, "ccds");

INSERT INTO source VALUES (5,  "IPI", 1, 'Y',30, 1, "");
INSERT INTO source VALUES (6,  "UniGene", 1, 'Y',30, 1, "");
INSERT INTO source VALUES (10, "RefSeq_peptide_predicted", 1, 'Y',20, 1, "");
INSERT INTO source VALUES (11, "RefSeq_dna_predicted", 1, 'Y',20, 1, "");
# Agilent probes
INSERT INTO source VALUES (2700, 'AgilentProbe', 1, 'Y', 50, 1, "");
INSERT INTO source VALUES (2701, 'AgilentCGH', 1, 'Y', 50, 1, "");

# Other sources - used to create dependent xrefs, but not to download from

INSERT INTO source VALUES (1010, 'EMBL', 1, 'N', 30, 1, "");
INSERT INTO source VALUES (1021, 'MIM_GENE', 1, 'N', 30, 1, "");
INSERT INTO source VALUES (1022, 'MIM_MORBID', 1, 'N', 30, 1, "");
INSERT INTO source VALUES (1030, 'PDB', 1, 'N', 30, 1, "");
INSERT INTO source VALUES (1040, 'protein_id', 1, 'N', 30, 1, "");
INSERT INTO source VALUES (1050, 'PUBMED', 1, 'N', 30, 1, "");
INSERT INTO source VALUES (1060, 'MEDLINE', 1, 'N', 30, 1, "");
INSERT INTO source VALUES (1070, 'GO', 1, 'Y',60, 1, "");
INSERT INTO source VALUES (1080, 'MarkerSymbol', 1, 'Y',30, 1, "");

# note for hugo 1090 has a lower order ie. it will be parsed after the 1091
# but has a higher priority as it should take precedense.
# 1090 needs 1091 to be loaded first to obtain descriptions etc.

INSERT INTO source VALUES (1090, 'HUGO', 1, 'Y',30, 1, "havana"); # use havana
INSERT INTO source VALUES (1093, 'HUGO', 1, 'Y',51, 2, "ccds"); # use CCDS
INSERT INTO source VALUES (1091, 'HUGO', 1, 'Y',29, 3, "refseq"); # use refseq
INSERT INTO source VALUES (1092, 'HUGO', 1, 'N',29, 4, "uniprot"); # use uniprot
INSERT INTO source VALUES (1094, 'HUGO', 1, 'N',29, 5, "entrezgene"); # use entrezgene



INSERT INTO source VALUES (1200, 'RGD', 1, 'Y',30, 1, "");
INSERT INTO source VALUES (1300, 'Interpro', 1, 'Y', 30, 1, "");
INSERT INTO source VALUES (1400, 'ZFIN_ID', 1, 'Y', 30, 1, "");

#INSERT INTO source VALUES (2400, 'WormBase', 1, 'Y',50, 1, "");
INSERT INTO source VALUES (2400, 'wormpep_id', 1, 'Y', 50, 1, "");
INSERT INTO source VALUES (2410, 'wormbase_gene', 1, 'N',50, 1, "");
INSERT INTO source VALUES (2420, 'wormbase_transcript', 1, 'N', 50, 1, "");
INSERT INTO source VALUES (2430, 'wormbase_locus', 1, 'N', 50, 1, "");
INSERT INTO source VALUES (2440, 'wormbase_pseudogene', 1, 'N', 50, 1, "");

# drosphila melanogster sources 
INSERT INTO source VALUES (2500, 'flybase_gff', 1, 'Y', 50, 1, "");
INSERT INTO source VALUES (2510, 'flybase_gene_id', 1, 'N', 50, 1, "");
INSERT INTO source VALUES (2520, 'flybase_transcript_id', 1, 'N', 50, 1, "");
INSERT INTO source VALUES (2530, 'flybase_polypeptide_id', 1, 'N', 50, 1, "");
INSERT INTO source VALUES (2540, 'flybase_annotation_id', 1, 'N', 50, 1, "");
INSERT INTO source VALUES (2550, 'flybase_synonym', 1, 'N', 50, 1, "");
INSERT INTO source VALUES (2560, 'flybase_name', 1, 'N', 50, 1, "");
INSERT INTO source VALUES (2561, 'FlyBaseName_gene', 1, 'N', 50, 1, "");
INSERT INTO source VALUES (2562, 'FlyBaseName_transcript', 1, 'N', 50, 1, "");
INSERT INTO source VALUES (2563, 'FlyBaseName_translations', 1, 'N', 50, 1, "");
INSERT INTO source VALUES (2570, 'gadfly_gene_cgid', 1, 'N', 50, 1, "");
INSERT INTO source VALUES (2571, 'gadfly_transcript_cgid', 1, 'N', 50, 1, "");
INSERT INTO source VALUES (2572, 'gadfly_translation_cgid', 1, 'N', 50, 1, "");

# ciona intestinalis sources
INSERT INTO source VALUES (2601, 'cint_jgi_v1', 1, 'Y', 50, 1, "");
INSERT INTO source VALUES (2602, 'cint_jgi_v2', 1, 'Y', 50, 1, "");
INSERT INTO source VALUES (2610, 'cint_aniseed_jgi_v1', 1, 'Y', 50, 1, "");
INSERT INTO source VALUES (2611, 'cint_aniseed_jgi_v2', 1, 'Y', 50, 1, "");

# anopheles gambiae sources
# predicted versions of Uniprot/EMBL - created as dependent in appropriate parse
INSERT INTO source VALUES (2857, "ARRAY_EMBL_MMC1_20k_v1", 1, 'Y', 40,1, "" );
INSERT INTO source VALUES (2858, "ARRAY_EMBL_MMC2_12k_v1", 1, 'Y', 40, 1, "");
INSERT INTO source VALUES (2859, "ARRAY_JHU_GAM3_0_21k_v1", 1, 'Y', 40, 1, "");
INSERT INTO source VALUES (2860, "ARRAY_LIV_GAMDETOX_0_25k_v1", 1, 'Y', 40, 1, "");
INSERT INTO source VALUES (2801, "Uniprot/SWISSPROT_predicted", 1, 'N',20, 1, "");
INSERT INTO source VALUES (2802, "Uniprot/SPTREMBL_predicted", 1, 'N',20, 1, "");
INSERT INTO source VALUES (2810, "EMBL_predicted", 1, 'N',20, 1, "");
INSERT INTO source VALUES (2840, 'protein_id_predicted', 1, 'N', 30, 1, "");
INSERT INTO source VALUES (2850, "Celera_Gene", 1, 'Y', 40, 1, "");
INSERT INTO source VALUES (2852, "Celera_Pep", 1, 'Y', 40, 1, "");
INSERT INTO source VALUES (2854, "Celera_Trans", 1, 'Y', 40, 1, "");
INSERT INTO source VALUES (2856, "Anopheles_symbol", 1, 'Y', 40, 1, "");

# aedes aegypti
INSERT INTO source VALUES (7159, 'AedesGenBank', '1','Y','30', 1, "") ;
INSERT INTO source VALUES (2861, "ARRAY_JHU_AEG1.0_20k_v1", 1, 'Y', 40, 1, "");
INSERT INTO source VALUES (2862, "ARRAY_LIV_AEGDETOX_0_25k_v1", 1, 'Y', 40, 1, "");
INSERT INTO source VALUES (2863, "ARRAY_ND_TIGRTC_9_6K_v1", 1, 'Y', 40, 1, "");
INSERT INTO source VALUES (2864, "ARRAY_UCR_GillMgMT_0_2K_v2", 1, 'Y', 40, 1, "");

# Uniprot alternative splice
INSERT INTO source VALUES (3000, "Uniprot/Varsplic", 1, 'Y',20, 1, "");

# Xenopus Jamboree peptides
INSERT INTO source VALUES (3100, "Xenopus_Jamboree", 1, 'Y',20, 1, "");

# ncRNA sources
INSERT INTO source VALUES (4000, "ncRNA", 1, 'Y',70, 1, "");
INSERT INTO source VALUES (4010, "RFAM", 1, 'N',70, 1, "");
INSERT INTO source VALUES (4020, "miRBase", 1, 'N',70, 1, "");
INSERT INTO source VALUES (4025, "miRBase_predicted", 1, 'N',70, 1, "");

# temporary (?) source for mapping Havana OTT transcripts to Ensembl transcripts
INSERT INTO source VALUES (5000, 'OTTT', 1, 'Y', 50, 1, "");

# Illumina
INSERT INTO source VALUES (5011, 'Illumina_V1', 1, 'Y', 50, 1, "");
INSERT INTO source VALUES (5012, 'Illumina_V2', 1, 'Y', 50, 1, "");

# Codelink
INSERT INTO source VALUES (5020, 'Codelink', 1, 'Y', 50, 1, "");

# IMGT sources for Ig/TCR segments
INSERT into source VALUES(4100, 'IMGT/GENE-DB', '200641-6', 'Y', 40, 1, "");
INSERT into source VALUES(4101, 'IMGT/LIGM-DB', '200641-6', 'N', 10, 1, "");



################################################################################
# Files to fetch data from


# --------------------------------------------------------------------------------
### DROSOPHILA TESTING GFF

#enrtexgene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 7227, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");


# chr 2L	
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gff/dmel-2L-r4.2.1.gff.gz', now(), now(), "Flybase_dmel_GFFv3_Parser");

# chr 2R
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gff/dmel-2R-r4.2.1.gff.gz', now(), now(), "Flybase_dmel_GFFv3_Parser");

# chr 3L
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gff/dmel-3L-r4.2.1.gff.gz', now(), now(), "Flybase_dmel_GFFv3_Parser");

# chr 3R
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gff/dmel-3R-r4.2.1.gff.gz', now(), now(), "Flybase_dmel_GFFv3_Parser");

# chr 4
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gff/dmel-4-r4.2.1.gff.gz', now(), now(), "Flybase_dmel_GFFv3_Parser");

# chr X 
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gff/dmel-X-r4.2.1.gff.gz', now(), now(), "Flybase_dmel_GFFv3_Parser");

## HETERO-CHROMATIN

# 2h 
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current_hetchr/gff/dmel-2h-hetr32b2.gff.gz', now(), now(), "Flybase_dmel_GFFv3_Parser");

# 3h 
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current_hetchr/gff/dmel-3h-hetr32b2.gff.gz', now(), now(), "Flybase_dmel_GFFv3_Parser");

# 4h 
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current_hetchr/gff/dmel-4h-hetr32b2.gff.gz', now(), now(), "Flybase_dmel_GFFv3_Parser");

# U 
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current_hetchr/gff/dmel-U-hetr32b2.gff.gz', now(), now(), "Flybase_dmel_GFFv3_Parser");

# Xh 
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current_hetchr/gff/dmel-Xh-hetr32b2.gff.gz', now(), now(), "Flybase_dmel_GFFv3_Parser");

# Yh	 
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (2500, 7227, 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current_hetchr/gff/dmel-Yh-hetr32b2.gff.gz', now(), now(), "Flybase_dmel_GFFv3_Parser");


# Uniprot for drosophila
 
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1,7227,'ftp://ftp.ebi.ac.uk/pub/databases/integr8/uniprot/proteomes/17.D_melanogaster.dat.gz',now(),now(),'UniProtParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1,7227,'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz',now(),now(),'UniProtParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (2,7227,'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz',now(),now(),'UniProtParser');

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 7227, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

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

#        EntrzGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1110, 9606, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

##       uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 9606, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', '', now(), now(), "UniProtParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 9606, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', '', now(), now(), "UniProtParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 9606, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 9606,'ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4, 9606,'ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.rna.fna.gz', '', now(), now(), "RefSeqParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (14, 9606,'LOCAL:RefSeq_dna/CCDS_to_REFSEQ', '', now(), now(), "RefSeq_CCDSParser");


##       GO
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1070, 9606,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/gene_association.goa_human.gz', '', now(), now(), "GOParser");

##       HUGO
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1090, 9606,'LOCAL:HUGO/HUGO_TO_ENSG', '', now(), now(), "HUGO_ENSGParser");

##       HUGO
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1093, 9606,'LOCAL:HUGO/CCDS_TO_HUGO', '', now(), now(), "HUGO_CCDSParser");

# lower priority HUGO data  1091-> refseq  1092->uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1091, 9606,'http://www.gene.ucl.ac.uk/cgi-bin/nomenclature/gdlw.pl?title=Genew+output+data&col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_prev_sym&col=gd_aliases&col=md_prot_id&col=gd_pub_refseq_ids&col=md_eg_id&status=Approved&status=Approved+Non-Human&status_opt=3&=on&where=&order_by=gd_hgnc_id&limit=&format=text&submit=submit&.cgifields=&.cgifields=status&.cgifields=chr', '', now(), now(), "HUGOParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1300, 9606,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', '', now(), now(), "InterproParser");

# mim data obtained from reseq and uniprot files
###      MIM 
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1020, 9606,'ftp://ftp.ncbi.nih.gov/repository/OMIM/omim.txt.Z', '', now(), now(), "MIMParser");

##      IPI
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (5, 9606,'ftp://ftp.ebi.ac.uk/pub/databases/IPI/current/ipi.HUMAN.fasta.gz', '', now(), now(), "IPIParser");

##      CCDS
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2000, 9606,'LOCAL:CCDS/CCDS.txt', '', now(), now(), "CCDSParser");


##      Agilent ( zipfile in ensembl-personal/ianl ) 
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2700, 9606,'LOCAL:AgilentProbe/HumanExpression.fasta', '', now(), now(), "AgilentParser");
#INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2701, 9606,'LOCAL:AgilentCGH/HumanCGH.fasta', '', now(), now(), "AgilentParser");

##      UniGene
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (6, 9606,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Homo_sapiens/Hs.seq.uniq.gz ftp://ftp.ncbi.nih.gov/repository/UniGene/Homo_sapiens/Hs.data.gz', '', now(), now(), "UniGeneParser");

##      ncRNAs presently inhouse ( points to file with dumped xrefs from transfer_ncRNA.pl) .
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4000, 9606,'LOCAL:ncRNA/ncRNA.txt', '', now(), now(), "ncRNAParser");


##      OTTT
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (5000, 9606,'LOCAL:OTTT/OTTT.txt', '', now(), now(), "OTTTParser");

## Illumina version 1
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (5011, 9606, 'LOCAL:Illumina_V1/Human_WG-6.fa', '', now(), now(), "FastaParser");
## Illumina version 2
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (5012, 9606, 'LOCAL:Illumina_V2/Human-6_v2_sequence.csv', '', now(), now(), "IlluminaParser");

# Codelink
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (5020, 9606,'LOCAL:Codelink/GEHC_Amersham_CodeLink_HWG_300026.fna', '', now(), now(), "CodelinkParser");

##      TcR/Ig gene segments presently inhouse ( points to file with dumped direct xrefs )
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4100, 9606,'LOCAL:Segments/Seg_xrefs.txt', '', now(), now(), "SegmentParser");

# --------------------------------------------------------------------------------

###MOUSE
##      EmtrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1110, 10090, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

##      uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 10090, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', '', now(), now(), "UniProtParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 10090, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', '', now(), now(), "UniProtParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 10090, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 10090,'ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4, 10090,'ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.rna.fna.gz', '', now(), now(), "RefSeqParser");

##      mgd (MGI -- MarkerSymbol)
# need to add  ftp://ftp.informatics.jax.org/pub/reports/MRK_Sequence.rpt
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1080, 10090,'ftp://ftp.informatics.jax.org/pub/reports/MRK_SwissProt_TrEMBL.rpt ftp://ftp.informatics.jax.org/pub/reports/MRK_Synonym.sql.rpt', '', now(), now(), "MGDParser");

##      GO 
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1070, 10090,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/MOUSE/gene_association.goa_mouse.gz', '', now(), now(), "GOParser");
# Jax file
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1070, 10090,'http://www.geneontology.org/cgi-bin/downloadGOGA.pl/gene_association.mgi.gz', '', now(), now(), "GOParser");


##      IPI
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (5, 10090,'ftp://ftp.ebi.ac.uk/pub/databases/IPI/current/ipi.MOUSE.fasta.gz', '', now(), now(), "IPIParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1300, 10090,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', '', now(), now(), "InterproParser");

##      UniGene
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (6, 10090,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Mus_musculus/Mm.seq.uniq.gz ftp://ftp.ncbi.nih.gov/repository/UniGene/Mus_musculus/Mm.data.gz', '', now(), now(), "UniGeneParser");

##      Agilent ( zipfile in ensembl-personal/ianl ) 
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2700, 10090,'LOCAL:AgilentProbe/MouseExpression.fasta', '', now(), now(), "AgilentParser");

##      ncRNAs presently inhouse ( points to file with dumped xrefs from transfer_ncRNA.pl) .
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4000, 10090,'LOCAL:ncRNA/ncRNA.txt', '', now(), now(), "ncRNAParser");

# Illumina V1
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (5011, 10090,'LOCAL:Illumina_V1/Mouse-6_V1.fa', '', now(), now(), "FastaParser");

# Codelink
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (5020, 10090,'LOCAL:Codelink/GEHC_Amersham_CodeLink_MWG_300033.fna', '', now(), now(), "CodelinkParser");

# TcR/Ig gene segments presently inhouse ( points to file with dumped direct xrefs )
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4100, 10090,'LOCAL:Segments/Seg_xrefs.txt', '', now(), now(), "SegmentParser");



# --------------------------------------------------------------------------------
### RAT
##      EntrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 10116, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

##      uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 10116, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', '', now(), now(), "UniProtParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 10116, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', '', now(), now(), "UniProtParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 10116, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 10116,'ftp://ftp.ncbi.nih.gov/refseq/R_norvegicus/mRNA_Prot/rat.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4, 10116,'ftp://ftp.ncbi.nih.gov/refseq/R_norvegicus/mRNA_Prot/rat.rna.fna.gz', '', now(), now(), "RefSeqParser");

##      GO 
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1070, 10116,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/RAT/gene_association.goa_rat.gz', '', now(), now(), "GOParser");

##  RGD
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1200, 10116,'ftp://rgd.mcw.edu/pub/data_release/GENES', '', now(), now(), "RGDParser");

##  IPI
##      IPI
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (5, 10116,'ftp://ftp.ebi.ac.uk/pub/databases/IPI/current/ipi.RAT.fasta.gz', '', now(), now(), "IPIParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1300, 10116,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', '', now(), now(), "InterproParser");

##      UniGene
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (6, 10116,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Rattus_norvegicus/Rn.seq.uniq.gz ftp://ftp.ncbi.nih.gov/repository/UniGene/Rattus_norvegicus/Rn.data.gz', '', now(), now(), "UniGeneParser");

##      Agilent ( zipfile in ensembl-personal/ianl ) 
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2700, 10116,'LOCAL:AgilentProbe/RatExpression.fasta', '', now(), now(), "AgilentParser");

##      ncRNAs presently inhouse ( points to file with dumped xrefs from transfer_ncRNA.pl) .
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4000, 10116,'LOCAL:ncRNA/ncRNA.txt', '', now(), now(), "ncRNAParser");

# Codelink
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (5020, 10116,'LOCAL:Codelink/GEHC_Amersham_CodeLink_RWG_300031.fna', '', now(), now(), "CodelinkParser");

# --------------------------------------------------------------------------------

### Zebrafish
##      EmtrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 7955, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

##      uniprot

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 7955, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', '', now(), now(), "UniProtParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 7955, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', '', now(), now(), "UniProtParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 7955, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

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

##      Agilent ( original in /ecs4/work5/is1/zv6/oligo/agilent/AGILENT_G2518A.fa ) 
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2700, 7955,'LOCAL:AgilentProbe/ZebrafishExpression.fasta', '', now(), now(), "AgilentParser");

##      ncRNAs presently inhouse ( points to file with dumped xrefs from transfer_ncRNA.pl) .
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4000, 7955,'LOCAL:ncRNA/ncRNA.txt', '', now(), now(), "ncRNAParser");

# --------------------------------------------------------------------------------
### Chicken
##      EntrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 9031, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

##      uniprot

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 9031, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', '', now(), now(), "UniProtParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 9031, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', '', now(), now(), "UniProtParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 9031, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

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

##      ncRNAs presently inhouse ( points to file with dumped xrefs from transfer_ncRNA.pl) .
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4000, 9031,'LOCAL:ncRNA/ncRNA.txt', '', now(), now(), "ncRNAParser");


# --------------------------------------------------------------------------------
### XENOPUS
##      EmtrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 8364, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

# uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 8364, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', '', now(), now(), "UniProtParser");

#uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 8364, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', '', now(), now(), "UniProtParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 8364, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");


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

##      Xenopus Jamboree
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3100, 8364,'LOCAL:/ecs4/work4/mc2/Xenopus/jamboree/names/names_seq.fa', '', now(), now(), "XenopusJamboreeParser");

##      ncRNAs presently inhouse ( points to file with dumped xrefs from transfer_ncRNA.pl) .
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4000, 8364,'LOCAL:ncRNA/ncRNA.txt', '', now(), now(), "ncRNAParser");

# --------------------------------------------------------------------------------

### Dog

##      EntrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 9615, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

#        uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 9615, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', '', now(), now(), "UniProtParser");

#        uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 9615, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', '', now(), now(), "UniProtParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 9615, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

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

##      ncRNAs presently inhouse ( points to file with dumped xrefs from transfer_ncRNA.pl) .
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4000, 9615,'LOCAL:ncRNA/ncRNA.txt', '', now(), now(), "ncRNAParser");

# --------------------------------------------------------------------------------

### C elegans
#       EntrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 6239, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

#        uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 6239, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', '', now(), now(), "UniProtParser");

#        uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 6239, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', '', now(), now(), "UniProtParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 6239, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

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

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2400, 6239, 'ftp://ftp.sanger.ac.uk/pub/databases/wormpep/wormpep160/wormpep.table160', '', now(), now(), "WormPepParser");

##   Stable ID xref transfer - note use of WormbaseDatabaseStableIDParser
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2400, 6239, 'mysql:ensembldb.ensembl.org::caenorhabditis_elegans_core_39_150a:anonymous', '', now(), now(), "WormbaseDatabaseStableIDParser");



# --------------------------------------------------------------------------------

#### HoneyBee
##      EntrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 7460, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

## Uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 7460, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', '', now(), now(), "UniProtParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 7460, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', '', now(), now(), "UniProtParser");

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 7460, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

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


#### SeaSquirt Ciona intestinalis
##      EntrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 7719, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

## Uniprot
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (1, 7719, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (2, 7719, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 7719, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

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
VALUES (2601, 7719,'ftp://ftp.jgi-psf.org/pub/JGI_data/Ciona/v1.0/ciona.prot.fasta.gz',now(),now(),'JGI_ProteinParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (2602, 7719,'ftp://ftp.jgi-psf.org/pub/JGI_data/Ciona/v2.0/FM1.aa.fasta.gz',now(),now(),'JGI_ProteinParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (2610, 7719,'ftp://ftp.jgi-psf.org/pub/JGI_data/Ciona/v1.0/ciona.prot.fasta.gz',now(),now(),'JGI_ProteinParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (2611, 7719,'ftp://ftp.jgi-psf.org/pub/JGI_data/Ciona/v2.0/FM1.aa.fasta.gz',now(),now(),'JGI_ProteinParser');

##      ncRNAs presently inhouse ( points to file with dumped xrefs from transfer_ncRNA.pl) .
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4000, 7719,'LOCAL:ncRNA/ncRNA.txt', '', now(), now(), "ncRNAParser");



#### Ciona savignyi 

##      EntrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 51511, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

## Uniprot
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (1, 51511, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (2, 51511, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 51511, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

##       refseq ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/ *gpff
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 51511,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other1.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 51511,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other2.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 51511,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other3.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 51511,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other4.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 51511,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other5.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 51511,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other6.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 51511,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other7.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser"); 

##       refseq ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/ *rna*fna*
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 51511,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other1.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 51511,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other2.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 51511,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other3.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 51511,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other4.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 51511,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other5.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 51511,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other6.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 51511,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other7.rna.fna.gz',  now(), now(), "RefSeqParser"); 


INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,51511,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate1.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,51511,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate2.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,51511,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate3.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,51511,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate4.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,51511,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate5.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,51511,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate6.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,51511,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate7.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,51511,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate8.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,51511,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate9.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,51511,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate10.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,51511,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate11protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,51511,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate12.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,51511,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate13.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,51511,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate14.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (3,51511,'ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate15.protein.gpff.gz',now(),now(),'RefSeqGPFFParser');


##       GO
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(1070, 51511,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz', now(), now(), "GOParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(1300, 51511,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', now(), now(), "InterproParser");

##      UniGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(6, 51511,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Ciona_savignyi/Csa.seq.uniq.gz ftp://ftp.ncbi.nih.gov/repository/UniGene/Ciona_savignyi/Csa.data.gz',  now(), now(), "UniGeneParser");


##      ncRNAs presently inhouse ( points to file with dumped xrefs from transfer_ncRNA.pl) .
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4000, 51511,'LOCAL:ncRNA/ncRNA.txt', '', now(), now(), "ncRNAParser");



# --------------------------------------------------------------------------------

###### Anopheles

##      EntrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 7165, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

# Uniprot
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1,7165,'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz',now(),now(),'UniProtParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (2,7165,'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz',now(),now(),'UniProtParser');

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 7165, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

# Interpro
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1300,7165,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz',now(),now(),'InterproParser');

# Unigene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (6,7165,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Anopheles_gambiae/Aga.seq.uniq.gz  ftp://ftp.ncbi.nih.gov/repository/UniGene/Anopheles_gambiae/Aga.data.gz',now(),now(),'UniGeneParser');

# Celera
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (2852, 7165,'LOCAL:Celera_Pep/consensus-proteins_xref.fasta', now(), now(), "CeleraProteinParser");

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (2854, 7165,'LOCAL:Celera_trans/consensus-transcripts_xref.fasta', now(), now(), "CeleraTranscriptParser");

# Anopheles symbols
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (2856, 7165,'LOCAL:Anopheles_symbol/GeneName_translation_UniqID.txt', now(), now(), "AnophelesSymbolParser");

# GO
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1070, 7165,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz', now(), now(), "GOParser");

# VectorBase Arrays
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (2857, 7165,'LOCAL:Anopheles_gambiae/EMBL/EMBL_MMC1_20k_v1.txt', now(), now(), "VbDirectParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (2858, 7165,'LOCAL:Anopheles_gambiae/EMBL/EMBL_MMC2_12k_v1.txt', now(), now(), "VbDirectParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (2859, 7165,'LOCAL:Anopheles_gambiae/JHU/JHU_GAM3.0_21k_v1.txt', now(), now(), "VbDirectParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (2860, 7165,'LOCAL:Anopheles_gambiae/LIV/LIV_GAMDETOX_0.25k_v1.txt', now(), now(), "VbDirectParser");

###### Aedes 

##      EntrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 7159, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

# Uniprot
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1,7159,'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz',now(),now(),'UniProtParser');

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (2,7159,'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz',now(),now(),'UniProtParser');

INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 7159, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

# Interpro
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1300,7159,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz',now(),now(),'InterproParser');

# Unigene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (6,7159,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Aedes_aegypti/Aae.seq.uniq.gz ftp://ftp.ncbi.nih.gov/repository/UniGene/Aedes_aegypti/Aae.data.gz',now(),now(),'UniGeneParser');

# Aedes GenBank
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (7159, 7159,'LOCAL:AedesGenBank/Aedes_proteinID.fa', now(), now(), "AedesGenBankParser");

# GO
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1070, 7159,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz', now(), now(), "GOParser");

# Vectorbase Arrays
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (2861, 7159,'LOCAL:Aedes_aegypti/JHU/JHU_AEG1.0_20k_v1.txt', now(), now(), "VbDirectParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (2862, 7159,'LOCAL:Aedes_aegypti/LIV/LIV_AEGDETOX_0.25k_v1.txt', now(), now(), "VbDirectParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (2863, 7159,'LOCAL:Aedes_aegypti/ND/ND_TIGRTC_9.6k_v1.txt', now(), now(), "VbDirectParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (2864, 7159,'LOCAL:Aedes_aegypti/UCR/UCR_GillMgMT_0.2k_v2.txt', now(), now(), "VbDirectParser");


# --------------------------------------------------------------------------------


#### Cow

##      EntrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 9913, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

## Uniprot
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (1, 9913, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (2, 9913, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 9913, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

##       refseq protein
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 9913,'ftp://ftp.ncbi.nih.gov/genomes/Bos_taurus/protein/protein.gbk.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4, 9913,'ftp://ftp.ncbi.nih.gov/genomes/Bos_taurus/RNA/rna.gbk.gz', '', now(), now(), "RefSeqGPFFParser");

##       GO
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1070, 9913,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz', now(), now(), "GOParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1300, 9913,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', now(), now(), "InterproParser");

##      UniGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (6, 9913,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Bos_taurus/Bt.seq.uniq.gz ftp://ftp.ncbi.nih.gov/repository/UniGene/Bos_taurus/Bt.data.gz',  now(), now(), "UniGeneParser");

##      ncRNAs presently inhouse ( points to file with dumped xrefs from transfer_ncRNA.pl) .
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4000, 9913,'LOCAL:ncRNA/ncRNA.txt', '', now(), now(), "ncRNAParser");

# -----------------------------------------------------------------------------------

#### Monodelphis

##      EntrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 13616, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

## Uniprot
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (1, 13616, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (2, 13616, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 13616, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian1.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian2.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian3.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian4.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian5.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian6.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian7.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian8.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian9.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian10.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian11.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian12.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian13.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian14.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian15.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian16.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian17.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian18.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian19.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian20.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian21.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian22.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian23.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian24.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian25.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian26.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian27.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian28.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian29.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian30.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian31.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian32.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian33.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian34.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian35.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian36.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian37.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian38.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian39.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian40.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian41.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian42.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian43.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian44.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian45.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian46.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian47.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian48.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian49.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian50.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian51.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian52.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian53.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian54.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian55.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian56.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian57.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian58.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian59.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian60.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian61.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian62.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian1.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian2.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian3.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian4.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian5.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian6.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian7.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian8.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian9.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian10.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian11.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian12.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian13.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian14.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian15.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian16.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian17.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian18.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian19.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian20.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian21.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian22.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian23.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian24.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian25.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian26.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian27.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian28.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian29.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian30.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian31.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian32.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian33.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian34.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian35.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian36.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian37.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian38.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian39.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian40.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian41.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian42.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian43.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian44.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian45.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian46.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian47.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian48.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian49.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian50.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian51.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian52.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian53.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian54.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian55.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian56.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian57.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian58.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian59.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian60.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian61.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 13616,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian62.rna.fna.gz',  now(), now(), "RefSeqParser");

##       GO
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(1070, 13616,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz', now(), now(), "GOParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(1300, 13616,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', now(), now(), "InterproParser");

##      ncRNAs presently inhouse ( points to file with dumped xrefs from transfer_ncRNA.pl) .
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4000, 13616,'LOCAL:ncRNA/ncRNA.txt', '', now(), now(), "ncRNAParser");


# ------------------------------------------------------------------------------
#### Fugu

##      EntrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 31033, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

## Uniprot
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (1, 31033, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (2, 31033, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 31033, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 31033,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other1.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 31033,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other2.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 31033,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other3.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 31033,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other4.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 31033,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other5.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 31033,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other6.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 31033,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other7.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 31033,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other1.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 31033,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other2.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 31033,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other3.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 31033,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other4.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 31033,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other5.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 31033,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other6.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 31033,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other7.rna.fna.gz',  now(), now(), "RefSeqParser");

##       GO
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(1070, 31033,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz', now(), now(), "GOParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(1300, 31033,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', now(), now(), "InterproParser");


##      UniGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(6, 31033,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Takifugu_rubripes/Tru.seq.uniq.gz ftp://ftp.ncbi.nih.gov/repository/UniGene/Takifugu_rubripes/Tru.data.gz',  now(), now(), "UniGeneParser");

##      ncRNAs presently inhouse ( points to file with dumped xrefs from transfer_ncRNA.pl) .
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4000, 31033,'LOCAL:ncRNA/ncRNA.txt', '', now(), now(), "ncRNAParser");


# -----------------------------------------------------------------------------------
#### Yeast

##      EntrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 4932, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1250, 4932, 'ftp://genome-ftp.stanford.edu/pub/yeast/gene_registry/registry.genenames.tab', now(), now(), "SGDParser");



## Uniprot
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (1, 4932, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (2, 4932, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 4932, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi1.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi2.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi3.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi4.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi5.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi6.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi7.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi8.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi9.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi10.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi11.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi12.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi1.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi2.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi3.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi4.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi5.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi6.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi7.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi8.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi9.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi10.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi11.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 4932,'ftp://ftp.ncbi.nih.gov/refseq/release/fungi/fungi12.rna.fna.gz',  now(), now(), "RefSeqParser");

##       GO
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(1070, 4932,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz', now(), now(), "GOParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(1300, 4932,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', now(), now(), "InterproParser");



####rhesus macaque

##      EntrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 9544, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

## Uniprot
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (1, 9544, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (2, 9544, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', now(), now(), "UniProtParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian1.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian2.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian3.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian4.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian5.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian6.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian7.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian8.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian9.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian10.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian11.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian12.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian13.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian14.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian15.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian16.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian17.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian18.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian19.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian20.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian21.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian22.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian23.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian24.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian25.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian26.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian27.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian28.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian29.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian30.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian31.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian32.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian33.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian34.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian35.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian36.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian37.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian38.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian39.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian40.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian41.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian42.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian43.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian44.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian45.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian46.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian1.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian2.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian3.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian4.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian5.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian6.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian7.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian8.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian9.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian10.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian11.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian12.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian13.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian14.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian15.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian16.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian17.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian18.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian19.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian20.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian21.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian22.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian23.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian24.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian25.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian26.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian27.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian28.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian29.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian30.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian31.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian32.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian33.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian34.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian35.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian36.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian37.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian38.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian39.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian40.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian41.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian42.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian43.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian44.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian45.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9544,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian46.rna.fna.gz',  now(), now(), "RefSeqParser");

##       GO
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(1070, 9544,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz', now(), now(), "GOParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(1300, 9544,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', now(), now(), "InterproParser");

##      ncRNAs presently inhouse ( points to file with dumped xrefs from transfer_ncRNA.pl) .
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4000, 9544,'LOCAL:ncRNA/ncRNA.txt', '', now(), now(), "ncRNAParser");


##------------------------------------------------------------------------------
##### elephant

##      EntrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 9785, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

## Uniprot
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (1, 9785, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (2, 9785, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 9785, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

##       GO
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1070, 9785,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz', now(), now(), "GOParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1300, 9785,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', now(), now(), "InterproParser");

##      ncRNAs presently inhouse ( points to file with dumped xrefs from transfer_ncRNA.pl) .
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4000, 9785,'LOCAL:ncRNA/ncRNA.txt', '', now(), now(), "ncRNAParser");


##------------------------------------------------------------------------------
##### armadillo

##      EntrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 9361, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

## Uniprot
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (1, 9361, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (2, 9361, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 9361, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

##       GO
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1070, 9361,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz', now(), now(), "GOParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1300, 9361,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', now(), now(), "InterproParser");


##      ncRNAs presently inhouse ( points to file with dumped xrefs from transfer_ncRNA.pl) .
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4000, 9361,'LOCAL:ncRNA/ncRNA.txt', '', now(), now(), "ncRNAParser");


##------------------------------------------------------------------------------
##### rabbit

##      EntrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 9986, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

## Uniprot
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (1, 9986, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (2, 9986, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 9986, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

##       GO
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1070, 9986,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz', now(), now(), "GOParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1300, 9986,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', now(), now(), "InterproParser");


##      ncRNAs presently inhouse ( points to file with dumped xrefs from transfer_ncRNA.pl) .
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4000, 9986,'LOCAL:ncRNA/ncRNA.txt', '', now(), now(), "ncRNAParser");

##------------------------------------------------------------------------------
##### tenrec

##      EntrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 9371, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

## Uniprot
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (1, 9371, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (2, 9371, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 9371, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

##       GO
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1070, 9371,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz', now(), now(), "GOParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1300, 9371,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', now(), now(), "InterproParser");

##------------------------------------------------------------------------------
##### cat

##      EntrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 9685, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

## Uniprot
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (1, 9685, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (2, 9685, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 9685, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

##       GO
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1070, 9685,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz', now(), now(), "GOParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1300, 9685,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', now(), now(), "InterproParser");

##------------------------------------------------------------------------------
### STICKLEBACK

##      EntrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 69293, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

##       uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1,69293 , 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', '', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2,69293 , 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', '', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000,69293 , 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");
##       GO
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1070,69293 ,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz', now(), now(), "GOParser");
##      Interpro
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (1300,69293 ,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', now(), now(), "InterproParser");
# refseq
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (3,69293 ,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other1.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (3,69293 ,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other2.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (3,69293 ,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other3.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (3,69293 ,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other4.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (3,69293 ,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other5.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (3,69293 ,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other6.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (3,69293 ,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other7.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (4, 69293,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other1.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (4, 69293,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other2.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (4, 69293,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other3.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (4, 69293,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other4.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (4, 69293,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other5.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (4, 69293,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other6.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (4, 69293,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other7.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES (4, 69293,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/sequences.fasta.gz',  now(), now(), "RefSeqParser");

##      ncRNAs presently inhouse ( points to file with dumped xrefs from transfer_ncRNA.pl) .
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4000, 69293,'LOCAL:ncRNA/ncRNA.txt', '', now(), now(), "ncRNAParser");


##------------------------------------------------------------------------------
##### medaka

##      EntrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 8090, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

## Uniprot
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (1, 8090, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (2, 8090, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 8090, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 8090,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other1.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 8090,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other2.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 8090,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other3.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 8090,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other4.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 8090,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other5.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 8090,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other6.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 8090,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other7.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 8090,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other1.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 8090,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other2.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 8090,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other3.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 8090,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other4.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 8090,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other5.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 8090,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other6.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (4, 8090,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_other7.rna.fna.gz',  now(), now(), "RefSeqParser");

##       GO
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(1070, 8090,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz', now(), now(), "GOParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(1300, 8090,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', now(), now(), "InterproParser");


##      UniGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(6, 8090,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Oryzias_latipes/Ola.seq.uniq.gz ftp://ftp.ncbi.nih.gov/repository/UniGene/Oryzias_latipes/Ola.data.gz',  now(), now(), "UniGeneParser");


##-----------------------------------------------
##### microbat

##      EntrezGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser)\
VALUES (1110, 59463, 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz', now(), now(), "EntrezGeneParser");

## Uniprot
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (1, 59463, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (2, 59463, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 59463, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian1.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian10.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian11.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian12.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian13.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian14.protein.gpff.gz", now() , now() , "RefSeqGPFFParser" ) ;
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian15.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian16.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian17.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian18.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian19.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian2.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian20.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian21.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian22.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian23.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian24.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian25.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian26.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian27.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian28.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian29.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian3.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian30.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian31.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian32.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian33.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian34.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian35.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian36.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian37.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian38.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian39.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian4.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian40.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian41.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian42.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian43.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian44.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian45.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian46.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian47.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian48.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian49.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian5.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian6.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian7.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian8.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian9.protein.gpff.gz", now(), now(), "RefSeqGPFFParser" );


INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian1.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian10.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian11.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian12.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian13.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian14.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian15.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian16.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian17.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian18.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian19.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian2.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian20.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian21.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian22.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian23.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian24.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian25.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian26.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian27.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian28.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian29.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian3.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian30.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian31.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian32.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian33.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian34.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian35.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian36.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian37.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian38.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian39.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian4.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian40.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian41.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian42.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian43.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian44.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian45.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian46.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian47.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian48.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian49.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian5.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian6.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian7.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian8.rna.fna.gz", now(), now(), "RefSeqParser" );
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES \
  (4, 59463, "ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_other/vertebrate_mammalian9.rna.fna.gz", now(), now(), "RefSeqParser" );

##       GO
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(1070, 59463,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz', now(), now(), "GOParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(1300, 59463,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', now(), now(), "InterproParser");


##      UniGene
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(6, 59463,'ftp://ftp.ncbi.nih.gov/repository/UniGene/Myotis_lucifugus/Mlu.seq.uniq.gz ftp://ftp.ncbi.nih.gov/repository/UniGene/Myotis_lucifugus/Mlu.data.gz',  now(), now(), "UniGeneParser");
##Note: this UniGene entry does not yet exists

##      ncRNAs presently inhouse ( points to file with dumped xrefs from transfer_ncRNA.pl) .
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4000, 59463,'LOCAL:ncRNA/ncRNA.txt', '', now(), now(), "ncRNAParser");


##-----------------------------------------------
#### Platypus

## Uniprot
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (1, 9258, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
 (2, 9258, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.dat.gz', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3000, 9258, 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz', '', now(), now(), "UniProtVarSplicParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian1.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian2.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian3.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian4.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian5.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian6.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian7.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian8.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian9.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian10.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian11.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian12.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian13.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian14.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian15.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian16.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian17.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian18.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian19.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian20.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian21.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian22.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian23.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian24.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian25.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian26.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian27.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian28.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian29.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian30.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian31.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian32.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian33.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian34.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian35.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian36.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian37.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian38.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian39.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian40.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian41.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian42.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian43.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian44.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian45.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(3, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian46.protein.gpff.gz',  now(), now(), "RefSeqGPFFParser");

INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian1.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian2.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian3.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian4.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian5.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian6.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian7.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian8.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian9.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian10.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian11.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian12.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian13.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian14.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian15.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian16.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian17.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian18.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian19.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian20.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian21.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian22.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian23.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian24.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian25.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian26.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian27.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian28.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian29.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian30.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian31.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian32.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian33.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian34.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian35.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian36.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian37.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian38.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian39.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian40.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian41.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian42.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian43.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian44.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian45.rna.fna.gz',  now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(4, 9258,'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian46.rna.fna.gz',  now(), now(), "RefSeqParser");

##       GO
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(1070, 9258,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz', now(), now(), "GOParser");

##      Interpro
INSERT INTO source_url (source_id, species_id, url, file_modified_date, upload_date, parser) VALUES\
(1300, 9258,'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', now(), now(), "InterproParser");



##      ncRNAs presently inhouse ( points to file with dumped xrefs from transfer_ncRNA.pl) .
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (4000, 8090,'LOCAL:ncRNA/ncRNA.txt', '', now(), now(), "ncRNAParser");


##------------------------------------------------------------------------------


