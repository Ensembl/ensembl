# Populate the appropriate tables in an xref metadata database


################################################################################
# SPECIES

INSERT INTO species (taxonomy_id, name, aliases) VALUES (9606,  'homo_sapiens',            'human,hsapiens,homosapiens');
INSERT INTO species (taxonomy_id, name, aliases) VALUES (10090, 'mus_musculus',            'mouse,mmusculus,musmusculus');
INSERT INTO species (taxonomy_id, name, aliases) VALUES (10116, 'rattus_norvegicus',       'rat,rnovegicus,rattusnorvegicus');
INSERT INTO species (taxonomy_id, name, aliases) VALUES (31033, 'fugu_rubripes',           'pufferfish,fugu,frubripes,fugurubripes');
INSERT INTO species (taxonomy_id, name, aliases) VALUES (7165,  'anopheles_gambiae',       'mosquito,anopheles,agambiae,anophelesgambiae');
INSERT INTO species (taxonomy_id, name, aliases) VALUES (7227,  'drosophila_melanogaster', 'drosophila,dmelongaster,drosophilamelanogaster' );
INSERT INTO species (taxonomy_id, name, aliases) VALUES (6239,  'caenorhabditis_elegans',  'elegans,celegans,caenorhabditiselegans');
INSERT INTO species (taxonomy_id, name, aliases) VALUES (6238,  'caenorhabditis_briggsae', 'briggsae,cbriggsae,caenorhabditisbriggsae');
INSERT INTO species (taxonomy_id, name, aliases) VALUES (7955,  'danio_rerio',             'zebrafish,danio,drerio,daniorerio' );
INSERT INTO species (taxonomy_id, name, aliases) VALUES (9598,  'pan_troglodytes',         'chimp,chimpanzee,ptroglodytes,pantroglodytes');
INSERT INTO species (taxonomy_id, name, aliases) VALUES (9031,  'gallus_gallus',           'chicken,chick,ggallus,gallusgallus' );
INSERT INTO species (taxonomy_id, name, aliases) VALUES (99883, 'tetraodon_nigroviridis',  'tetraodon,tnigroviridis,tetraodonnigroviridis');
INSERT INTO species (taxonomy_id, name, aliases) VALUES (9913,   'bos_taurus',             'cow,btaurus,bostaurus');
INSERT INTO species (taxonomy_id, name, aliases) VALUES (9615,   'canis_familaris',        'dog,doggy,cfamiliaris,canisfamiliaris');

################################################################################
# SOURCES - types of data we can read

# "High level" sources that we will also download from (via source_url)

INSERT INTO source VALUES (1, "UniProtSwissProt", 1, 'Y',1);
INSERT INTO source VALUES (2, "UniProtSPTrEMBL", 1, 'Y',1);
INSERT INTO source VALUES (3, "RefSeq", 1, 'Y',1);

# Other sources - used to create dependent xrefs, but not to upload from

INSERT INTO source VALUES (1010, 'EMBL', 1, 'N', 2);
INSERT INTO source VALUES (1020, 'MIM', 1, 'N', 2);
INSERT INTO source VALUES (1030, 'PDB', 1, 'N', 2);
INSERT INTO source VALUES (1040, 'protein_id', 1, 'N', 2);
INSERT INTO source VALUES (1050, 'PUBMED', 1, 'N', 2);
INSERT INTO source VALUES (1060, 'MEDLINE', 1, 'N', 2);

INSERT INTO source VALUES (1070, 'GO', 1, 'Y',2);
INSERT INTO source VALUES (1080, 'MarkerSymbol', 1, 'Y',2);
INSERT INTO source VALUES (1090, 'HUGO', 1, 'Y',2);

################################################################################
# Files to fetch data from

# --------------------------------------------------------------------------------
# UniProt (SwissProt & SPTrEMBL)

# Note currently no UniProt data for fugu, anopheles, c.briggsae or chicken.


###HUMAN
##       uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 1,'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/9606.SPC', '', now(), now(), "UniProtParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 1,'ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot/human.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");

##       refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 1,'ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot/human.rna.fna.gz', '', now(), now(), "RefSeqParser");

##       GO
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1070, 1,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/gene_association.goa_human.gz', '', now(), now(), "GOParser");

##       HUGO
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1090, 1,'http://www.gene.ucl.ac.uk/public-files/nomen/ens1.txt', '', now(), now(), "HUGOParser");




###MOUSE
##      uniprot
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 2, 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/10090.SPC', '', now(), now(), "UniProtParser");

##      refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 2,'ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");

##      refseq
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 2,'ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.rna.fna.gz', '', now(), now(), "RefSeqParser");

##      mgd (MGI -- MarkerSymbol)
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1080, 2,'ftp://ftp.informatics.jax.org/pub/reports/MRK_SwissProt_TrEMBL.rpt', '', now(), now(), "MGDParser");

##      GO 
INSERT INTO source_url (source_id, species_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1070, 2,'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/MOUSE/gene_association.goa_mouse.gz', '', now(), now(), "GOParser");

################################################################################

