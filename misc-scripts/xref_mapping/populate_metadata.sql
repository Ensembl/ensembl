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

INSERT INTO source VALUES (1, "UniProtSwissProt", 1, 'Y');
INSERT INTO source VALUES (2, "UniProtSPTrEMBL", 1, 'Y');
INSERT INTO source VALUES (3, "RefSeq", 1, 'Y');

# Other sources - used to create dependent xrefs, but not to upload from

INSERT INTO source VALUES (1000, 'EMBL', 1, 'N');
INSERT INTO source VALUES (1001, 'PUBMED', 1, 'N');
INSERT INTO source VALUES (1002, 'MEDLINE', 1, 'N');
INSERT INTO source VALUES (1003, 'GO', 1, 'N');

################################################################################
# Files to fetch data from

# --------------------------------------------------------------------------------
# UniProt (SwissProt & SPTrEMBL)

# Note currently no UniProt data for fugu, anopheles, c.briggsae or chicken.
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/9606.SPC', '', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/10090.SPC', '', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/10116.SPC', '', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/7227.SPC', '', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/6239.SPC', '', now(), now(), "UniProtParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/9031.SPC', '', now(), now(), "UniProtParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/9598.SPC', '', now(), now(), "UniProtParser");

# --------------------------------------------------------------------------------
# RefSeq - release/ and cumulative/ directories, for protein and mRNA
## XXXXXXx
# release/rna
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian1.rna.fna.gz', '', now(), now(), "RefSeqParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian3.rna.fna.gz', '', now(), now(), "RefSeqParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian3.rna.fna.gz', '', now(), now(), "RefSeqParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian4.rna.fna.gz', '', now(), now(), "RefSeqParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian5.rna.fna.gz', '', now(), now(), "RefSeqParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian6.rna.fna.gz', '', now(), now(), "RefSeqParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian7.rna.fna.gz', '', now(), now(), "RefSeqParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian8.rna.fna.gz', '', now(), now(), "RefSeqParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian9.rna.fna.gz', '', now(), now(), "RefSeqParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian10.rna.fna.gz', '', now(), now(), "RefSeqParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian11.rna.fna.gz', '', now(), now(), "RefSeqParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian13.rna.fna.gz', '', now(), now(), "RefSeqParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian13.rna.fna.gz', '', now(), now(), "RefSeqParser");
##
### release/protein
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian1.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian3.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian3.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian4.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian5.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian6.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian7.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian8.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian9.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian10.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian11.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian13.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian13.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
##
### cumulative
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/cumulative/rscu.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (3, 'ftp://ftp.ncbi.nih.gov/refseq/cumulative/rscu.fna.gz', '', now(), now(), "RefSeqParser");

################################################################################

