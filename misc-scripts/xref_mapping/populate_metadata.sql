# Populate the appropriate tables in an xref metadata database


################################################################################
# SPECIES

INSERT INTO species (taxonomy_id, name) VALUES (9606,  'homo_sapiens');
INSERT INTO species (taxonomy_id, name) VALUES (10090, 'mus_musculus');
INSERT INTO species (taxonomy_id, name) VALUES (10116, 'rattus_norvegicus');
INSERT INTO species (taxonomy_id, name) VALUES (31033, 'fugu_rubripes');
INSERT INTO species (taxonomy_id, name) VALUES (7165,  'anopheles_gambiae');
INSERT INTO species (taxonomy_id, name) VALUES (7227,  'drosophila_melanogaster');
INSERT INTO species (taxonomy_id, name) VALUES (6239,  'caenorhabditis_elegans');
INSERT INTO species (taxonomy_id, name) VALUES (6238,  'caenorhabditis_briggsae');
INSERT INTO species (taxonomy_id, name) VALUES (7955,  'danio_rerio');
INSERT INTO species (taxonomy_id, name) VALUES (9598,  'pan_troglodytes');
INSERT INTO species (taxonomy_id, name) VALUES (9031,  'gallus_gallus');
INSERT INTO species (taxonomy_id, name) VALUES (99883, 'tetraodon_nigroviridis');
INSERT INTO species (taxonomy_id, name) VALUES (170,   'apis_mellifera');

################################################################################
# SOURCES - types of data we can read

# "High level" sources that we will also download from (via source_url)

INSERT INTO source VALUES (1, "UniProtSwissProt", 1, 'Y');
INSERT INTO source VALUES (2, "RefSeq", 1, 'Y');

# Other sources - used to create dependent xrefs, but not to upload from

INSERT INTO source VALUES (1000, 'EMBL', 1, 'N');
INSERT INTO source VALUES (1001, 'PUBMED', 1, 'N');
INSERT INTO source VALUES (1002, 'MEDLINE', 1, 'N');
INSERT INTO source VALUES (1003, 'GO', 1, 'N');

################################################################################
# Files to fetch data from

# --------------------------------------------------------------------------------
# UniProt/SwissProt
# Note currently no UniProt/SwissProt data for fugu, anopheles, c.briggsae or chicken.
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/9606.SPC', '', now(), now(), "SwissProtParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/10090.SPC', '', now(), now(), "SwissProtParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/10116.SPC', '', now(), now(), "SwissProtParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/7227.SPC', '', now(), now(), "SwissProtParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/6239.SPC', '', now(), now(), "SwissProtParser");
##INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/9031.SPC', '', now(), now(), "SwissProtParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (1, 'ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/9598.SPC', '', now(), now(), "SwissProtParser");

# --------------------------------------------------------------------------------
# RefSeq - release/ and cumulative/ directories, for protein and mRNA

# release/rna
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian1.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian2.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian3.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian4.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian5.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian6.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian7.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian8.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian9.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian10.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian11.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian12.rna.fna.gz', '', now(), now(), "RefSeqParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian13.rna.fna.gz', '', now(), now(), "RefSeqParser");

# release/protein
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian1.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian2.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian3.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian4.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian5.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian6.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian7.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian8.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian9.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian10.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian11.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian12.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian13.protein.gpff.gz', '', now(), now(), "RefSeqGPFFParser");

# cumulative
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/cumulative/rscu.gpff.gz', '', now(), now(), "RefSeqGPFFParser");
INSERT INTO source_url (source_id, url, checksum, file_modified_date, upload_date, parser) VALUES (2, 'ftp://ftp.ncbi.nih.gov/refseq/cumulative/rscu.fna.gz', '', now(), now(), "RefSeqParser");

################################################################################

