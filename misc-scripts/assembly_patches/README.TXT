#Get the sequence data. Should be able to down load fasta later but for now do it manually.

mfetch -d embl -v fasta AC093223.3 AC145145.2 AC144587.2 AC142525.2 AC139494.2 AC118459.2 AC139500.2 AC108108.2 AC138827.3 AC146335.1 AC139463.2 AC138832.2 > seq.fasta


Get the mappings

In Genome workbench click on Patch Scaffolds and export agp.
to get hap -> contig

See file thjat was attached to get ref->hap mappings.


NOW use....
-----------

 wget ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/PATCHES/patch_release_1/FASTA/alt.scaf.fa.gz
gunzip alt.scaf.fa.gz

 wget ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/PATCHES/patch_release_1/alt_scaffold_placement.txt

wget ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/PATCHES/patch_release_1/AGP/alt.scaf.agp.gz
gunzip alt.scaf.agp.gz



