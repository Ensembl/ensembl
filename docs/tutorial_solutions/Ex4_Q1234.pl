# file - Ex4_Q1234.pl
# author - Dan.Andrews@sanger.ac.uk
# desc - answers questions 1, 2, 3 and 4 from Exercise 4 of the EnsEMBL api tutorial.

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use strict;

# Initialise a new database adaptor object.

my $host   = 'kaka.sanger.ac.uk';
my $user   = 'anonymous';
my $dbname = 'current';

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-host   => $host,
					     -user   => $user,
					     -dbname => $dbname);

# Set the static golden path type to a particular golden path assembly.
$db->static_golden_path_type('NCBI_28');

# Derive from the database adaptor a static golden path adaptor.
my $sgp = $db->get_StaticGoldenPathAdaptor;

# Derive a virtual contig by chromosome number, taking the first megabase.
my $vcontig = $sgp->fetch_VirtualContig_by_chr_start_end('1',1,1000000);

# Get an exhaustive list of gene objects from the virtual contig.
my @genes = $vcontig->get_all_Genes;

# Loop through the list of genes, printing some information about each.
foreach my $gene (@genes){

    # Print the ensembl id.
    print $gene->stable_id . "\n";

    # Check whether the gene is a known gene.
    if ($gene->is_known) {
	print "     This is a known gene.  Aliases for this gene from other databases are\:\n";

        # If the gene is known, follow DBLinks...

	my @db_links = $gene->each_DBLink;
	foreach my $db_link (@db_links){

	    # And for each link followed, print the common name for this gene.

	    print "     ". $db_link->display_id . "\n";
	}

    } else {print "     This is a novel gene.\n";} # Print a message if the gene is not previously known.

    # For this gene, grab a list of the transcripts and translate each of them. 
    foreach my $transcript ($gene->each_Transcript){

	# To get a translation, use the translate method of the transcript object.
	# This returns a sequence object and the string representation of the sequence
	# can be derived using the seq method.

	my $peptide = $transcript->translate;
	if($peptide->seq =~ /\*/){print 'Excuse me, but there appears to be a stop codon in your translation.\n'}
	print "          Transcript id " . $transcript->stable_id . "\n";
	print "          Translation " . $peptide->seq . "\n";
    }
}


