# file - Ex3_Q123.pl
# author - Dan.Andrews@sanger.ac.uk
# desc - answers questions 1, 2 and 3 from Exercise 3 of the EnsEMBL api tutorial.

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use strict;

# Once again, use a DBAdaptor object to connect to an EnsEMBL database.
my $host   = 'kaka.sanger.ac.uk';
my $user   = 'anonymous';
my $dbname = 'current';

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $host,
					    -user   => $user,
					    -dbname => $dbname);

# Get all the clone ids.
my @clone_ids = $db->get_all_Clone_id;

my $gene_count = 1;
my $num_genes = 50;
my $first_few_genes = 10;
my $multi_transcripts = 0;
my $exon_counter = 0;
my $gene_length_counter = 0;

# Using the clone ids as a means to retrieve contigs and thus genes, work through until 
# 50 genes and their transcripts have been examined.
CLONE:
foreach my $clone_id (@clone_ids){
    print "Looking at clone $clone_id\n";

    my $clone = $db->get_Clone($clone_id);
    foreach my $contig ($clone->get_all_Contigs){

	foreach my $gene ($contig->get_all_Genes){
	    if ($gene_count <= $num_genes){

		print "\tGene " . $gene->stable_id . "\n";

		# Retrieve all Transcript objects for a given gene using the each_Transcript method.
		my @transcripts = $gene->each_Transcript;

		if (scalar @transcripts > 1) {$multi_transcripts++;}

		foreach my $transcript (@transcripts){
		    print "\t\tTranscript " . $transcript->stable_id . "\n";

		    if ($gene_count <= $first_few_genes) {	
			# Using the translate method of a Transcript object will return a bioperl
			# Seq object from which the translated sequence can be derived.
			my $peptide = $transcript->translate;
			print "\t\t\tTranslation: " . $peptide->seq . ".\n";
		    }

		    # From a Transcript object it is possible to retrieve a corresponding
		    # array of Exon objects using the get_all_Exons method.
		    my @exons = $transcript->get_all_Exons;

		    $exon_counter += scalar @exons;
		    print "\t\t\tThis transcript has " . scalar @exons . " exons.\n";

		    foreach my $exon (@exons){
			my $exon_size = $exon->end - $exon->start;

			print "\t\t\t\tExon size is " . $exon_size . " bp\n";

			$gene_length_counter += $exon_size;
		    }
		}
		$gene_count++;

	    }else {last CLONE;}
	}
    }
} 

# Print a few statistics...

print "Out of " . ($gene_count-1) . " genes " . $multi_transcripts . " had more that one transcript.\n";

my $average_exons = $exon_counter / ($gene_count - 1);
my $average_length = $gene_length_counter / ($gene_count - 1);

print "On average, each gene was " . $average_length . " bp in length and had " . $average_exons . " exons.\n";


