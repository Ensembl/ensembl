# file - Ex4_Q56.pl
# author - Dan.Andrews@sanger.ac.uk
# desc - answers questions 5 and 6 from Exercise 4 of the EnsEMBL api tutorial.

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

# Grab a virtual contig of gene ENSG00000100259 and take 5000bp of the upstream and downstream sequence.
my $vcontig = $sgp->fetch_VirtualContig_of_gene('ENSG00000100259', 5000);

# Grab all the genes on this virtual contig.
my @genes = $vcontig->get_all_Genes;

# Check through the list of returned genes for the actual gene we want
# (there is only one gene, but a little paranoia is not altogether unwarranted).

foreach my $gene (@genes){
    if ($gene->stable_id == 'ENSG00000100259'){
	print "Found gene ENSG00000100259.\n";
	my $exon_counter = 0;

	# Work through a list of the exon objects for this gene.  Using the start and end
	# coordinates of the exons, print the flanking sequences.  We *know* that this gene
	# is on the positive strand and hence we dont need to faff about with exon direction.
	foreach my $exon ($gene->get_all_Exons){
	    $exon_counter++;
	    print "Exon number  " . $exon_counter . "\nExon start " . $exon->start . "\t\tEnd " . $exon->end . "\n";
	    my $up_start = $exon->start - 200;
	    my $down_end = $exon->end + 200;
	    print "Upstream 200 bp - " . $vcontig->subseq($up_start,$exon->start) . "\n";
	    print "Downstream 200 bp - " . $vcontig->subseq($exon->end,$down_end) . "\n";
	}
	last;
    }
}

