# file - Ex4_Q789.pl
# author - Dan.Andrews@sanger.ac.uk
# desc - answers questions 7, 8 and 9 from Exercise 4 of the EnsEMBL api tutorial.

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::Seq;
use Bio::SeqIO;
use strict;

# Set a filename to write introns to futher on - choose something for your system.
my $work_file = './introns.fasta';

# Grab a database adaptor object - set your own system settings here or use the default.

my $host = 'kaka.sanger.ac.uk';
my $user = 'anonymous';
my $dbname = 'current';

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-host => $host,
					     -user => $user,
					     -dbname => $dbname);

# Open an bioperl IO stream for later writing to file.
my $seqio = Bio::SeqIO->new(-format => 'fasta',
			    -file => ">$work_file"
#			    -fh => \*STDOUT 
			    );

# Set the identification of the golden path being used, and then get an adaptor to this
# static golden path database.
$db->static_golden_path_type('NCBI_28');
my $sgp = $db->get_StaticGoldenPathAdaptor;

# Get a virtual contig to work with.
my $vcontig = $sgp->fetch_VirtualContig_by_chr_start_end('1',1,1000000);

# Start looping through all the genes on the virtual contig, counting to 10.
my $gene_count = 0;

foreach my $gene ($vcontig->get_all_Genes){
    if ($gene_count <= 10) {
	$gene_count++;
	print "Gene " . $gene->stable_id . "\n";

	# Make sure to look at each transcript for each gene, not just the first.
	foreach my $transcript ($gene->each_Transcript){
	    my @starts = <undef>;
	    my @ends = <undef>;

	    print "     Transcript ". $transcript->stable_id ."\n";

	    # Derive the sequence objects representing the UTR sequences, if defined, and from
	    # these print their text sequences.
	    my $five_prime_utr = $transcript->five_prime_utr;
	    my $three_prime_utr = $transcript->three_prime_utr;
	    if (defined $five_prime_utr){
		print "          5\' UTR\: " . $five_prime_utr->seq . "\n";
	    } else {print "          No 5\' UTR was found for this transcript\n";}
	    if (defined $three_prime_utr){
		print "          3\' UTR\: " . $three_prime_utr->seq . "\n";
	    } else {print "          No 3\' UTR was found for this transcript\n";}

	    # There is no guarantee that the order of the exons as they are returned from the database
	    # will be sequential along the gene - so run a sort on the transcript object to fix this.
	    $transcript->sort;

	    # Work through each exon as a way to determining the INTRON start and end locations.
	    foreach my $exon ($transcript->get_all_Exons){
		print "          Exon " . $exon->stable_id . "\t" . $exon->strand . "\t" . $exon->start . "\t" . $exon->end . "\n";
		print "               200bp upstream\:   " . $vcontig->subseq((($exon->start)-200),$exon->start) . "\n";
		print "               200bp downstream\: " . $vcontig->subseq($exon->end, (($exon->end)+200)) . "\n";
		push @starts, ($exon->start);
		push @ends, ($exon->end);
	    } 
	    if (defined $starts[1]){              # Check whether the gene is on the reverse strand and if so reverse the coordinates.
		if ($starts[0] > $starts[1]){
		    @starts = reverse @starts;
		    @ends = reverse @ends;
		}
	    }
	    my @intron_starts = @ends;            # Final juggling exon start/stop locations to derive the intron coordinates.
	    pop @intron_starts;
	    my @intron_ends = @starts;
	    shift @intron_ends;

	    print "          Dumping introns for gene " . $gene->stable_id . " transcript " . 
		$transcript->stable_id . " to file " . $work_file . "\n\n";

	    # Writing each intron sequence to the seq IO output stream.
	    for(my $i = 0; $i <= $#intron_starts; $i++){
		my $intron_seq = Bio::Seq->new(-seq => $vcontig->subseq($intron_starts[$i], $intron_ends[$i]), 
					       -display_id => $gene->stable_id . " " . $transcript->stable_id . " intron " . ($i + 1));
		$seqio->write_seq($intron_seq);
	    }
	}
    } else {last;}  # After 10 genes we've had enough.
}






