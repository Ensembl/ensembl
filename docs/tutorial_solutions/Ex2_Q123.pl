# file - Ex2_Q123.pl
# author - Dan.Andrews@sanger.ac.uk
# desc - answers questions 1, 2 and 3 from Exercise 2 of the EnsEMBL api tutorial.

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use strict;

# Use a DBAdaptor object to connect to an EnsEMBL database.
my $host   = 'kaka.sanger.ac.uk';
my $user   = 'anonymous';
my $dbname = 'current';

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $host,
					    -user   => $user,
					    -dbname => $dbname);

# Grab an array of all clone ids
my @clone_list = $db->get_all_Clone_id;

# Iterate through the clone ids extracting the contigs comprising each clone.
my $contig_count = 1;
my $repeat_length_tally = 0;
my $contig_length_tally = 0;

CLONE:
foreach $clone_acc (@clone_list){
    my $clone = $db->get_Clone($clone_acc);

    foreach my $contig ($clone->get_all_Contigs){
	 if ($contig_count <= 100){

	     $contig_length_tally += $contig->length;           # Tally the lengths of the contigs.

             # Ask for the repeat features on the contig and work through them printing details.
	     my @repeat_feat = $contig->get_all_RepeatFeatures;
	     foreach my $repeat_feat (@repeat_feat){
		 $repeat_length_tally += $repeat_feat->length;  # Tally the lengths of the repeats
		 print $repeat_feat->gffstring . "\n";
	     }

	     $contig_count++;

	 } else {last CLONE;}
    }
}

print "Repetitive sequences comprise " . ($repeat_length_tally/$contig_length_tally)*100 . 
    "\% of the total sequence over the first 100 contigs.\n";


# Now working with clone AC005663

my %hits;
my $all_hits = 0;
my $sig_hits = 0;
my $arbitrary_significance_threshold = 1e-50;  ## Choose your own significance cut-off.

my $clone = $db->get_Clone('AC005663');

# Work through the contigs on the clone calling for all the similarity features on each contig.

foreach my $contig ($clone->get_all_Contigs){

    # Look at each similarity feature - count it and assess its significance.
    foreach my $sim_feature ($contig->get_all_SimilarityFeatures){
	if ($sim_feature->isa('Bio::EnsEMBL::FeaturePairI')){   ## Check the returned object *is* a FeaturePair object.
                                                                ## Occasionally a SeqFeature object is returned instead
                                                                ## (e.g. for cpg islands) and this object will neither 
                                                                ## have a hseqname nor a p_value method.
	    $all_hits++;
	    $hits{$sim_feature->hseqname}++;
	    if (($sim_feature->p_value < $arbitrary_significance_threshold)&  # <- Check the p-value here
		(defined $sim_feature->p_value)&                              # <- and here ignore messy (empty or NULL) entries
		($sim_feature->p_value ne 'NULL')){                           #  in the db.
		$sig_hits++;
	    }
	}
    }
}
    
print "In total, there were " . $all_hits . 
    " hits found for this clone.  These hits represented matches to " . 
    scalar (keys %hits) . " distinct sequences.  " . 
    ($sig_hits/$all_hits)*100 . " of all matches were significant at a threshold p-value of " . 
    $arbitrary_significance_threshold . ".\n";
















