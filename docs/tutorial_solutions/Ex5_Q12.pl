# file - Ex5_Q12.pl
# author - Dan.Andrews@sanger.ac.uk
# desc - answers questions 1 and 2 from Exercise 5 of the EnsEMBL api tutorial.

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use strict;

# Use a DBAdaptor object to connect to an EnsEMBL database.
my $host   = 'kaka.sanger.ac.uk';
my $user   = 'anonymous';
my $dbname = 'current';

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $host,
					    -user   => $user,
					    -dbname => $dbname);


# We need a protein adaptor such that we can later grab a protein object.
my $protein_adaptor = $db->get_Protein_Adaptor;

my $gene_count = 0;
my $genes_to_count = 50;
my %pfam_hits_hash;
my $toggle;
my @featureless;
my $wd40_count = 0;

# To find the first 100 genes we can perhaps use a list of clone ids that we 
# can progessively retrieve.
my @clone_list = $db->get_all_Clone_id;

# The next 12 or so lines are all fiddling around used to get a decent Protein object or two.
CLONE:
foreach my $clone_id (@clone_list){
    my $clone = $db->get_Clone($clone_id);
    print "Clone " . $clone->id . "\n";

    my @contigs = $clone->get_all_Contigs;
    foreach my $contig ($clone->get_all_Contigs){
	print "\tContig " . $contig->id . "\n";

	foreach my $gene ($contig->get_all_Genes){
	    $gene_count++;
	    print "\t\tGene " . $gene->stable_id . "\n"; 
	    if ($gene_count <= $genes_to_count){
		foreach my $transcript ($gene->each_Transcript){
		    print "\t\t\tTranscript " . $transcript->stable_id . "\n";

		    # Ok.  Finally we can try and get a Protein object.  There are two possible ways...
#		    my $translation = $transcript->translation;
#		    my $protein = $protein_adaptor->fetch_Protein_by_dbid($translation->dbID);

		    # Or, we can also use...
		    my $protein = $protein_adaptor->fetch_Protein_by_transcriptId($transcript->stable_id);

		    # Once we have a Protein object we can extract information about various features.
		    # You can get the features from the protein using...
		    # get_Family
		    # get_all_CoilsFeatures
		    # get_all_DBLinks
		    # get_all_DomainFeatures
		    # get_all_IntronFeatures
		    # get_all_LowcomplFeatures
		    # get_all_PfamFeatures
		    # get_all_PrintsFeatures
		    # get_all_PrositeFeatures
		    # get_all_SigpFeatures
		    # get_all_SnpsFeatures
		    # get_all_SuperfamilyFeatures
		    # get_all_TransmembraneFeatures
		    # among other methods...

		    print "\t\t\t\tPfam features\:\n";
		    my @pfam_features = $protein->get_all_PfamFeatures;  # Returns an array of Protein_FeaturePair objects
                                                                         # which are objects that inherit from FeaturePair.
                                                                         # These objects contain all sorts of useful
                                                                         # information (have a look).

		    foreach my $pfam_feat (@pfam_features){
			print "\t\t\t\t" . $pfam_feat->hseqname. "\n";   # <- One such useful piece of information

			# Store this useful information for later
			$pfam_hits_hash{$pfam_feat->hseqname}++;

			if (($pfam_feat->idesc) =~/WD40/i){
			    $wd40_count++;
			}
		    }

		    # Write something if there are no features.
		    if (scalar @pfam_features == 0){
			print "\t\t\t\tNone.\n";
		    }

		    my @all_feat = $protein->all_SeqFeature;
		    if (scalar @all_feat == 0){
			$toggle = 1;
		    } else {$toggle = 0;}
		}

		# Keep a record of any featureless genes/proteins.
		if ($toggle) {push @featureless, ($gene->stable_id);}
	    }else {last CLONE;}
	}
    }
}

# Print the information that has been gathered.

print "Pfam Domain\tOccurances\n";
foreach my $domain (keys %pfam_hits_hash){
    print $domain . "\t-\t" . $pfam_hits_hash{$domain} . "\n";
}

print "The following genes did not have features...\n";
foreach my $featureless_gene (@featureless){
    print $featureless_gene . "\n";
}
my $num_featureless = @featureless;
if ($num_featureless == 0) {print "None.\n";}

print "Saw " . $wd40_count . " instances of the WD40 domain in this sample of proteins.\n";


