# file - Ex6_Q12345.pl
# author - Dan.Andrews@sanger.ac.uk
# desc - answers questions 1, 2, 3, 4 and 5 from Exercise 6 of the EnsEMBL api tutorial.

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::ExternalData::EXONERATESQL::DBAdaptor;
use strict;


# Use a DBAdaptor object to connect to an EnsEMBL database.
my $host   = 'kaka.sanger.ac.uk';
my $user   = 'anonymous';
my $dbname = 'current';

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-host   => $host,
					     -user   => $user,
					     -dbname => $dbname);


# Also connect to the mouse database using a specialised db adaptor.
my $mouse_dbname = 'human_mouse_130';

my $mouse_db = Bio::EnsEMBL::ExternalData::EXONERATESQL::DBAdaptor->new(-host   => $host,
									-user   => $user,
									-dbname => $mouse_dbname);

# Get an adaptor object to create an 'external feature factory'
my $mouse_external_feature_factory = $mouse_db->get_ExonerateAdaptor;

# Add the external data adaptor to the main db adaptor object.
$db->add_ExternalFeatureFactory($mouse_external_feature_factory);

# Set the golden path assembly, create a static golden path adaptor and
# retrieve a virtual contig.
$db->static_golden_path_type('NCBI_28');

my $sgp = $db->get_StaticGoldenPathAdaptor;

my $vcontig = $sgp->fetch_VirtualContig_by_chr_start_end('1', 1, 1000000);

# Now that the external database is plugged in, get the features from this
# db that lie on our virtual contig.
my @mouse_features = $vcontig->get_all_ExternalFeatures;

print "There are " . scalar @mouse_features . " hits of mouse features to this portion of sequence.\n";

# Create an array of all exons to work with.
my @all_exons;
foreach my $gene ($vcontig->get_all_Genes){
    foreach my $exon ($gene->get_all_Exons){
	push (@all_exons, $exon);
    }
}

# Work through each exon and see if mouse features overlap with human exons.
my $exon2mousefeat_hit = 0;
foreach my $ext_feature ($vcontig->get_all_ExternalFeatures){
#    print "Mouse hit : " . $ext_feature->gffstring . "\n";
    foreach my $exon (@all_exons){
	if ($exon->overlaps($ext_feature)){
	    $exon2mousefeat_hit++;
	}
    }
}

# Print our findings.

print $exon2mousefeat_hit . " mouse hits overlap human exons.\n";
print "Judging from this data, this means that a mouse homolog possibly exists for " . ($exon2mousefeat_hit/(scalar @all_exons))*100 . "% of human exons.\n";

