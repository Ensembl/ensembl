# file - Ex1_Q123.pl
# author - Dan.Andrews@sanger.ac.uk
# desc - answers questions 1, 2 and 3 from Exercise 1 of the EnsEMBL api tutorial.

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use strict;

# Use a DBAdaptor object to connect to an EnsEMBL database.
my $host   = 'kaka.sanger.ac.uk';
my $user   = 'anonymous';
my $dbname = 'current';

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host =>$host,
					    -user =>$user,
					    -dbname => $dbname);

# Determine the number of clones present in the database.
my @clone_ids = $db->get_all_Clone_id;
my $num_clones = scalar @clone_ids;

print 'There are ' . $num_clones . ' clones in the ' . $dbname . ' database' . "\n";

# Loop through the first 100 clones and count the number of contigs.
my $total_contigs = 0;
my $count = 0;
my $cutoff = 100;

foreach my $clone_id (@clone_ids) {
  my $clone = $db->get_Clone($clone_id);     # Get a clone object using its id
  my @contigs = $clone->get_all_Contigs;     # Retrieve an array of all contig objects for a given clone.

  $total_contigs += scalar @contigs;

  $count++;
  last if $count>=$cutoff;
}

print 'On average, there are ' . ($total_contigs/$cutoff) . ' contigs per clone' . "\n";

# Open a bioperl SeqIO output stream.
my $seqio = Bio::SeqIO->new(-fh     => \*STDOUT,
			    -format => 'fasta');

# Loop through the to the last 10 clone ids and write out their repeatmasked sequences.
for(my $i = ($num_clones - 10); $i < $num_clones; $i++){ # $num_clones and @clone_ids are carried over from above. 
    my $clone = $db->get_Clone($clone_ids[$i]);
    my @contigs = $clone->get_all_Contigs;
    foreach my $contig (@contigs) {
	$seqio->write_seq($contig->get_repeatmasked_seq);
    }
}



































