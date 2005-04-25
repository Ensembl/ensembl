# Check how many genes have regulatory features associated

use strict;

use DBI;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

my ($host, $user, $pass, $port, $dbname, $file, $del);

GetOptions( "host=s",   \$host,
	    "user=s",   \$user,
	    "pass=s",   \$pass,
	    "port=i",   \$port,
	    "dbname=s", \$dbname,
	    "help",     \&usage);

my $db_adaptor = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host,
						    -user => $user,
						    -port => $port,
						    -pass => $pass,
						    -dbname => $dbname);

my $gene_adaptor = $db_adaptor->get_GeneAdaptor();

my @genes = @{$gene_adaptor->list_dbIDs()};

my $total_genes = 0;
my $genes_with_rf = 0;

$| = 1; # flush after every print
print "Genes analysed: ";

foreach my $gene_id (@genes) {

  my $gene = $gene_adaptor->fetch_by_dbID($gene_id);
  $total_genes++;
  my @features = @{$gene->get_all_regulatory_features(1)};
  $genes_with_rf++ if (@features);
  print $total_genes . " " if ($total_genes % 1000 == 0);

}

my $pc = (100 * $genes_with_rf) / $total_genes;

printf "\n\nTotal of %s genes, of which %d (%2.1f%%) have associated regulatory features.\n", $total_genes, $genes_with_rf, $pc;
