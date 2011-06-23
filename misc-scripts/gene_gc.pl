# Calculate per-gene GC content and store as gene attributes

use strict;
use DBI;

use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my ( $host, $user, $pass, $port, $dbpattern, $print);

$port = 3306;

GetOptions( "host|h=s",       \$host,
	    "user|u=s",       \$user,
	    "pass|p=s",       \$pass,
	    "port=i",         \$port,
	    "pattern=s",      \$dbpattern,
	    "print",          \$print,
	    "help" ,          \&usage
	  );


usage() if (!$host || !$dbpattern || !$user || !$pass);

# loop over databases
my $dsn = "DBI:mysql:host=$host";
$dsn .= ";port=$port" if ($port);

my $db = DBI->connect($dsn, $user, $pass);

my @dbnames = map {$_->[0] } @{$db->selectall_arrayref("show databases")};

for my $dbname (@dbnames) {

  next if ($dbname !~ /$dbpattern/);

  my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor('-host' => $host,
					       '-port' => $port,
					       '-user' => $user,
					       '-pass' => $pass,
					       '-dbname' => $dbname,
					       '-species' => $dbname);

  print STDOUT "$dbname\n";

  delete_existing($dba) if !($print);

  print STDOUT "Calculating Gene GC attributes\n";

  my $attribute_adaptor = $dba->get_AttributeAdaptor();

  my $genes = $dba->get_GeneAdaptor()->fetch_all();

  my $total_count = 0;

  while (my $gene = shift(@$genes)) {

    my $gc = $gene->feature_Slice()->get_base_count->{'%gc'};

    if (!$print) {
      # attribute types need to be propagated from production database to all dbs
      # if the type exists it won't be overwritten
      my $attribute = Bio::EnsEMBL::Attribute->new(-CODE        => 'GeneGC',
						   -NAME        => 'Gene GC',
						   -DESCRIPTION => 'Percentage GC content for this gene',
						   -VALUE       => $gc);
      my @attributes = ($attribute);
      $attribute_adaptor->store_on_Gene($gene->dbID, \@attributes);
 
      $total_count++; 

    } else {

      print $gene->stable_id() . " " . $gc . "\n";

    }

  }
  if (!$print) {
      print STDOUT "Written $total_count 'GeneGC' gene attributes to database $dbname on server $host.\n";
  }

}

# ----------------------------------------------------------------------

sub delete_existing {

  my $dba = shift;

  print STDOUT "Deleting existing 'GeneGC' gene attributes\n";
  my $dsth = $dba->dbc()->prepare("DELETE ga FROM gene_attrib ga, attrib_type at WHERE at.attrib_type_id=ga.attrib_type_id AND at.code='GeneGC'");
  $dsth->execute();

}


sub usage {
  my $indent = ' ' x length($0);
  print <<EOF; exit(0);

What does it do?

This script calculates per-gene GC content and stores it as gene attributes.
It deletes existing Gene GC attributes. Then fetches all genes in the
core db and calculates the %gc for each gene.

Input data: dna sequence 
Output tables: gene_attrib 


When to run it in the release cycle?

It can be run whenever the genes and sequence are stable, i.e. any time after 
the genebuild handover, but before the handover to Mart.


Which databases to run it on?

It needs to be run across all core databases for every release.


How long does it take?

It takes a total of about 10 hours to run for all core databases in normal queue,


Usage:

  $0 -h host [-port port] -u user -p password \\
  $indent -pattern pattern [-print] \\
  $indent [-help]  \\

  -h|host              Database host to connect to

  -port                Database port to connect to (default 3306)

  -u|user              Database username

  -p|pass              Password for user

  -pattern             Database name regexp

  -print               Just print, don't insert or delete attributes

  -help                This message


EOF

}
