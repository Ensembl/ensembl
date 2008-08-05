# Calculate per-gene GC content and store as gene attributes

use strict;
use DBI;

use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my ( $host, $user, $pass, $port, $dbpattern, $print);

$port = 3306;

GetOptions( "dbhost|host=s",       \$host,
	    "dbuser|user=s",       \$user,
	    "dbpass|pass=s",       \$pass,
	    "dbport|port=i",       \$port,
	    "dbpattern|pattern=s", \$dbpattern,
	    "print",               \$print,
	    "help" ,               \&usage
	  );


usage() if (!$host || !$dbpattern);

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

  print STDERR "$dbname\n";

  delete_existing($dba) if !($print);

  print STDERR "Calculating Gene GC attributes\n";

  my $attribute_adaptor = $dba->get_AttributeAdaptor();

  my $genes = $dba->get_GeneAdaptor()->fetch_all();

  while (my $gene = shift(@$genes)) {

    my $gc = $gene->feature_Slice()->get_base_count->{'%gc'};

    if (!$print) {

      my $attribute = Bio::EnsEMBL::Attribute->new(-CODE        => 'GeneGC',
						   -NAME        => 'Gene GC',
						   -DESCRIPTION => 'Percentage GC content for this gene',
						   -VALUE       => $gc);
      my @attributes = ($attribute);
      $attribute_adaptor->store_on_Gene($gene->dbID, \@attributes);

    } else {

      print $gene->stable_id() . " " . $gc . "\n";

    }

  }
}

# ----------------------------------------------------------------------

sub delete_existing {

  my $dba = shift;

  print STDERR "Deleting existing 'GeneGC' gene attributes\n";
  my $dsth = $dba->dbc()->prepare("DELETE ga FROM gene_attrib ga, attrib_type at WHERE at.attrib_type_id=ga.attrib_type_id AND at.code='GeneGC'");
  $dsth->execute();

}


sub usage {
  print <<EOF; exit(0);

Calculate per-gene GC content and store as gene attributes.

Usage: perl $0 <options>

  -host|dbhost       Database host to connect to

  -port|dbport       Database port to connect to (default 3306)

  -dbpattern         Database name regexp

  -user|dbuser       Database username

  -pass|dbpass       Password for user

  -print             Just print, don't insert or delete attributes

  -help              This message


EOF

}
