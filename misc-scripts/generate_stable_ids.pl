# Generate stable IDs for genes/transcripts/translations/exons that have none
# Start from current max stable ID + 1

use strict;
use DBI;
use Getopt::Long;

my ($host, $port, $dbname, $user, $pass, @types, $start);

GetOptions('user=s'       => \$user,
	   'pass=s'       => \$pass,
	   'host=s'       => \$host,
	   'port=i'       => \$port,
	   'dbname=s'     => \$dbname,
	   'types=s'      => \@types,
	   'start=s'      => \$start,
	   'help'         => sub { usage(); exit(0); });

@types = ('gene','transcript','translation','exon') if (!@types);
@types = split(/,/,join(',',@types));

if (!$user || !$host || !$dbname || !@types) {

  usage();
  exit(1);

}

my $dbi = DBI->connect( "DBI:mysql:host=$host:port=$port;database=$dbname", $user, $pass,
			{'RaiseError' => 1}) || die "Can't connect to database\n";

foreach my $type (@types) {

  my $table = $type . "_stable_id";
  my $sth;

  # get starting stable ID, either specified or current max
  my $new_stable_id = $start ? $start : get_highest_stable_id($dbi, $type);

  # get timestamp so all new stable IDs have the same created/modified dates
  $sth = $dbi->prepare("SELECT NOW()");
  $sth->execute();
  my $ts;
  if (my @row = $sth->fetchrow_array()) {
    $ts= $row[0];
  } else {
    die "Can't get timestamp\n";
  }

  # get a list of objects that don't currently have stable IDs assigned
  # and assign new ones, incrementing & preserving formatting as we go
  my $sql = "SELECT $type.${type}_id FROM $type LEFT JOIN $table sid ON $type.${type}_id=sid.${type}_id WHERE sid.stable_id IS NULL";
  $sth = $dbi->prepare($sql);
  $sth->execute();

  while (my @row = $sth->fetchrow_array()) {
    $new_stable_id = increment_stable_id($new_stable_id);
    print "INSERT INTO $table VALUES($row[0],\'$new_stable_id\',1,\'$ts\',\'$ts\');\n";
  }

}

# --------------------------------------------------------------------------------

sub increment_stable_id {

  my $stable_id = shift;

  my ($prefix,$suffix) = $stable_id =~ /([a-zA-Z]+)([0-9]+)/;

  return sprintf "%s%011d", $prefix, $suffix+1;

}

# --------------------------------------------------------------------------------

sub get_highest_stable_id {

  my ($dbi, $type) = @_;

  my $sid = $type . "_stable_id";

  my ($highest_from_current, $highest_from_archive);

  # get highest stable ID from the relevant table
  my $sth = $dbi->prepare("SELECT MAX(stable_id) FROM $sid");
  $sth->execute();
  if (my @row = $sth->fetchrow_array()) {
    $highest_from_current = $row[0];
  } else {
    die "Can't get max $type stable ID from $sid\n";
  }

  return $highest_from_current if ($type eq "exon");

  # and from relevant archive
  $sth = $dbi->prepare("SELECT MAX($sid) FROM gene_archive");
  $sth->execute();
  if (my @row = $sth->fetchrow_array()) {
    $highest_from_archive = $row[0];
  } else {
    die "Can't get max $type stable ID from gene_archive\n";
  }

  my $max = ($highest_from_current ge $highest_from_archive) ? $highest_from_current : $highest_from_archive;

  return $max;

}
# --------------------------------------------------------------------------------

sub usage {

  print << "EOF";

  generate_stable_ids.pl -user {user} -pass {password} -host {host} -port {port} -dbname {database} -types {gene,exon,transcript,translation} -start {first stable ID}

  Argument to -types is a comma-separated list of types of stable IDs to be produced.

  If the -types argument is ommitted, stable IDs are generated for all types (gene,transcript,translation,exon).

  Assigns stable IDs to objects that currently have none. Starting stable ID is found by incrementing the highest current stable ID for that type *or* by using -start argument.

  Note if -start is used you must generate gene, transcript, exon and translation stable IDs separately in order to specify the correct prefix for each.

  The parameter to -start should be the stable ID you wish to start from, e.g. ENSG00000000001 as this is used to obtain the prefix (e.g. ENSG).

  Produces SQL which can be run against the target database.

EOF

}

# --------------------------------------------------------------------------------
