# Figure out which seq_regions are non-redundant and set the appropriate
# attribute in seq_region_attrib

use strict;
use warnings;

use DBI;
use Getopt::Long;

my ($host, $port, $user, $password, $db, $verbose);
$host = "127.0.0.1";
$port = 5000;
$password = "";
$user = "ensro";

GetOptions ('host=s'      => \$host,
            'user=s'      => \$user,
            'password=s'  => \$password,
            'port=s'      => \$port,
            'db=s'        => \$db,
            'verbose'     => \$verbose,
            'help'        => sub { &show_help(); exit 1;} );

die "Host must be specified"           unless $host;
die "Database must be specified"       unless $db;

my $dbi = DBI->connect("dbi:mysql:host=$host;port=$port;database=$db", "$user", "$password", {'RaiseError' => 1})  || die "Can't connect to target DB";

# check that there is an entry in the attrib_type table for nonredundant
# if there is, cache the ID for later; if not, make one
my $nr_attrib_type_id;
my $sth = $dbi->prepare("SELECT attrib_type_id FROM attrib_type WHERE code='nonredundant'");
$sth->execute();
while ( my @row= $sth->fetchrow_array()) {
  $nr_attrib_type_id = $row[0];
}
$sth->finish();

if ($nr_attrib_type_id) {

  debug("Attribute with name nonredundant already set in attrib_type with attrib_type_id " . $nr_attrib_type_id);

} else {

  debug("Attribute with name nonredundant not set in attrib_type; adding ...");

  $sth = $dbi->prepare("INSERT INTO attrib_type (code, name) VALUES ('nonredundant', 'Non-redundant sequence region')");
  $sth->execute();

  $nr_attrib_type_id = $sth->{mysql_insertid};

  debug("Added nonredundant attribute with ID " . $nr_attrib_type_id);

}


# ----------------------------------------------------------------------
# Misc / utility functions

sub show_help {

  print "Usage: perl set_non_redundant_attribs.pl {options}\n";
  print "Where options are:\n";
  print "  --host {hostname} The database host.\n";
  print "  --user {username} The database user. Must have write permissions\n";
  print "  --password {pass} The password for user, if required.\n";
  print "  --port {folder}   The database port to use.\n";
  print "  --db {schema}     The name of the database\n";
  print "  --verbose         Print extra output information\n";

}

# ----------------------------------------------------------------------

sub debug {

  my $str = shift;

  print $str . "\n" if $verbose;

}
