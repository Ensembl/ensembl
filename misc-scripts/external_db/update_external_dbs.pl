#!/usr/local/ensembl/bin/perl -w

#
# updates the external db tables on all of the core databases on a given host
#

use strict;

use Getopt::Long;
use DBI;
use IO::File;

my ( $host, $user,        $pass,   $port, @dbnames,
     $file, $release_num, $master, $nonreleasemode, $force );

GetOptions( "dbhost|host=s",        \$host,
            "dbuser|user=s",        \$user,
            "dbpass|pass=s",        \$pass,
            "dbport|port=i",        \$port,
            "file=s",        \$file,
            "dbnames=s@",    \@dbnames,
            "release_num=i", \$release_num,
            "master=s",      \$master,
	    "nonreleasemode", \$nonreleasemode, 
            "force",         \$force );

$port ||= 3306;

$file ||= "external_dbs.txt";

usage("[DIE] Need a host") if(!$host);

#release num XOR dbname are required.
usage(   "[DIE] Need either both a release number and "
       . "database names or neither" )
  if ( ( $release_num && @dbnames ) || ( !$release_num && !@dbnames ) );

if(!$nonreleasemode){ 
  # master database is required
  usage("[DIE] Master database required") if (!$master);
}

my $dsn = "DBI:mysql:host=$host;port=$port";

my $db = DBI->connect( $dsn, $user, $pass, {RaiseError => 1} );

if($release_num) {
  @dbnames = map {$_->[0] } @{ $db->selectall_arrayref( "show databases" ) };

  #
  # filter out all non-core databases
  #
  @dbnames = grep {/^[a-zA-Z]+\_[a-zA-Z]+\_(core|est|estgene|vega|otherfeatures|cdna)\_${release_num}\_\d+[A-Za-z]?$/} @dbnames;
}

my @field_names = qw(external_db_id db_name release status dbprimary_acc_linkable display_label_linkable priority  db_display_name type);

my @types = qw(ARRAY ALT_TRANS MISC LIT PRIMARY_DB_SYNONYM ALT_GENE);

#
# make sure the user wishes to continue
#
print STDERR "Please make sure you've updated $file from CVS!\n";
print STDERR
  "The following databases will have their external_db tables "
  . "updated if necessary:\n  ";
print join( "\n  ", @dbnames );
print "\nContinue with update (yes/no)>  ";

my $input = lc(<STDIN>);
chomp($input);
if ($input ne 'yes') {
  print "external_db conversion aborted\n";
  exit();
}

#
# read all of the new external_db entries from the file
#
my $fh = IO::File->new();
$fh->open($file) or die("Could not open input file $file");
my @rows;
my %bad_lines;

while (my $row = <$fh>) {
  chomp($row);
  next if ($row =~ /^#/);    # skip comments
  next if ($row =~ /^$/);    # and blank lines
  next if ($row =~ /^\s+$/); # and whitespace-only lines
  my @a = split(/\t/, $row);
  push @rows, {
      'external_db_id'         => $a[0],
      'db_name'                => $a[1],
      'release'                => $a[2],
      'status'                 => $a[3],
      'dbprimary_acc_linkable' => $a[4],
      'display_label_linkable' => $a[5],
      'priority'               => $a[6],
      'db_display_name'        => $a[7],
      'type'                   => $a[8] };

  if ( $a[1] =~ /-/ ) {
    print STDERR "Database name "
      . $a[1]
      . " contains '-' characters "
      . "which will break Mart, "
      . "please replace them with '_' until Mart is fixed\n";
    exit(1);
  }

  # do some formatting checks
  my $blank;
  for (my $i=0; $i < scalar(@a); $i++) {
    if ($a[$i] eq '') {
      $bad_lines{$row} = $field_names[$i] . " - field blank - check all tabs/spaces in line";
    }
  }

  if ($a[1] =~ /\s/) {
    $bad_lines{$row} = "db_name field appears to contain spaces";
  }
  if ($a[1] =~ /^$/) {
    $bad_lines{$row} = "db_name field appears to be missing";
  }
  if ($a[1] =~ /^\s+$/) {
    $bad_lines{$row} = "db_name field appears to be blank";
  }
  if ($a[1] =~ /^\d+$/) {
    $bad_lines{$row} = "db_name field appears to be numeric - check formatting";
  }

  my $type_ok;
  foreach my $type (@types) {
    $type_ok = 1 if ($a[8] eq $type);
  }
  $bad_lines{$row} = "type field is " . $a[8] . ", not one of the recognised types" if (!$type_ok);

}
$fh->close();

if (%bad_lines) {
  print STDERR "Cannot parse the following line(s) from $file; check that all fields are present and are separated by one tab (not spaces). \n";
  print STDERR "Name of problem field, and the error is printed in brackets first\n\n";
  foreach my $row (keys %bad_lines) {
    print STDERR "[". $bad_lines{$row} . "]" .  " $row\n";
  }
  exit(1);
}

# Load into master database
  if(!$nonreleasemode){
    load_database($db, $master, @rows);
  }
  # Check each other database in turn
  # Load if no extra rows in db that aren't in master
  # Warn and skip if there are

  foreach my $dbname (@dbnames) {

    print STDERR "Looking at $dbname ... \n";
    if ($force || $nonreleasemode) {

      print STDERR "Forcing overwrite of external_db table in "
	. "$dbname from $file\n";
      load_database( $db, $dbname, @rows );

    } elsif (compare_external_db($db, $master, $dbname)) {

      print STDERR "$dbname has no additional rows. "
	. "Overwriting external_db table from $file\n";
      load_database( $db, $dbname, @rows );

    } else {

      print STDERR "$dbname has extra rows "
	. "that are not in $file, skipping\n";
	
    }

  }

print STDERR "Updates complete\n";



sub load_database {
  my ($db, $dbname, @rows) = @_;

  $db->do("USE $dbname");

  # Save all existing release information from the table.
  my $sth = $db->prepare(
    qq( SELECT  external_db_id, db_release
        FROM    external_db) );
  $sth->execute();

  my %saved_release;
  while ( my ( $id, $release ) = $sth->fetchrow_array() ) {
    if ( defined($release) && $release ne '1' ) {
      $saved_release{$id} = $release;
    }
  }
  $sth->finish();

  # Delete the existing table
  $sth = $db->prepare('DELETE FROM external_db');
  $sth->execute();
  $sth->finish();

  # Populate the table with data from the file (using the saved release
  # information)
  $sth = $db->prepare(
    qq( INSERT INTO external_db (
          external_db_id,
          db_name,
          db_release,
          status,
          dbprimary_acc_linkable,
          display_label_linkable,
          priority,
          db_display_name,
          type) VALUES (?,?,?,?,?,?,?,?,?)) );

  foreach my $row (@rows) {
    my $id = $row->{'external_db_id'};

    $sth->execute( $id,
                   $row->{'db_name'}, (
                     exists( $saved_release{$id} )
                     ? $saved_release{$id}
                     : $row->{'release'}
                   ),
                   $row->{'status'},
                   $row->{'dbprimary_acc_linkable'},
                   $row->{'display_label_linkable'},
                   $row->{'priority'},
                   $row->{'db_display_name'},
                   $row->{'type'} );
  }

  $sth->finish();
}


# return true if the tables are the same, undef if not
sub compare_external_db {

  my ($db, $master, $dbname) = @_;

  my $same = 1;

  # check each row in $dbname against each row in $master
  # only compare ID and name since we're only aiming to catch extra rows in $dbname
  $db->do("use $dbname");

  my $sth = $db->prepare(qq {SELECT d.external_db_id, d.db_name
			     FROM $dbname.external_db d
			     LEFT JOIN $master.external_db m
			     ON (d.external_db_id=m.external_db_id AND d.db_name=m.db_name)
			     WHERE m.external_db_id IS NULL OR m.db_name IS NULL });
  $sth->execute();

  while (my ($id, $external_db_name) = $sth->fetchrow_array) {
    print "$dbname has external_db entry for $external_db_name (ID $id) which is not present in $master\n";
    $same = undef;

  }

  $sth->finish();

  return $same;

}

sub usage {
  my $error = shift;

  print STDERR <<EOF;
  $error

  Usage: $0 options

    -host       hostname
    -user       username
    -pass       password
    -port       port_of_server optional
    -master     the name of the master database to load the file into
    -force      force update, even if there are rows in the database
                that are not in the file
    -release    the release of the database to update used to match
                database names, e.g. 13
    -file       the path of the file containing the insert statements
                of the entries of the external_db table.  Default is
                'external_dbs.txt'
    -dbnames    the names of the database to update.  If not provided
                all of the core databases matching the release arg
                will be updated.  Either -dbnames or -release must
                be specified, but not both.  Multiple dbnames can be
                provided.
    -nonreleasemode Does not require master schema and forces the update.


  Examples:

  # Update two databases

  ./update_external_dbs.pl -host ecs1c -file external_dbs.txt \\
    -user ensadmin -pass secret -master master_schema_14 \\
    -dbnames homo_sapiens_core_14_33 -dbnames mus_musculus_core_14_30

  # Update all Core databases for release 14

  ./update_external_dbs.pl -host ens-staging -file external_dbs.txt \\
    -user ensadmin -pass secret -release 42 -master master_schema_42

  If the databases to be updated contain rows that are not in the file,
  a warning will be given and the database in question skipped, unless
  -force is used.

  This program will not overwrite the db_release column of any table.

EOF
    exit;
} ## end sub usage
