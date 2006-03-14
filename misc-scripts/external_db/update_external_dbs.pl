#
# updates the external db tables on all of the core databases on a given host
#


use strict;

use Getopt::Long;
use DBI;
use IO::File;

my ( $host, $user, $pass, $port,@dbnames, $file, $release_num);

GetOptions( "host=s", \$host,
	    "user=s", \$user,
	    "pass=s", \$pass,
	    "port=i", \$port,
	    "file=s", \$file,
            "dbnames=s@", \@dbnames,
	    "release_num=i", \$release_num
	  );

#both host and file are required
usage() if(!$host || !$file);
#release num XOR dbname are required
usage() if(($release_num && @dbnames) || (!$release_num && !@dbnames));

$port ||= 3306;

my $dsn = "DBI:mysql:host=$host;port=$port";

my $db = DBI->connect( $dsn, $user, $pass, {RaiseError => 1} );

if($release_num) {
  @dbnames = map {$_->[0] } @{ $db->selectall_arrayref( "show databases" ) };
  
  #
  # filter out all non-core databases
  #
  @dbnames = grep {/^[a-zA-Z]+\_[a-zA-Z]+\_(core|est|estgene|vega)\_${release_num}\_\d+[A-Za-z]?$/} @dbnames;
}


#
# make sure the user wishes to continue
#
print STDERR "The following databases will be external_db updated:\n  ";
print join("\n  ", @dbnames);
print "\ncontinue with update (yes/no)>  ";

my $input = lc(<STDIN>);
chomp($input);
if($input ne 'yes') {
  print "external_db conversion aborted\n";
  exit();
}

#
# read all of the new external_db entries from the file
#
my $fh = IO::File->new();
$fh->open($file) or die("could not open input file $file");
my @rows;
my $row;
while($row = <$fh>) {
  chomp($row);
  my @a = split(/\t/, $row);
  push @rows, {'external_db_id'         => $a[0],
               'db_name'                => $a[1],
               'release'                => $a[2],
               'status'                 => $a[3],
	       'dbprimary_acc_linkable' => $a[4],
	       'display_label_linkable' => $a[5],
	       'priority'               => $a[6],
	       'db_display_name'        => $a[7]};
}
$fh->close();

foreach my $dbname (@dbnames) {
  print STDERR "updating $dbname\n";
  $db->do("use $dbname");
  my $sth = $db->prepare('DELETE FROM external_db');
  $sth->execute();
  $sth->finish();

  $sth = $db->prepare('INSERT INTO external_db (external_db_id, db_name,
                                                db_release, status, dbprimary_acc_linkable,
                                                display_label_linkable, priority,
                                                db_display_name)
                       VALUES (?,?,?,?,?,?,?,?)');

  foreach my $row (@rows) {
    $sth->execute($row->{'external_db_id'},
		  $row->{'db_name'},
		  $row->{'release'},
		  $row->{'status'},
		  $row->{'dbprimary_acc_linkable'},
		  $row->{'display_label_linkable'},
		  $row->{'priority'},
		  $row->{'db_display_name'});
  }

  $sth->finish();
}

print STDERR "updates complete\n";


sub usage {
  print STDERR <<EOF

             Usage: update_external_db options
 Where options are: -host hostname 
                    -user username 
                    -pass password 
                    -port port_of_server optional
                    -release the release of the database to update used to 
                             match database names.  e.g. 13
                    -file the path of the file containing the insert statements
                          of the entries of the external_db table
                    -dbnames db1
                          the names of the database to update. if not provided
                          all of the core databases matching the release arg
                          will be updated.  Either -dbnames or -release must
                          be specified, but not both.  Multiple dbnames can
                          be provided.

 E.g.:

  #update 2 databases
  perl update_external_dbs.pl -host ecs1c -file external_dbs.txt -user ensadmin -pass secret -dbnames homo_sapiens_core_14_33 -dbnames mus_musculus_core_14_30

  #update all core databases for release 14
  perl update_external_dbs.pl -host ecs2d -file external_dbs.txt -user ensadmin -pass secret -release 14

EOF
;
  exit;
}
