#
# updates the unmapped_reason tables on all of the core databases on a given host
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
            "dbnames=s@", \@dbnames, # either provide -dbnames or -release  
	    "release_num=i", \$release_num
	  );

#both host and file are required
usage() if(!$host || !$file);
#release num XOR dbname are required
if(($release_num && @dbnames) || (!$release_num && !@dbnames)) {  
  print "\nYou can't use -dbnames <> and -release_num <> options at the same time\n" ;
  sleep(3) ;  
  usage()  ; 
} 

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
print STDERR "The following databases will be unmapped_reason updated:\n  ";
print join("\n  ", @dbnames);
print "\ncontinue with update (yes/no)>  ";

my $input = lc(<STDIN>);
chomp($input);
if($input ne 'yes') {
  print "unmapped_reason conversion aborted\n";
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
  push @rows, {'unmapped_reason_id'    => $a[0],
               'summary_description'   => $a[1],
	       'full_description'      => $a[2]};
}
$fh->close();

foreach my $dbname (@dbnames) {
  print STDERR "updating $dbname\n";
  $db->do("use $dbname");
  my $sth = $db->prepare('DELETE FROM unmapped_reason');
  $sth->execute();
  $sth->finish();

  $sth = $db->prepare('INSERT INTO unmapped_reason (unmapped_reason_id,summary_description, full_description)
                       VALUES (?,?,?)');

  foreach my $row (@rows) {
    print $row->{'unmapped_reason_id'}."\n";
    $sth->execute($row->{'unmapped_reason_id'},
		  $row->{'summary_description'},
		  $row->{'full_description'});
  }

  $sth->finish();
}

print STDERR "updates complete\n";


sub usage {
  print STDERR <<EOF

             Usage: update_unmapped_reason_id options
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
  perl update_unmapped_reasons.pl -host ecs1c -file unmapped_reason.txt -user ensadmin -pass secret -dbnames homo_sapiens_core_14_33 -dbnames mus_musculus_core_14_30

  #update all core databases for release 14
  perl update_unmapped_reasons.pl -host ecs2d -file unmapped_reason.txt -user ensadmin -pass secret -release 14

EOF
;
  exit;
}
