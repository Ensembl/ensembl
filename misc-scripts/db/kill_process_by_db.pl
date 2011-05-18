# kill processes by db pattern
# Author: Monika Komorowska
# Date : May 2011


use strict;
use DBI;

use Getopt::Long;

my ( $host, $user, $pass, $port, $dbpattern );

GetOptions( "host|h=s", \$host,
	    "user|u=s", \$user,
	    "pass|p=s", \$pass,
	    "port|P=i", \$port,
	    "dbpattern|pattern=s", \$dbpattern,
	  );

if( !$host || !$user || !$pass || !$dbpattern ) {
  usage();
}

my $dsn = "DBI:mysql:host=$host";
if( $port ) {
  $dsn .= ";port=$port";
}

my $db = DBI->connect( $dsn, $user, $pass );

my $sth=$db->prepare("SHOW FULL PROCESSLIST") || die $DBI::err.": ".$DBI::errstr;

$sth->execute || die DBI::err.": ".$DBI::errstr;

my $ref;

print "Id \t User \t Host \t db \t Command \t Time \t State \t Info\n";
my @proc_ids;

while ( $ref = $sth->fetchrow_hashref() ) {
    if ( $$ref{'db'} =~ /$dbpattern/ ) {
 	push (@proc_ids, $$ref{'Id'});
    	print "$$ref{'Id'} \t $$ref{'Host'}  \t $$ref{'db'}  \t $$ref{'Command'}  \t $$ref{'Time'}  \t $$ref{'State'}  \t $$ref{'Info'}\n";
    }
}

my $proc_id_count = @proc_ids;
if ($proc_id_count > 0) {
	print "Are you sure you want to kill process(es) listed above? (y/n)\n";
} else {
	print "No processes found for db pattern $dbpattern\n";
}
my $decision = <>;

if ($decision == 'y') {
  my $killed_count = 0;
  foreach my $proc_id (@proc_ids) {
  	if ( $db->do("KILL $proc_id") ) {
		$killed_count ++;
	} else { print $DBI::errstr; }
  }

  print "$killed_count procesess were killed\n";

}


sub usage {
  print STDERR <<EOF

             Usage: kill_process_by_db options
 Where options are: -host hostname 
                    -user username 
                    -pass password 
                    -port port_of_server optional
                    -dbpattern regular expression that the database name has to match
EOF
;
  exit;
}

