# apply given sql statement to all databases on a given host
# you can specify a name pattern for the database
# it displays results for select statements
# Author: Arne Stabenau
#  Date : 21.02.2003



use strict;
use DBI;

use Getopt::Long;

my ( $host, $user, $pass, $port, $expression, $dbpattern );

GetOptions( "host=s", \$host,
	    "user=s", \$user,
	    "pass=s", \$pass,
	    "port=i", \$port,
	    "expr=s", \$expression,
	    "dbpattern=s", \$dbpattern
	  );

if( !$host ) {
  usage();
}

my $dsn = "DBI:mysql:host=$host";
if( $port ) {
  $dsn .= ";port=$port";
}

my $db = DBI->connect( $dsn, $user, $pass );

my @dbnames = map {$_->[0] } @{ $db->selectall_arrayref( "show databases" ) };

for my $dbname ( @dbnames ) {
  if( $dbpattern ) {
    if( $dbname !~ /$dbpattern/ ) {
      next;
    }
  }

  $db->do( "use $dbname" );
  print STDERR "$dbname\n";
  if( $expression =~ /^\s*select/i ||
      $expression =~ /^\s*show/i ||
      $expression =~ /^\s*describe/i ) {
    my $res = $db->selectall_arrayref( $expression );
    my @results = map { join( " ", @$_ ) } @$res ;
    for my $result ( @results ) {
      print STDERR "  Result: ",$result,"\n";
    }
  } else {
    $db->do( $expression );
    print STDERR "  done.\n";
  }
}


sub usage {
  print STDERR <<EOF
             Usage: multidb_sql options
 Where options are: -host hostname 
                    -user username 
                    -pass password 
                    -port port_of_server optional
                    -dbpattern regular expression that the database name has to match
                    -expr sql statement you want to execute. Has to be a select.
EOF
;
}

