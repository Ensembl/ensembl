#
# These script runs all density creation scripts on all
# databases that match input pattern. It calls all configured
# scripts of this directory with host, user,pass,port and dbname 
# parameter.
#
use strict;
use Getopt::Long;
use DBI;

#
# configure this hash and command line parameter section for 
#  adding new scripts
#

my %density_scripts = 
  ( "gene"   => "gene_density_calc.pl",
    "gc"     => "percent_gc_calc.pl", 
    "repeat" => "repeat_coverage_calc.pl",
    "stats"  => "seq_region_stats.pl",
    "snp"    => "snp_lite_density.pl" );

my @density_scripts = ();

my ( $host, $user, $pass, $port, $dbpattern, $snpstats, $genestats, $stats );

GetOptions( "host=s", \$host,
	    "user=s", \$user,
	    "pass=s", \$pass,
	    "port=i", \$port,
      "genestats", \$genestats,
      "snpstats", \$snpstats,
	    "gene", sub { push( @density_scripts, $density_scripts{'gene'} ); },
	    "gc", sub { push( @density_scripts, $density_scripts{'gc'} ); },
	    "repeat", sub { push( @density_scripts, $density_scripts{'repeat'} ); },
	    "stats", \$stats,
	    "snp", sub { push( @density_scripts, $density_scripts{'snp'} ); },
	    "dbpattern=s", \$dbpattern
	  );

if( !$host || !$user ) {
  usage();
}


if($stats) {
  if($genestats || $snpstats) {
    print STDERR "Only specify one of -stats,  -snpstats or -genestats\m";
    usage();
  }
  push @density_scripts, $density_scripts{'stats'};
} elsif($genestats || $snpstats) {
  my $arg = '';
  $arg .= ' -genestats ' if($genestats);
  $arg .= ' -snpstats ' if($snpstats);
  my $line = $density_scripts{'stats'} . $arg;
  push @density_scripts, $line;
}


my $dsn = "DBI:mysql:host=$host";
my $args = "-host $host ";
if( $user ) {
  $args .= "-user $user ";
}

if( $port ) {
  $args .= "-port $port ";
  $dsn .= ";port=$port";
}

if( $pass ) {
  $args .= "-pass $pass ";
}

if( $host ) {
  $args .= "-host $host ";
}


my $db = DBI->connect( $dsn, $user, $pass );

my @dbnames = map {$_->[0] } @{ $db->selectall_arrayref( "show databases" ) };

for my $dbname ( @dbnames ) {
  if( $dbpattern ) {
    if( $dbname !~ /$dbpattern/ ) {
      next;
    }
  } else {
    print STDERR "$dbname\n";
    next;
  }

  if( ! @density_scripts ) {
    print "Matches $dbname.\n";
  } else {
    for my $scriptname ( @density_scripts ) {
      print( "/usr/local/ensembl/bin/perl $scriptname $args -dbname $dbname\n" );
      system( "/usr/local/ensembl/bin/perl $scriptname $args -dbname $dbname" );
    }
  }
}


sub usage {
  print STDERR "densities_multi_db.pl runs \n  ";
  print STDERR <<EOF;

 with parameters -host [hostname] -user [username with write permit] 
                 -pass [password] -port [portnumber]
 on all databases that match -dbpattern [some_db_regexp]
 on -gene, -repeat, -snp, -stats/-genestats/-snpstats, -gc.

 Without -dbpattern it lists all available databases and does nothing.

EOF
  exit();
}
