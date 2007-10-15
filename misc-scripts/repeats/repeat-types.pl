#
# Repeat classification script
# 
# This script is used to do the repeat classification for web display 
# on newer v32 databases.    
#

use strict;

use DBI;
use Getopt::Long;

my ( $host, $user, $pass, $port, $expression, $dbpattern, $help );

GetOptions( "dbhost|host=s", \$host,
	    "dbuser|user=s", \$user,
	    "dbpass|pass=s", \$pass,
	    "dbport|port=i", \$port,
	    "dbname|dbpattern=s", \$dbpattern,
      "help", \$help
	  );

if($help) {
  usage();
}

if( !$host ) {
  print STDERR "-host argument is required\n";
  usage();
}

if( !$dbpattern ) {
  print STDERR "-dbpattern argument is required\n";
  usage();
}

my $dsn = "DBI:mysql:host=$host";
if( $port ) {
  $dsn .= ";port=$port";
}

my $dbh = DBI->connect( $dsn, $user, $pass );

my @dbnames = map {$_->[0] } @{ $dbh->selectall_arrayref( "show databases" ) };

my @dbs = grep {$_ =~ /$dbpattern/} @dbnames;
die("Haven't found any real dbs from $dbpattern") if(!@dbs);
foreach my $db (@dbs) {
  warn "using $db";
  $dbh->do("use $db");

  print STDERR "  Setting repeat types\n";

  my %mappings = (
    'Low_Comp%' => 'Low complexity regions',
    'LINE%'	=> 'Type I Transposons/LINE',
    'SINE%'	=> 'Type I Transposons/SINE',
    'DNA%'	=> 'Type II Transposons',
    'LTR%'	=> 'LTRs',
    'Other%'	=> 'Other repeats',
    'Satelli%'	=> 'Satellite repeats',
    'Simple%'	=> 'Simple repeats',
    'Other%'	=> 'Other repeats',
    'Tandem%'	=> 'Tandem repeats',
    'TRF%'	=> 'Tandem repeats',
    'Waterman'	=> 'Waterman',
    'Recon'	=> 'Recon',
    'Tet_repeat'	=> 'Tetraodon repeats',
    'MaskRegion'	=> 'Mask region',
    'dust%' => 'Dust',
    'Unknown%'	=> 'Unknown',
    '%RNA'	=> 'RNA repeats',
  );
  foreach (keys %mappings) { 
    $dbh->do(qq(update repeat_consensus set repeat_type = '$mappings{$_}' where repeat_class like '$_')); 
  }

  # type all remaining repeats as unknown
  $dbh->do(qq(update repeat_consensus set repeat_type = 'Unknown' where repeat_type = ''));
  $dbh->do(qq(update repeat_consensus set repeat_type = 'Unknown' where repeat_type = null));
}

print STDERR "All done.\n";

$dbh->disconnect;


sub usage {
  print STDERR <<EOF

This program classifies the repeats stored in a core database into some
somewhat sensible categories.  It does this through a combination of a
repeat.txt file extracted from RepeatMasker repeat libraries and through
some simple pattern matching of the repeat names.

usage: perl repeat-types.pl  [-user <user>] [-port <port>] [-pass <pass>]
               -host <host> -dbpattern <regexp>

example: perl repeat-types.pl -user ensadmin -pass secret -host ecs1g \\
             -port 3306 -dbpattern '^homo_sapiens_(core|vega)_20_34c$'

EOF
;
  exit;
}
