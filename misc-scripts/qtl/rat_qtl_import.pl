# read the qtl information
# create a list of marker ids from rgd
# find them in rat database
# create a hash to ensembl marker ids.

# upload the stuff in qtl table


use strict;
use DBI;

use Getopt::Long;

my ( $host, $user, $pass, $port, $dbname, $qtlfile );
my $verbose = 0;

GetOptions( "host=s", \$host,
	    "user=s", \$user,
	    "pass=s", \$pass,
	    "port=i", \$port,
	    "dbname=s", \$dbname,
	    "qtlfile=s", \$qtlfile,
	    "verbose", \$verbose
	  );

if( !$host ) {
  usage();
}

my $dsn = "DBI:mysql:host=$host;dbname=$dbname";
if( $port ) {
  $dsn .= ";port=$port";
}

my $db = DBI->connect( $dsn, $user, $pass );


open( FH, $qtlfile );
my $first = 1;
my @all_qtl;

while( <FH> ) {
  chomp;
  # skip header line
  if( $first ) {
    $first = 0;
    next;
  }

  my @cols = split( "\t", $_ );
  
  push( @all_qtl, \@cols );
}

my %rat_marker_hash;

for my $qtl ( @all_qtl ) {
  $rat_marker_hash{ $qtl->[2] } = 0;
  $rat_marker_hash{ $qtl->[3] } = 0;
  $rat_marker_hash{ $qtl->[4] } = 0;
}

# the rgd to ensembl table
my $rgd_ens_map = $db->selectall_arrayref( "select name, marker_id from marker_synonym where source = 'rgd' or source = 'rgdgene'" );

my %rgd_ens_map = map { $_->[0], $_->[1] } @$rgd_ens_map;


# need source_database, source_primary, trait, lod, flank1, flank2, peak
QTL:
for my $qtl ( @all_qtl ) {
  my @ens_marker_ids = ();
  my $tmp;

  for my $i ( 2,3,4 ) {
    if( $qtl->[$i] ) {
      if( ! exists $rgd_ens_map{ $qtl->[$i] } ) {
#	next QTL;
	debug("Marker RGD:".$qtl->[$i]." not found for Qtl ".$qtl->[ );
	push( @ens_marker_ids, "\\N" );
      } else {
	push( @ens_marker_ids, $rgd_ens_map{ $qtl->[$i] } );
      }
    } else {
      push( @ens_marker_ids, "\\N" );
    }
  }

  print join( "\t", "rat genome database", $qtl->[0], $qtl->[5], $qtl->[1]?$qtl->[1]:"\\N",
	      @ens_marker_ids ),"\n";
	      
}

sub debug {
  my $message = shift;
  print STDERR ( $message,"\n" ) if $verbose;
}

 


sub usage {
  exit;
}
