# this script tales filled qtl table and tries to calculate 
# qtl_features. It uses marker_feature table.

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Map::DBSQL::QtlAdaptor;

use Getopt::Long;


my ( $host, $user, $pass, $port, $dbname, $verbose);

$verbose = 0;

GetOptions( "host=s", \$host,
	    "user=s", \$user,
	    "pass=s", \$pass,
	    "port=i", \$port,
	    "dbname=s", \$dbname,
	    "verbose", \$verbose
	  );


if( !$host ) {
  usage();
}


my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new
  (
   -host => $host,
   -user => $user,
   -pass => $pass,
   -dbname => $dbname,
   -port => $port 
  );




my $qtlAdaptor = $db->get_QtlAdaptor();
my $qtls = $qtlAdaptor->fetch_all;

QTL:
for my $qtl ( @$qtls ) {
  my $marker;
  my $features;
  my @positions;

  $marker = $qtl->flank_marker_1;
  if( $marker ) {
    $features = $marker->get_all_MarkerFeatures();
    if( scalar( @$features ) <= 1 ) {
      push( @positions, @$features );
    } else {
      debug( "Non unique flanking marker ".$marker->display_MarkerSynonym->name()." for qtl ".
	     $qtl->source_primary_id());

      # position of marker not unique, dont place qtl
      next;
    }
  }

  $marker = $qtl->flank_marker_2;
  if( $marker ) {
    $features = $marker->get_all_MarkerFeatures();
    if( scalar( @$features ) <= 1 ) {
      push( @positions, @$features );
    } else {
      # position of marker not unique, dont place qtl
      debug( "Non unique flanking marker ".$marker->display_MarkerSynonym->name()." for qtl ".
	     $qtl->source_primary_id());

      next;
    }
  }

  $marker = $qtl->peak_marker;
  if( $marker ) {
    $features = $marker->get_all_MarkerFeatures();
    if( scalar( @$features ) <= 1 ) {
      push( @positions, @$features );
    } else {
      # position of marker not unique, dont place qtl
      debug( "Non unique peak marker ".$marker->display_MarkerSynonym->name()." for qtl ".
	     $qtl->source_primary_id());
      next;
    }
  }

  
  if( scalar( @positions )) {
    my ( $chr_name, $chr_start, $chr_end, $chr_id );
    
    for my $feature ( @positions ) {
      my $empty_slice = Bio::EnsEMBL::Slice->new
	( 
	 -empty => 1,
	 -adaptor => $db->get_SliceAdaptor(),
	 -assembly_type => $db->assembly_type()
	);
      $feature->transform( $empty_slice );

      if( $chr_name && $feature->contig->chr_name() ne $chr_name ) {
	# inconsistent chromsome skip this one
	next QTL;
      } 
      $chr_name = $feature->contig->chr_name();
      $chr_id = $feature->contig->get_Chromosome()->dbID();

      if( ! defined $chr_start ) {
	$chr_start = $feature->start();
      } elsif( $feature->start() < $chr_start ) {
	$chr_start = $feature->start();
      }
      if( ! defined $chr_end ) {
	$chr_end = $feature->end();
      } elsif( $feature->end() > $chr_end ) {
	$chr_end = $feature->end();
      }
    }
    
    if( $chr_end - $chr_start < 1000 ) {
      debug( "Qtl ".$qtl->source_primary_id()." is shorter than 1kb." );
      $chr_end += 10_000;
      $chr_start -= 10_000;
    } elsif( $chr_end - $chr_start > 20_000_000 ) {
      debug( "Qtl ".$qtl->source_primary_id()." covers more than 20MB" );
      next;
    }
    print join( "\t", ( $chr_id, $chr_start, $chr_end, $qtl->dbID(), 50 )),"\n";
  }
}

sub debug {
  my $string = shift;
  print STDERR ($string,"\n") if $verbose;
}


sub usage {
  print STDERR <<EOF

  Usage: perl qtl_feature_calculation
          -host hostname
          -user username
          -pass password
          -port portnumber of EnsEMBL SQL server
          -dbname name of EnsEMBL database
          -verbose optional more output

EOF
;
  exit();
}
