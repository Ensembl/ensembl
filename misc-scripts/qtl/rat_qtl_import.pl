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


# the rgd to ensembl table
my $rgd_ens_map = $db->selectall_arrayref
     ("select name, marker_id from marker_synonym " .
      "where source = 'RGD' or source = 'RGD_NUM' or source = 'rgd' or  source = 'rgdgene'" );

my %rgd_ens_map = map { $_->[0], $_->[1] } @$rgd_ens_map;

# need source_database, source_primary, trait, lod, flank1, flank2, peak

my $qtl_id = 0;
my $qtl_syn_id = 0;

QTL:
for my $qtl ( @all_qtl ) {
  my ($qtl_rgd_id, $species, $qtl_symbol,$qtl_name,
      $peak_offset, $chromosome, $lod_score, $p_value, $variance,
      $fmark1_rgd_id, $fmark1_symbol, $fmark2_rgd_id, $fmark2_symbol,
      $pmark_rgd_id, $pmark_symbol, $trait_name, $subtrait_name,
      $trait_desc, $curated_ref_rgd_id, $curated_ref_pubmed_id,
      $uncurated_ref_pubmed_id, $ratmap_id, $locus_link_id, @other) =
        @$qtl;

  if(!$qtl_name) {
    debug("QTL $qtl_rgd_id has no name, discarding");
    next QTL;
  }


  my $fmark1_id = $rgd_ens_map{$fmark1_rgd_id};
  my $fmark2_id = $rgd_ens_map{$fmark2_rgd_id};
  my $pmark_id = $rgd_ens_map{$pmark_rgd_id};

  if(!$fmark1_id) {
    debug("Flanking Marker1 RGD:$fmark1_rgd_id not found for Qtl $qtl_symbol");
    $fmark1_id = 'null';
  }

  if(!$fmark2_id) {
    debug("Flanking Marker2 RGD:$fmark2_rgd_id not found for Qtl $qtl_symbol");
    $fmark2_id = 'null';
  }

  if(!$pmark_id) {
    debug("Peak Marker RGD:$pmark_rgd_id not found for Qtl $qtl_symbol");
    $pmark_id = 'null';
  }

  $lod_score ||= 'null';

  $qtl_id++;

  $qtl_name    = $db->quote($qtl_name);
  $qtl_symbol  = $db->quote($qtl_symbol);

  print "INSERT INTO qtl (qtl_id, trait, lod_score, flank_marker_id_1,"
                       . "flank_marker_id_2, peak_marker_id)\n" .
        "     VALUES ($qtl_id, $qtl_name, $lod_score, $fmark1_id, $fmark2_id, "
     .  "$pmark_id);\n";

  $qtl_syn_id++;
  print "INSERT INTO qtl_synonym (qtl_synonym_id, qtl_id, source_database, " .
                                 "source_primary_id)\n" .
   "     VALUES ($qtl_syn_id, $qtl_id, 'rat genome database', $qtl_rgd_id);\n";

  if($ratmap_id) {
    $ratmap_id   = $db->quote($ratmap_id);
    $qtl_syn_id++;
    print "INSERT INTO qtl_synonym (qtl_synonym_id, qtl_id, source_database,"
                                 ."  source_primary_id)\n" .
          "     VALUES ($qtl_syn_id, $qtl_id, 'ratmap', $ratmap_id);\n";
  }

}

$db->disconnect();

sub debug {
  my $message = shift;
  print STDERR ( $message,"\n" ) if $verbose;
}


sub usage {
  print STDERR "usage:\n" .
               "  perl rat_qtl_import.pl -host <host> -user <user> " .
               "-dbname <dbname> -qtlfile <qtlfile> [-pass <pass>] " .
               "[-port <port>] [-verbose] > qtl.sql\n" .
               "  cat qtl.sql | mysql ...\n";
  exit;
}
