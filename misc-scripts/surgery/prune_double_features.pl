use DBI;
use strict;

my $dbh = DBI->connect( "DBI:mysql:host=ensrv3;database=homo_sapiens_core_120", "ensro", "");

open( PRUNED, ">/scratch1/ensembl/arne/feature.txt" ) or die;


my $sth = $dbh->prepare( "select internal_id from contig" );
$sth->execute;
my $prune = 0;
my $read = 0;
my @contigs;

while( my $arrref = $sth->fetchrow_arrayref() ) {
  push( @contigs, $arrref->[0] );
}

$sth = $dbh->prepare( "select  id, contig, seq_start, seq_end, score, strand,
                                 analysis, name, hstart, hend, hid, evalue, perc_id, phase, end_phase
                            from feature where contig = ? ");

for( my $i = 0; $i <= $#contigs; $i++ ) {
  my $contig = $contigs[$i];
  my %linehash = ();
  
  $sth->execute( $contig );
  while( my $arrref = $sth->fetchrow_arrayref ) {
    $read++;

    my $hashkey = join( "\t", ((@$arrref)[1..14]) );
#    print $hashkey,"\n";

    if( exists $linehash{$hashkey} ) {
      $prune++;
    } else {
      $linehash{$hashkey} = $arrref->[0];
      print PRUNED $arrref->[0],"\t$hashkey\n";
    }
  }
  if(( $i +1 ) % 500 == 0 ) {
    print STDERR "$i/",$#contigs," contigs, $prune/$read features\n";
  }
}

print STDERR "$prune/$read features\n";

$dbh->disconnect;
