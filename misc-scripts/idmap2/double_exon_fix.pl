# this script removes exon duplicates from the ensembl database
# An Exon duplicate is an exon whose start end strand contig phase information
# is the same to another one.

# step one: gather all exon information
# step two build mappings from internal ids to internal ids
# step three change exon_id in exon_transcript 
# remove superfluous exon records 
# remove exon_stable_id records and record which stable_ids were
#  included from which others.

use DBI;
use strict;

my ( %exonHash, %stickies );
my ( $rowsRead, $nrExons );
my %exonTranslate;

my $dbh = DBI->connect( "DBI:mysql:host=ensrv3;database=ensembl110_new_schema_2", "ensadmin", "ensembl");

my $sth = $dbh->prepare( "select exon_id, contig_id, seq_start, seq_end, strand, phase, end_phase, sticky_rank from exon order by exon_id, sticky_rank desc" );
$sth->execute();





while( my $hr = $sth->fetchrow_hashref ) {
  $rowsRead++;
  my $hashkey = join( "\t", ( $hr->{'contig_id'}, $hr->{'seq_start'}, $hr->{'seq_end'},
			      $hr->{'strand'}, $hr->{'phase'}, $hr->{'end_phase'} ));
  push( @{$exonHash{$hashkey}},[ $hr->{'exon_id'}, $hr->{'sticky_rank'} ] ); 
  if( $hr->{'sticky_rank'} > 1 ) {
    $stickies{$hr->{'exon_id'}} = 1;
  }
}

$nrExons = scalar( keys %exonHash );

print ( "Exons in table: ", $rowsRead,"\n" );
print ( "Nonredundant: ",$nrExons,"\n" );


# build equivalence classes
my @equClass;
my %allExons;

my $i = 0;
for my $hashkey ( keys %exonHash ) {
  $equClass[$i] = { 'position' => $hashkey,
		    'exon_id' => undef };

  for my $exon ( @{$exonHash{$hashkey}} ) {
    $allExons{$exon->[0]}->[$exon->[1]-1] = $i; 
  }
  $i++;
}

# make new equclasses by combining old ones
my %combinedEqu;

for my $exonid ( keys %allExons ) {
  if( scalar( @{$allExons{$exonid}} ) > 1 ) {
    my $newEqKey = join( "-",@{$allExons{$exonid}} );
    my $newEq;
    # print STDERR "New equclass made $newEqKey\n"; 

    if( ! exists $combinedEqu{$newEqKey} ) {
      $equClass[$i] = { 'position' => $newEqKey,
			'exon_id' => undef };
      $combinedEqu{$newEqKey} = $i;
      $newEq = $i;
      $i++;
    } else {
      $newEq = $combinedEqu{$newEqKey};
    }

    $allExons{$exonid} = [$newEq];
  }
}


# now give away exon_ids for equClasses
for my $exonid ( keys %allExons ) {
  my $equNo = $allExons{$exonid}->[0];
    if( defined $equClass[$equNo]->{'exon_id'} ) {
      $exonTranslate{$exonid} = $equClass[$equNo]->{'exon_id'}
    } else {
      $equClass[$equNo]->{'exon_id'} = $exonid;
    }
}

print "Translation entries: ",scalar( keys %exonTranslate ),"\n";

for my $exon ( keys %exonTranslate ) {
  print "Map:\t",$exon,"\t",$exonTranslate{$exon},"\n";
}

# dumping_data

$dbh->disconnect;


exit;

sub dumping_data {
  # dump and filter exon, exon_transcript, exon_stable_id
  # should do supporting feature I guess ...
  # this is ONE OFF, 
  $sth = $dbh->prepare( "select * from exon" );
  $sth->execute();
  
  open ( EXON, ">exon.txt" ) or die;
  open ( EXONT, ">exon_transcript.txt" ) or die;
  open ( EXONS, ">exon_stable_id.txt" ) or die;
  
  while( my $arr = $sth->fetchrow_arrayref() ) {
    if( ! defined $exonTranslate{$arr->[0]} ) {
      print EXON join( "\t", @$arr ),"\n";
    }
  }
  
  
  $sth = $dbh->prepare( "select * from exon_transcript" );
  $sth->execute();
  while( my $arr = $sth->fetchrow_arrayref() ) {
    if( defined $exonTranslate{$arr->[0]} ) {
      $arr->[0] = $exonTranslate{$arr->[0]};
    }
    print EXONT join( "\t", @$arr ),"\n";
  }
  
  
  $sth = $dbh->prepare( "select * from exon_stable_id" );
  $sth->execute();
  while( my $arr = $sth->fetchrow_arrayref() ) {
    if( ! defined $exonTranslate{$arr->[0]} ) {
      print EXONS join( "\t", @$arr ),"\n";
    }
  }
  
  close( EXON );
  close( EXONT );
  close( EXONS );
}


