use Data::Dumper;

# change name of input files further down. 
# they contain sorted alphabetically by clone data
# on contigs. Finds out, which clones to delete and 
# which to reload.


open (  ENS, "f17Dump/cloneContig.srt" );
open ( NCBI, "ncbiSortClones.txt" );
open ( WRI, ">clonesToLoad.txt" );
open ( DEL, ">clonesToDelete.txt" );

$cloneEns = nextInFile( \*ENS );
$cloneNCBI = nextInFile( \*NCBI );

while( defined $cloneEns && defined $cloneNCBI ) {
  $cmp = cmpClones( $cloneEns, $cloneNCBI );  
  # print Dumper( $cloneEns, $cloneNCBI ),"\n";
  # print $cloneEns->[0][0],"\n";
  # print $cloneNCBI->[0][0],"\n";
  # print $cmp,"\n";
  
  if( $cmp eq "equal" ) {
    $cloneEns = nextInFile( \*ENS );
    $cloneNCBI = nextInFile( \*NCBI );
  } elsif( $cmp eq "update" ) {
    # print Dumper( $cloneEns, $cloneNCBI ),"\n";
    $update += scalar( @$cloneNCBI );
    
    writeNcbi( $cloneNCBI );
    deleteEns( $cloneEns );
    $cloneEns = nextInFile( \*ENS );
    $cloneNCBI = nextInFile( \*NCBI );
  } elsif( $cmp eq "lower" ) {
    # ncbi clonename is bigger
    deleteEns( $cloneEns );
    $cloneEns = nextInFile( \*ENS );
  } elsif( $cmp eq "higher" ) {
    # ensembl clonename is bigger
    writeNcbi( $cloneNCBI );
    $cloneNCBI = nextInFile( \*NCBI );
  } 
}
    
if( defined $cloneNCBI ) {
  while( defined $cloneNCBI ) {
    writeNcbi( $cloneNCBI );
    $cloneNCBI = nextInFile( \*NCBI );
  }
} elsif( defined $cloneEns ) {
  while( defined $cloneEns ) {
    deleteEns( $cloneEns );
    $cloneEns = nextInFile( \*ENS );
  }
}

close WRI;
close DEL;
close ENS;
close NCBI;

print( "Update contigs were $update\n" );

sub nextInFile {
  my $fh = shift;
  my $first = 1;
  my ( $lastpos, $clonename );
  my $result;
#  print STDERR ( "nextNcbi called \n" );
  
  while( <$fh> ) {
    my @arr = split( /[\t ]/,$_ );
    
    if( $first ) {
      $lastpos = tell( $fh );
      $clonename = $arr[0];
      push( @$result, \@arr );
      $first = 0;
    } else {
      if( $arr[0] eq $clonename ) {
        push( @$result, \@arr );
	$lastpos = tell( $fh );
      } else {
        seek( $fh, $lastpos, 0 );
	last;
      }
    }
  }
  return $result;
  print "$result\n";
}


# equal
# update
# first lower
# second lower
sub cmpClones {
  my $clone1 = shift;
  my $clone2 = shift;

  my $cmp;
  if(( $cmp = ( $clone1->[0][0] cmp $clone2->[0][0] )) != 0 ) {
    if( $cmp == -1 ) {
      return "lower";
    } else {
      return "higher";
    }
  }
  if( $clone1->[0][1] != $clone2->[0][1] ) {
    return "update";
  }
  my @clone1Sorted = sort { $a->[3] <=> $b->[3] } @$clone1;
  my @clone2Sorted = sort { $a->[3] <=> $b->[3] } @$clone2;
  my ( $contig1, $contig2 );
  
  if( scalar( @clone1Sorted ) != scalar( @clone2Sorted )) {
    return "update";
  }
  
  while( $#clone1Sorted > -1 ) {
    $contig1 = shift( @clone1Sorted );
    $contig2 = shift( @clone2Sorted );
    if( $contig1->[3] != $contig2->[3] ||
        $contig1->[4] != $contig2->[4] ) {
      return "update";
    }
  }
  return "equal";
}

sub deleteEns {
  my $clone = shift;
  print DEL $clone->[0][0],"\n";
}

sub writeNcbi {
  my $clone = shift;
  print WRI $clone->[0][0],"\n";
}

