# parse jim kent clonePos file
# check if all clones are in map db

# check if some clones in map db are not in jims file
use DBI;


sub reading_the_db {
  my $db = DBI->connect
    ( "DBI:mysql:database=maps100_arne;host=ecs1a.sanger.ac.uk",
      "ensadmin","" );

  my $sth = $db->prepare( "
   SELECT embl_id, clone_name, organisation, contig_name
     FROM Fpc_Clone cl, Fpc_Contig co
    WHERE cl.contig_id = co.contig_id
  " );

  $sth->execute();

  while( $arrRef = $sth->fetchrow_arrayref ) {
    $clones{$arrRef->[0]} =  [ $arrRef->[1], $arrRef->[2], $arrRef->[3]];
  }
}

reading_the_db();

open( JIM, shift ) or die( "No Jim Kent File" );

my $first = 1;

while( <JIM> ) {
  if( $first ) {
    $first = 0;
    next;
  }
  $line = $_;
  split( '\t', $_ );
  my $cname = $_[0];
  my $cstart = $_[4];
  my $cend = $_[5];
  my $cchrom = $_[3];
  
  # leaving out random contigs
  if( $cchrom =~ /random/ ) {
    next;
  }
  $cchrom =~ s/^chr//g;

  ( $cname ) = ( $cname =~ /([^.]*)/);
#  if( defined $clones{$cname} ) {
#    my $contig = $clones{$cname};
#    $contigs{$contig} = 2;
#  }
  push( @jclones, { 'embl' => $cname, 
		    'chrom' => $cchrom, 
		    'start' => $cstart, 
		    'end' => $cend } );
  if( $cstart =~ /chrom/ ) {
    print STDERR $line;
    print STDERR join( " ", ($cname,$cchrom, $cstart, $cend ));
    print STDERR "\n";
  }
}


# if it has info from fpc db, set it
foreach $clone ( @jclones ) {
  my $fpcinfo = $clones{$clone->{'embl'}};
  if( defined $fpcinfo ) {
    $clone->{'name'} = $fpcinfo->[0];
    $clone->{'orga'} = $fpcinfo->[1];
    $clone->{'contig'} = $fpcinfo->[2];
  }
}

# compile fpc_contigs from the clones
# foreach $clone ( @jclones ) {
#  my $contig = $clone->{'contig'};
#  if( defined $contig ) {
#    if( defined $fpcContigs{$contig} ) {
#      $fpcContig = $fpcContigs{$contig};
#    } else {
#      $fpcContig = { 'start' => 1000000000,
#		     'end' => 0,
#		   };
#      $fpcContigs{$contig} = $fpcContig;
#    }
#
#    $fpcContig->{'start'} = $clone->{'start'} < $fpcContig->{'start'} ?
#      $clone->{'start'} : $fpcContig->{'start'};
#    $fpcContig->{'end'} = $clone->{'end'} > $fpcContig->{'end'} ?
#      $clone->{'end'} : $fpcContig->{'end'};
#    push( @{$fpcContig->{'clones'}}, $clone );
#  }
#}


# compile overlap contigs
@overlapContigs = ();
foreach $clone ( @jclones ) {
  my $oc = { 'start' => $clone->{'start'},
	     'end' =>$clone->{'end'},
	     'chrom' =>$clone->{'chrom'},
	     'clones' => [ $clone ]
	   };
  push( @overlapContigs, $oc );
}

@newOverlaps = sort {
  if( $a->{'chrom'} cmp $b->{'chrom'}  ) {
    $a->{'chrom'} cmp $b->{'chrom'};
  } else {
    $a->{'start'} <=> $b->{'start'};
  } } @overlapContigs;


@overlapContigs = @newOverlaps;
@newOverlaps = ();

$no = shift( @overlapContigs );

while( @overlapContigs ) {
  $oc = shift( @overlapContigs );

  if( $oc->{'chrom'} eq $no->{'chrom'} ) {
    if( $oc->{'start'} <= $no->{'end'}+1 &&
	$oc->{'end'}+1 >= $no->{'start'} ) {
      if( $no->{'start'} > $oc->{'start' } ) {
	$no->{'start'} = $oc->{'start'};
      }
      if( $no->{'end'} < $oc->{'end' } ) {
	$no->{'end'} = $oc->{'end'};
      }
      push( @{$no->{'clones'}}, @{$oc->{'clones'}} );
      next;
    }
  }
  push( @newOverlaps, $no );
  $no = $oc;
}
push( @newOverlaps, $no );

@overlapContigs = @newOverlaps;

# check if all contigs agree in the fpc contig
# (or dont disagree)
$unknown = 1;

foreach $oc ( @overlapContigs ) {
  $fpcContig = undef;

  foreach $clone ( @{$oc->{'clones'}} ) {
    if( defined $fpcContig && defined $clone->{'contig'} &&
	$fpcContig ne $clone->{'contig'} ) {
      print STDERR "Contig mismatch clone ", $clone->{'embl'},"\n";
    }
    if( ! defined $fpcContig && defined $clone->{'contig'} ) {
      $fpcContig = $clone->{'contig'};
    }
  }

  if( defined $fpcContig ) {
    $oc->{'contig'} = $fpcContig;
  } else {
    $oc->{'contig'} = "ukn".$unknown;
    $unknown++;
  }
  foreach $clone ( @{$oc->{'clones'}} ) {
    $clone->{'contig'} = $oc->{'contig'};
  }
}

# merge adjacent contigs with same fpc name
@newOverlaps = ();
$no = shift( @overlapContigs );
while( @overlapContigs ) {
  $oc = shift( @overlapContigs );
  if( $no->{'contig'} eq $oc->{'contig'} ) {
    if( $no->{'start'} > $oc->{'start' } ) {
      $no->{'start'} = $oc->{'start'};
    }
    if( $no->{'end'} < $oc->{'end' } ) {
      $no->{'end'} = $oc->{'end'};
    }
    push( @{$no->{'clones'}}, @{$oc->{'clones'}} );
    next;
  } else {
    push( @newOverlaps, $no );
    $no = $oc;
  }
}

push( @newOverlaps, $no );

@overlapContigs = @newOverlaps;

$contigNum = 1;
open( FPC_CONTIG, ">fpc_contig.txt" ) or die( "Cant write" );
open( FPC_CLONE, ">fpc_clone.txt" ) or die( "Cant write" );

foreach $oc ( @overlapContigs ) {
  print FPC_CONTIG join( "\t", ($contigNum, $oc->{'contig'}, $oc->{'start'}+1,$oc->{'end'}-$oc->{'start'}+1, $oc->{'chrom'} )),"\n";
  foreach $clone ( @{$oc->{'clones'}} ) {
    print FPC_CLONE join
      ( "\t", ( $clone->{'embl'}, $clone->{'name'},
		$clone->{'orga'}, $contigNum, 0, 
		$clone->{'end'}-$clone->{'start'}+1,
		$clone->{'end'}-$clone->{'start'}+1,
		$clone->{'start'}-$oc->{'start'}+1 )),"\n";
  }
  $contigNum++;
}

close( FPC_CLONE );
close( FPC_CONTIG );

__END__


print "We have ",scalar( @overlapContigs ), " overlapContigs.\n";

foreach $oc ( @overlapContigs ) {
  print $oc->{'chrom'}," ",$oc->{'start'}," ",$oc->{'end'},"\n";
  foreach $clone ( @{$oc->{'clones'}} ) {
    while( ($k,$v) = each %$clone ) {
      print " ",$k," ",$v,"\n";
    }
    print"\n";
  }
  print "\n";
}




__END__


open ( FPC, shift );

while( <FPC> ) {
  split( '\t', $_ );
  $cname = $_[0];
  if( $cname ) {
    $fclones{$cname} = 1;
  }
}

foreach $clone ( @jclones ) {
  if( ! defined $fclones{$clone->[0]} ) {
    print $clone->[0]," undefined in fpc.\n" ;
  }
}

