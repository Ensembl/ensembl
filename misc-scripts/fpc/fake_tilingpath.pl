# read in current fpc data for clone information
# ( name, institute, accession )

# get from static_golden_path contig clone 
# contigstart, clonename, offset, length
# fpc contigname

# project each contig in golden path
# can the clones project conflict free?

use DBI;
use strict;

my ( %contig, %clones );

sub reading_clone_info {
  my $db = DBI->connect
    ( "DBI:mysql:database=homo_sapiens_maps_130;host=ecs1d.sanger.ac.uk",
      "ensadmin","ensembl" );

  my $sth = $db->prepare( "
   SELECT embl_id, clone_name, organisation
     FROM Fpc_Clone cl
  " );

  $sth->execute();

  while( my $arrRef = $sth->fetchrow_arrayref ) {
    $clones{$arrRef->[0]} =  [ $arrRef->[1], $arrRef->[2] ];
  }
  $db->disconnect();
}


sub reading_assembly {
  my $db = DBI->connect
    ( "DBI:mysql:database=homo_sapiens_core_130;host=ecs1d.sanger.ac.uk",
      "ensro", "" );
  my $sth = $db->prepare( "
    SELECT sgp.fpcctg_name, sgp.chr_name, sgp.chr_start, sgp.chr_end, sgp.raw_start,
           sgp.raw_end, sgp.raw_ori, co.offset, co.length, cl.embl_id
      FROM static_golden_path sgp, contig co, clone cl
     WHERE sgp.raw_id = co.internal_id
       AND cl.internal_id = co.clone
  " );
  
  $sth->execute();

  while( my $hashRef = $sth->fetchrow_hashref() ) {
    my ( $chr, $ctg_start, $ctg_end, $fpc_ctg );
    if( $hashRef->{raw_ori} == 1 ) {
      $ctg_start = $hashRef->{chr_start} - $hashRef->{raw_start} - 1;
      $ctg_end = $ctg_start + $hashRef->{length} - 1;
    } else {
      $ctg_start = $hashRef->{chr_start} - ( $hashRef->{length} - $hashRef->{raw_end} );
      $ctg_end = $ctg_start + $hashRef->{length} - 1;
    }
    $chr = $hashRef->{chr_name};
    $fpc_ctg = $hashRef->{fpcctg_name};

    $contig{$hashRef->{embl_id}."-".$hashRef->{offset}} = [ $chr, $fpc_ctg, $ctg_start, $ctg_end, $hashRef->{embl_id} ];
  }

  $db->disconnect();
 }


# find out if all clones in golden path have info in 
# fpc tables. (info = organisation and name)
sub all_clones_info {
  my %cloneinfo_missing;

  for my $key ( keys %contig ) {
    my $clone = $contig{$key}->[4];
    if( ! exists $clones{$clone} ) {
      $cloneinfo_missing{$clone} = 1;
    }
  }

  for my $clone ( keys %cloneinfo_missing ) {
    print "$clone has no info.\n";
  }
  print scalar( keys( %cloneinfo_missing )), " clones have no info.\n";

}


# find chrstart end for the clones
# min max of all contigs in the clone
sub project_clones {
  my ( %clone_project );
  for my $key ( keys %contig ) {
    my $clone = $contig{$key}->[4];
    if( exists $clone_project{$clone} ) {
      my ( $chr, $fpcctg, $start, $end ) = @{$clone_project{$clone}};
      if(( $chr ne $contig{$key}->[0] ) ||
	 ( $fpcctg ne $contig{$key}->[1] )) {
	print STDERR "Bugger!\n $clone placed in diffferent contigs\n";
      }
      $start = $contig{$key}->[2] < $start ? $contig{$key}->[2] : $start;
      $end = $contig{$key}->[3] > $end ? $contig{$key}->[3] : $end;
      $clone_project{$clone} = [ $chr, $fpcctg, $start, $end, $clone ];
    } else {
      $clone_project{$clone} = [ $contig{$key}->[0], $contig{$key}->[1],
				 $contig{$key}->[2], $contig{$key}->[3], $clone ];
    }
  }
  return \%clone_project;
}



# the ensembl tables fpc_clone and fpc_contig
# need to be filled with the clone projection
sub fpc_output_format {
  my $clone_project = shift;
  my %contigCrash;
  my $hiLengthClones = 0;

  my @sortClones = sort { ( $a->[0] cmp $b->[0] ) ||
			    ( $a->[2] <=> $b->[2] ) } ( values %$clone_project );
  
  my ( $ctgNo, $ctgName, $ctgStart, $ctgEnd, $chr );
  $ctgName = undef;
  $ctgNo = 0;

  for my $clone ( @sortClones ) {
    if( $ctgName ne $clone->[1] ) {
      # new contig
      if( defined $ctgName ) {
	print join( "\t", ( "ctg", $ctgNo, $ctgName, $ctgStart, $ctgEnd-$ctgStart+1, $chr )),"\n";
	$contigCrash{$ctgName} = 1;
      }
      $ctgNo++;
      $ctgName = $clone->[1];
      $ctgStart = $clone->[2];
      $ctgEnd = $clone->[3];
      $chr = $clone->[0];
      if( exists $contigCrash{$ctgName} ) {
	print STDERR "CRASH: $ctgName overlapped with another contig\n";
      }
    }
    
    my ( $cloneName, $org, $embl );
    $embl = $clone->[4];
    if( exists $clones{$embl} ) {
      $cloneName = $clones{$embl}->[0];
      $org = $clones{$embl}->[1];
    } else {
      $cloneName = "?";
      $org = "?";
    }

    print join( "\t", ( "clone", $embl, $cloneName, $org, $ctgNo, 0, $clone->[3]-$clone->[2]+1,
			$clone->[3]-$clone->[2]+1, $clone->[2]-$ctgStart+1 )),"\n";
    if( $clone->[3] - $clone->[2] > 500000 ) {
      print STDERR "$embl has length ",$clone->[3]-$clone->[2]+1,"\n";
      $hiLengthClones++;
    }
  }
  print STDERR "There are $hiLengthClones clones longer than 500K in the assembly.\n";
}

reading_clone_info();
reading_assembly();


all_clones_info();
my $clone_project = project_clones();
  
# need output in fpc-map database format
fpc_output_format( $clone_project );





__END__


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

