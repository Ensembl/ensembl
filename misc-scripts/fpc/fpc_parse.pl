# input, the slightly edited ALLmap file from an Email
# output
# information about Fpc_Clone table entries and
# Fpc_Contig table entries

$allfile = shift;

$mingap = 150000;

open( FH, $allfile ) or die( "Couldnt open input file" );
open( CO, ">fpc_contig.txt" ) or 
  die( "Couldnt open fpc_contig.txt output file." );
open( CL, ">fpc_clone.txt" ) or 
  die( "Coulnt open fpc_clone.txt output file." );

$cloneHash = {};
$contigHash = {};
$contigNum = 0;

%chrlength = ( "1", 263000058, "10", 143999990, "11", 144000020, "12", 142999999, "13", 113999999, "14", 108999990, "15", 106000020, "16", 98000000, "17", 92000000, "18", 85000029, "19", 66999999, "2", 254999998, "20", 72000009, "21", 50000000, "22", 50000000, "3", 214000030, "4", 203000002, "5", 194000011, "6", 183000009, "7", 171000011, "8", 155000010, "9", 144999990, "X", 164000011, "Y", 50810739 );


while( <FH> ) {
  /^\#/ && next;
  $clone = {};

  if( /^start.*\.human\.(\S+) (\S+)/ ) {
    $chrom = $1;
    if( $2 eq "ORDERED" ) {
      $ordered = 1;
    } else {
      $ordered = 0;
    }
    next;
  }

  if( /^@\s*$/ ) {
    if( defined $lastContig ) {
      # contig_done( $chrom, $lastContig, $ordered );
      $lastContig = undef;
    }
    next;
  }

  if( /^end.*\.human\./ ) {
    next;
  }

  if( /^@\s+(\S+)\s+(\S+)/ ) {
    # contig done
    $contigHash->{$1}->{start} = $2;
    # contig_done( $chrom, $lastContig, $ordered, $1, $2 );
    $lastContig = undef;
    next;
  }

  if( /^@/  ) {
    # no start defined probably no name either
    $lastContig = undef;
    next;
  }

  chomp;
  if( /\t/ ) {
    split( '\t', $_ );
    map { s/^\s*//g;s/\s*$//g } @_;
  } else {
    split;
  }


  if( $_[0] =~ /^COMMITTED/) {
    $embl = 0;
  } else { 
    $embl = 1;
    $clone->{embl} = $_[0];
    if( defined $cloneHash->{$_[0]} ) {
      print("Double clone $_[0]\n" )
    }
    $cloneHash->{$_[0]} = 1;
  }

  if(!( $_[1] =~ /UNK/ || $_[1] =~ /multi/ )) {
    $clone->{name}  = $_[1];
  }
  $clone->{institute} = $_[2];
  $clone->{contig} = $_[3];
  $lastContig = $_[3];

  $clone->{start}  = $_[4]*1000;
  if( $embl ) {
    $clone->{length} = $_[6];
  }
  submit_clone( $chrom, $ordered, $clone );
  $clone = undef;
}

for (keys %$contigHash ) {
  if( $contigHash->{$_}->{ordered} ) {
    push( @{$chrom{$contigHash->{$_}->{chrom}}}, $contigHash->{$_} );
  }
}

for $chrom (keys %chrom) {
#  my @newlist = sort {$a->{start} <=> $b->{start}} @{$chrom{$chrom}};
#  $chrom{$chrom} = \@newlist;
#  recalc_starts_standard( $chrom, $chrom{$chrom} );
  recalc_eye_candy( $chrom, $chrom{$chrom} );
}


for (keys %$contigHash ) {
  $contigRef = $contigHash->{$_};
  if( ! defined $contigRef->{name} ) {
    print STDERR ("Key $_ empty contig.\n" );
    next;
  }
  print_Clones_in_Contig( $contigRef );

  print CO ( $contigRef->{num},"\t" );
  print CO ( $contigRef->{name},"\t" );
  if( ! defined $contigRef->{start_bp} ) {
    print CO ( "-1\t" );
  } else {
    print CO ( $contigRef->{start_bp},"\t" );
  }
  print CO ( $contigRef->{length},"\t" );
  print CO ( $contigRef->{chrom},"\n" );
}

sub submit_clone {
  my $chrom = shift;
  my $ordered = shift;
  my $clone = shift;
  
  # filter for double clones
  if( defined $cloneHash{$clone->{embl}."-".$clone->{name}} ) {
    return;
  } else {
    $cloneHash{$clone->{embl}."-".$clone->{name}} = 1;
  }

  if( ! defined $contigHash->{$clone->{contig}} ) {
    $contigNum++;
    $contigHash->{$clone->{contig}}->{num} = $contigNum;
    $contigHash->{$clone->{contig}}->{name} = $clone->{contig};
    
    $contigHash->{$clone->{contig}}->{chrom} = $chrom;
    $contigHash->{$clone->{contig}}->{ordered} = $ordered;
  }

  push( @{$contigHash->{$clone->{contig}}->{clones}}, $clone );

  my $contiglen = $clone->{start} + $clone->{length};
  if( !defined $contigHash->{$clone->{contig}}->{length} ) {
    $contigHash->{$clone->{contig}}->{length} = $contiglen;
  } else {
    if( $contiglen > $contigHash->{$clone->{contig}}->{length} ) {
      $contigHash->{$clone->{contig}}->{length} = $contiglen;
    }
  }

  return;
}

sub print_Clones_in_Contig {
  my $contig = shift;

  for $clone ( @{$contig->{clones}} ) {
    print CL ( $clone->{embl},"\t" );
    print CL ( $clone->{name},"\t" );
    print CL ( $clone->{institute},"\t" );
    print CL ( $contigHash->{$clone->{contig}}->{num},"\t" );
    print CL ( $clone->{length},"\t" );
    print CL ( $clone->{start},"\n" );
  }
}



sub recalc_starts_standard {
  my $chrom = shift;
  my $contigListRef = shift;
  my @contigList = @$contigListRef;
  my $maxlen = $chrlength{"$chrom"};
  print "Chromosome $chrom $maxlen\n";
  my $lastContig = $contigList[ $#contigList ];

  
  my $factor =( $maxlen - $lastContig->{length} ) / $lastContig->{start};
  
  for my $contigRef ( @contigList ) {
    $contigRef->{start_bp} = int( $contigRef->{start}*$factor );
#    print( $contigRef->{start_bp}," ", 
#	   ($contigRef->{start_bp} + $contigRef->{length}), "\n" );
  }
}

sub recalc_eye_candy {
  my $chrom  =shift;
  my $contigListRef = shift;

  my @contigList = @$contigListRef;
  my $maxlen = $chrlength{"$chrom"};

  print "Chromosome $chrom $maxlen\n";

  # length we have
  my $sumLength = 0 ;
  for my $contigRef ( @contigList ) {
    $sumLength += $contigRef->{length};
  }

  my $gap = $maxlen - $sumLength;
  if( $gap < 0 ) {
    # die( "Chromosome length estimate for $chrom is tooo low, please update.\n" );
    print( "Chromosome $chrom is ",-$gap, " too short, no gaps\n" );
    $gap = 0;
  }

  if( $#contigList > 0 ) {
    $allgap = $maxlen - $sumLength;
    $gap = int( $allgap / $#contigList );
    if( $gap < $mingap ) {
      $needgap = $mingap * $#contigList - $allgap;
      
      fake_clonestarts( $contigListRef, $needgap );
      $gap = $mingap;

    }
  }

  my $start = 0;
  for my $contigRef ( @contigList ) {
    $contigRef->{start_bp} = $start;
    $start += ( $gap + $contigRef->{length} );
  }
}
 
# The contigs are too long and cover too much chromosome
# so we make wrong clonestarts within the contigs to make
# everything shorter


sub fake_clonestarts {
  my $contigListRef = shift;
  my $neededGap = shift;

  my @contigList = @$contigListRef;

  my ( $maxStartPositions, $clone, $end );

  for my $contig ( @contigList ) {
    my $contigLength = 0;
    my $endCloneStart = 0;
    for my $clone ( @{$contig->{clones}} ) {
      $end = $clone->{start} + $clone->{length};
      if( $end > $contigLength ) {
	$contigLength = $end;
	$endCloneStart = $clone->{start};
      }
    }
    $maxStartPositions += $endCloneStart;
  }

  my $factor = ($maxStartPositions - $neededGap) / $maxStartPositions;
  print "Set Factor to $factor.\n";

  for my $contig ( @contigList ) {
    my $contigLength = 0;
    for my $clone ( @{$contig->{clones}} ) {
      $clone->{start} = int( $factor * $clone->{start});
      $end = $clone->{start} + $clone->{length};
      if( $end > $contigLength ) {
	$contigLength = $end;
      }
    }
    
    $contig->{length} = $contigLength;
  }
}





