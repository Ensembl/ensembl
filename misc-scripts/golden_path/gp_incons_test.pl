# find inconsistencies in .agp files.
# report when two fragments in agp overlap the same embl sequence.

# argument: the agp subdirectory. A find is done over it.


use File::Find;

$agp_dir =shift;

sub append_agp {
  /\.agp$/ &&
    push( @filelist, $File::Find::name );
}

sub read_agp {
  $filename = shift;
  local *FH;
  open( FH, $filename ) or die( "Cant open file $filename." );
  my %clonestore;
  my $clone;

  while( <FH> ) {
    /^#/ && next;
    @list = split( /\t/, $_ );
    if( $list[4] =~ /[^N]/ ) {
      $readLines++;
      ( $clone ) = $list[5] =~ /^([^\.]+)\./;
      $start = $list[6];
      $end = $list[7];
      if( $list[1] > $list[2] ) {
	print $_;
	$reverse++;
      }
      next;
      push( @{$clonestore{$clone}}, [ $start, $end, $_] );
    }
  }
  close FH;
  return;
  my ( $s1, $s2, $e1, $e2, $l1, $l2 );

  # check overlaps in %clonestore
  for my $clone ( keys %clonestore ) {
    my @frags = @{$clonestore{$clone}};

    for( my $i = 0; $i <= $#frags; $i++ ) {
      for( my $j = $i+1; $j <= $#frags; $j++ ) {
	( $s1, $e1, $l1 ) = @{$frags[$i]};
	( $s2, $e2, $l2 ) = @{$frags[$j]};
	
	if(( $s2 >= $s1 && $s2 <=$e1 ) ||
	   ( $s1 >= $s2 && $s1 <=$e2 )) {
	  print $l1,$l2,"\n";
	  $collisionCount++
	}
      }
    }
  }
}


find( \&append_agp, $agp_dir );
print STDERR ( "Reading ", $#filelist+1, " files.\n" );

$count = 0;

for $filename (@filelist) {
  read_agp( $filename );
  $count++;
  if( $count % 100 == 0 ) {
    # last;
    print STDERR ( "$count read.\n" );
  }
}

print( "Reverse count $reverse\n" );
print( "\nFound $collisionCount contradicting line pairs\n" );
print( "in $readLines read lines.\n" );

exit;

