# Script to parse the golden path files ".agp"
# into ensembl overlap table

# takes three arguments.first is contig file
# this has rows containing tab delimited
# contig_name \t offset in embl \t length of the contig \t dna_id in ensembl

# second argument is the directory where the .agp files live.
# A find for all files ending in .agp is run over it.

# third argument is a file with invalid clonenames
# one per line

use File::Find;


$contig_file = shift;
$agp_dir =shift;
$bogus_clones_file = shift;

open( CF, $contig_file ) or die( "Couldnt open contigfile." );
print STDERR ( "Start reading contig.txt.\n" );

while( <CF> ) {
  my ( $contig, $offset, $length, $dna_id ) = split;
  ( $clone ) = $contig =~ /^([^\.]+)\./;
  $clone_offset_hash{$clone}->{$contig} = [ $offset, $length, $dna_id ];
}
close CF;

print STDERR ( "Contig.txt read.\n" );

open( BC, $bogus_clones_file ) or die( "Couldnt open $bogus_clones_file" );
print STDERR ( "Reading bogus clones\n" );
while( <BC> ) {
  chomp;
  $bogus_clones{$_} = 1;
}
close( BC );


sub append_agp {
  /\.agp$/ &&
    push( @filelist, $File::Find::name );
}

sub read_agp {
  my $filename = shift;
  local *FH;
  open( FH, $filename ) or die( "Cant open file $filename." );
  my ( $start, $end );
  my $firstFrag = 1;
  my $lastFrag = undef;
  my $thisFrag = undef;
  print ( "$filename\n" );
  while( <FH> ) {
    /^#/ && next;
    chomp;
    @list = split( /\t/, $_ );
    if( $list[4] =~ /[^N]/ ) {
      $readLines++;
      ( $clone ) = $list[5] =~ /^([^\.]+)\./;
      $start = $list[6];
      $end = $list[7];
      ( $contig, $dna_id, $contig_start,
	$contig_end ) = check_line( $clone, $start, $end, $_ );
      if( !defined $contig ) {
	# this fragment is invalid for some reason
	# skip it
	# print( "Bang\n" );
	$firstFrag = 1;
	$lastFrag = {};
	$thisFrag = {};
	next;
      }
      if( $list[8] eq '-' ) {
	( $contig_start, $contig_end ) = ( $contig_end, $contig_start );
      }

      $thisFrag->{contig_name} = $contig;
      $thisFrag->{dna_id}  = $dna_id;
      $thisFrag->{start} = $contig_start;
      $thisFrag->{end} = $contig_end;
      $thisFrag->{strand} = $list[8];
      $thisFrag->{gap} = 1;

      if( !$firstFrag ) {
	build_overlap( $lastFrag, $thisFrag );
      }

      $lastFrag = $thisFrag;
      $thisFrag = {};
      if( $firstFrag )  {
	$fragmentCount++;
      }
      $firstFrag = 0;
    } else {
      # gap
      # if( $firstFrag ) {
      # print STDERR ("Contig ", $list[0], " started with gap\n" );
      # close( FH );
      # return;
      # }
      $lastFrag->{gap} = $list[5];
    }
    $start = 0;
    
  }
  close FH;
}

sub check_line {
  my ( $clone, $start, $end, $line ) = @_;
  my $clonehashref = $clone_offset_hash{$clone};
  my ( $offset, $length, $dna_id );

  if( $bogus_clones{$clone} ) {
    return( undef, undef, undef, undef );
  }

  for my $contig ( keys %{$clonehashref} ) {
    ( $offset, $length, $dna_id ) = @{$clonehashref->{$contig}};
    # match check
    if( $start >= $offset && ($end <= ( $offset+$length-1))) {
      return ($contig, $dna_id, $start - $offset + 1, $end - $offset + 1 );
    }
  }
  return ( undef ,undef, undef, undef );
}

sub build_overlap {
  my ( $lastFrag, $thisFrag) = @_;
  print( $lastFrag->{dna_id}, "\t", $thisFrag->{dna_id}, "\t",
	 $lastFrag->{end}, "\t", $thisFrag->{start},"\t" );
  print( $lastFrag->{gap},"\t","ucsc\t" );

  my $type = $lastFrag->{strand}.$thisFrag->{strand};

  if( $type =~ /\+\+/ ) {
    print( "right2left\n" );
  } elsif( $type =~ /\+-/ ) {
    print( "right2right\n" );
  } elsif( $type =~ /-\+/ ) {
    print( "left2left\n" );
  } elsif( $type =~ /--/ ) {
    print( "left2right\n" );
  } else {
    print( "?$type?\n" );
  }
}


sub clone_hash_string {
  my $clonename = shift;
  my $hashref = $clone_offset_hash{$clonename};
  my $result;
  my @keylist_sorted;

  $result .= "EnsEMBL Freeze Clone $clonename\n";
  @keylist_sorted = sort { $hashref->{$a}[0] <=> $hashref->{$b}[0] } ( keys %{$hashref} );
  for my $contig ( @keylist_sorted ) {
    $result .= "$contig ";
    my ( $offset, $length );
    ( $offset, $length ) = @{$hashref->{$contig}};
    $result .= "$offset ".($offset+$length-1)."\n";
  }

  return $result;
}



find( \&append_agp, $agp_dir );
print STDERR ( "Reading ", $#filelist+1, " files.\n" );

$count = 0;

for $filename (@filelist) {
  read_agp( $filename );
  $count++;
  if( $count % 100 == 0 ) {
  #  last;
    print STDERR ( "$count read.\n" );
  }
}

print STDERR ("Counted $fragmentCount fragments.\n" );

exit;
# fragmented
print( "The folllowing lines fragmented or mismatched contigs:\n\n" );
for my $clone ( keys %fragmented ) {
  for my $contig ( keys %{$fragmented{$clone}} ) {
    for my $line ( @{$contigUsage{$contig}} ) {
      $bogusLines++;
      print $line;
    }
  }

  if( defined $mismatched{$clone} ) {
    for my $line ( @{$mismatched{$clone}} ) {
      print $line;
      $bogusLines++;
    }
  }

  print ( "-----------\n" );
  print clone_hash_string( $clone );;
  print ( "-----------\n\n" );
}

for my $clone ( keys %mismatched ) {
  if( ! defined $fragmented{$clone} ) {
    for my $line ( @{$mismatched{$clone}} ) {
      print $line;
      $bogusLines++;
    }
  }    

  print ( "-----------\n" );
  print clone_hash_string( $clone );;
  print ( "-----------\n\n" );
}


print qq
(From $readLines read lines, $known_clone lines were in Freeze.
$bogusLines of them were bogus.
$mismatched lines were mismatches.
$unknown_clone lines had unknown clones.
);

exit;

