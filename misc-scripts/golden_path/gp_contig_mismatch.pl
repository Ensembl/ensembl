# check innconsistencies between contig offsets and golden path
# when a fragment of the golden path overlaps ends of contigs in 
# a way that its not contained completely by the sequenced contig.

# takes three arguments.first is contig file
# this has rows containing tab delimited
# contig_name \t offset in embl \t length of the contig \t dna_id in ensembl

# second argument is the directory where the .agp files live.
# A find for all files ending in .agp is run over it.

use File::Find;


$contig_file = shift;
$agp_dir =shift;

open( CF, $contig_file ) or die( "Couldnt open contigfile." );
print STDERR ( "Start reading contig.txt.\n" );

while( <CF> ) {
  my ( $contig, $offset, $length  ) = split;
  ( $clone ) = $contig =~ /^([^\.]+\.\d+)\./;
  $clone_offset_hash{$clone}->{$contig} = [ $offset, $length ];
}
close CF;

print STDERR ( "Contig.txt read.\n" );

find( \&append_agp, $agp_dir );
print STDERR ( "Reading ", $#filelist+1, " files.\n" );

$count = 0;

for $filename (@filelist) {
  read_agp( $filename );
  $count++;
  if( $count % 100 == 0 ) {
    print STDERR ( "$count read.\n" );
  }
}

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

    print ( "-----------\n" );
    print clone_hash_string( $clone );;
    print ( "-----------\n\n" );
  }    
}


print qq
(From $readLines read lines, $known_clone lines were in Freeze.
$bogusLines of them were bogus.
$mismatched lines were mismatches.
$unknown_clone lines had unknown clones.
);


sub append_agp {
  /\.agp$/ &&
    push( @filelist, $File::Find::name );
}


sub read_agp {
  $filename = shift;
  local *FH;
  open( FH, $filename ) or die( "Cant open file $filename." );
  
  while( <FH> ) {
    /^#/ && next;
    @list = split( /\t/, $_ );
    if( $list[4] =~ /[^N]/ ) {
      $readLines++;
      ( $clone ) = $list[5] =~ /^([^\.]+\.\d+)$/;
      $start = $list[6];
      $end = $list[7];
      if( defined( $clone_offset_hash{$clone} )) {

	$known_clone++;
	check_line( $clone, $start, $end, $_ );
      } else {
	$unknown_clone++;
	print "Unknown $clone\n";
      }
    }
  }
  close FH;
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


sub check_line {
  my ( $clone, $start, $end, $line ) = @_;
  my $clonehashref = $clone_offset_hash{$clone};
  my ( $offset, $length );

  for my $contig ( keys %{$clonehashref} ) {
    ( $offset, $length ) = @{$clonehashref->{$contig}};
    # match check
    if( $start >= $offset && ($end <= ( $offset+$length-1))) {
      if( defined $contigUsage{$contig} ) {
	push( @{$contigUsage{$contig}}, $line );
	$fragmented{$clone}->{$contig} = 1;
	$double_used_contigs++;
      } else {
	push( @{$contigUsage{$contig}}, $line );
      }
    }
    # partly match == mismatch
    if((( $start < $offset ) && ( $end >= $offset )) ||
       (( $start <= ($offset+$length-1)) && ( $end > ($offset+$length-1)))) {
      push( @{$mismatched{$clone}}, $line );
      $mismatched++;
      return;
    }
  }
}
