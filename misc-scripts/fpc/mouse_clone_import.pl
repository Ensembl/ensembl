use strict;

my ( %contigHash );
my $nextContig = 1;
my ( @sortedClones, @clones );

open( FH, "mouse_clone_pos.txt" ) or die( "No input" );;

while( <FH> ) {
  my @cols = split( "\t", $_ );

  push( @clones, { 'name' => $cols[0],
		   'chrom' => $cols[1],
		   'start' => $cols[2],
		   'end' => $cols[3],
		   'size' => $cols[3]-$cols[2]+1} );
}
@sortedClones = sort { ( $a->{'chrom'} cmp $b->{'chrom'} ) ||
			 ( $a->{'start'} <=> $b->{'start'} ) } @clones;
@clones = ();

my $currentChrom = undef;
my $contig = undef;;
my $lastClone = undef;

open ( CO, ">Fpc_Contig.txt" );
open ( CL, ">Fpc_Clone.txt" );

print "Sorted clones: ",scalar( @sortedClones),"\n";

for my $clone ( @sortedClones ) {
  if( ( ! defined $lastClone ) || 
      ( $lastClone->{'chrom'} ne $clone->{'chrom'} ) ||
      ( $lastClone->{'end'}+1 <= $clone->{'start'} )) {
    
    if( defined $contig ) {
      print CO ( join ( "\t", 
			( $contig->{'no'},
			  "ctg".$contig->{'no'},
			  $contig->{'start'},
			  $contig->{'end'} - $contig->{'start'} +1,
			  $contig->{'chrom'} )),"\n" );
    }		
  
    $contig = { 'start' => $clone->{'start'},
		'end' => $clone->{'end'},
		'chrom' =>$clone->{'chrom'},
		'no' => $nextContig++ };
  } 
  print CL ( join ( "\t",
		    $clone->{'name'},
		    $clone->{'name'},
		    "SC",
		    $contig->{'no'},
		    0,
		    $clone->{'size'},
		    $clone->{'size'},
		    $clone->{'start'}-$contig->{'start'}+1 )),"\n";
  if( $clone->{'end'} > $contig->{'end'} ) {
    $contig->{'end'} = $clone->{'end'};
  }

  $lastClone = $clone;
}

close( CL );
close( CO );
