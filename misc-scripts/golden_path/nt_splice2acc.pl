

use strict;

my $mouse_ctg = shift;
my $mouse_ctg_agp = shift;

my %ctg; # nt to ctg mapping
my %nt;  # ctg to nt mapping
my %nt_agp; # for each NT, an array of hashes, each hash with a line from the agp

if( !defined $mouse_ctg_agp ) {
  die "Must call mouse_ctg mouse_ctg_agp < mouse_agp file";
}


open(CTG,$mouse_ctg) || die "could not open $mouse_ctg $!";


#
# Load up the ctg <-> nt mapping file, ensuring that the mapping is consistent
#

while(<CTG>) {
  my ($nt,$ctg) = split;

  if( defined $ctg{$nt} ) {
    if( $ctg{$nt} ne $ctg ) {
      die "NT contig $nt has two ctgs, $ctg{$nt} and $ctg";
    }
  } else {
    $ctg{$nt} = $ctg;
    $nt{$ctg} = $nt;
  }
}

open(AGP,$mouse_ctg_agp) || die "could not open $mouse_ctg $!";

#
# Load up the NT agp file
#

while( <AGP> ) {
  #Mm1_25680             1     175226   1  F  AC034108.24         1   175226  +
  my ($ctg,$start,$end,$ord,$tag,$acc,$rstart,$rend,$rori) = split;

  if( !defined $nt{$ctg} ) {
    die "For contig $ctg, no NT contig assigned";
  }
  my $nt = $nt{$ctg};

  if( !defined $nt_agp{$nt} ) {
    $nt_agp{$nt} = [];
  }

  my $line = {};
  $line->{'start'} = $start;
  $line->{'end'}   = $end;
  $line->{'ord'}   = $ord;
  $line->{'tag'}   = $tag;
  $line->{'acc'}   = $acc;
  $line->{'rstart'} = $rstart;
  $line->{'rend'}   = $rend;
  $line->{'rori'}   = $rori;

  push(@{$nt_agp{$nt}},$line);
}  



#
# Run down the agp file, keeping track of chromosome splits (order changes) 
# when an NT is mentioned, splice in the correct lines
#

my $current_chr = "";
my $order;

while( <> ) {
  #chr1	1	3000000	1	N	3000000	clone	no
  #chr1	3000001	3002129	2	W	contig_188925	1	2129	+

  my ($chr,$start,$end,$ord,$tag,$acc,$rstart,$rend,$rori) = split;

  if( $chr ne $current_chr ) {
    $current_chr = $chr;
    $order = 1;
  }

  if( $tag eq 'N' || $acc =~ /^contig/ ) {
    print "$chr\t$start\t$end\t$order\t$tag\t$acc\t$rstart\t$rend\t$rori\n";
    $order++;
    next;
  }

  if( $acc !~ /^NT/ ) {
    die "A non NT contig, no contig line called $_\n";
  }



  $acc =~ s/\.\d+$//g;
  
  if( !defined $nt_agp{$acc} ) {
    die "No accession defined for $acc";
  }
  
  #
  # Handle trimming now, as we can do this for both + and -
  # strands the same.
  #

  my @lines;

  foreach my $line ( @{$nt_agp{$acc}} ) {
    # deal with the start trim
    if( $line->{'end'} < $rstart ) {
      # discard this whole line!
      next;
    } 
    if( $line->{'start'} < $rstart ) {
      # need to trim
      $line->{'rstart'} = $rstart - $line->{'start'} +1;
    }

    # end trim
    if( $line->{'start'} > $rend ) {
      # discard this line
      next;
    }

    if( $line->{'end'} > $rend ) {
      #
      $line->{'rend'} = $line->{'rend'} - ($line->{'end'} - $rend);
    }

    push(@lines,$line);
  }


  if( $rori eq '+' ) {
    # easy - splice in - assumming perfect abutting!

    my $rel_len = 0;
    
    
    foreach my $sub ( @lines ) {
     
      my $cstart = $start+$rel_len;
      my $cend   = $start+$rel_len + ($sub->{'rend'} - $sub->{'rstart'}+1)-1;

      $rel_len += ($sub->{'rend'} - $sub->{'rstart'}+1);

      print "$chr\t$cstart\t$cend\t$order\t",$sub->{'tag'},"\t",$sub->{'acc'},"\t",$sub->{'rstart'},"\t",$sub->{'rend'},"\t",$sub->{'rori'},"\n";
      $order++;
    }
  } else {
    my $rel_len = 0;

    foreach my $sub ( reverse @lines ) {
      my $cstart = $start+$rel_len;
      my $cend   = $start+$rel_len + ($sub->{'rend'} - $sub->{'rstart'}+1)-1;

      $rel_len += ($sub->{'rend'} - $sub->{'rstart'}+1);

      
      my $cori;
      if( $sub->{'rori'} eq '+' ) {
	$cori = '-';
      } else {
	$cori = '+';
      }

      print "$chr\t$cstart\t$cend\t$order\t",$sub->{'tag'},"\t",$sub->{'acc'},"\t",$sub->{'rstart'},"\t",$sub->{'rend'},"\t",$cori,"\n";
      $order++;
    }

  }
}
