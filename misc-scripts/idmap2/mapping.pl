# this script produces idmapping for genes, transcripts, translations and exons
# It needs Database orig with stable_ids set
# it needs Database target where the stable ids are to be set
# it produces log output

# it uses direct sql for speed. All queries are in SQL.pm
# it uses LSF for some jobs I guess ...

use DBI;
use SQL;


my $sdbh = DBI->connect( "DBI:mysql:host=geordy;port=3310;database=arne_ens080", "ensro" );

# my $tdbh = 

print STDERR "Start: ",scalar(localtime()),"\n";
my $oExonInfo = &SQL::orig_exon_information( $sdbh );
# my $nExonInfo = &SQL::target_exon_information( $tdbh );
print STDERR "Finish: ",scalar(localtime()),"\n";

print "Count: ",scalar( @$oExonInfo ),"\n";

direct_mapping( $oExonInfo, undef );



$sdbh->disconnect();
exit;


# extract information (lot of data)
# direct mappings
# contig version update



# similar as direct mapping
# checks if there is an old exon, and a new exon on same clone but
# updated version

sub contig_version_update {
  

}





# takes two exon lists and marks old exons with exact new exons matching
# marks new exons as matched as well
# does the direct approach. Exons are equal, when same position on unchanged contig

sub direct_mapping {
  my ( $old, $new ) = @_;
  
  my ( $sold, $snew );

  print STDERR "Direct Mapping started: ",scalar( localtime() ),"\n";

  @{$sold} = sort { exon_direct_compare( $a, $b ) } @$old;
  @{$snew} = sort { exon_direct_compare( $a, $b ) } @$new;

  print STDERR "Sorted: ",scalar(localtime()),"\n";
  
  my ( $new_iter, $old_iter );
  my ( $old_exon, $new_exon );

  $new_iter = 0;
  $old_iter = 0;

  while( 1 ) {
    $old_exon = $sold->[$old_iter];
    $new_exon = $snew->[$new_iter];
    my  $exon_compare_result;
    if( $old_exon->{'clone_id'} == $new_exon->{'clone_id'} &&
	$old_exon->{'clone_version'} == $new_exon->{'clone_version'} &&
	$old_exon->{'exon_strand'} == $new_exon->{'exon_strand'} &&
	$old_exon->{'contig_offset'} == $new_exon->{'contig_offset'} ) {
      $exon_compare_result = 0;
    } else {
      $exon_compare_result = exon_direct_compare( $old_exon, $new_exon );
    }
    
    if(  $exon_compare_result == 0 ) {
      # direct map possible
      
      # find overlap between the two exons. Need 95% of old exon
      # covered

      my ( $overlap ) = 
      sort { $a <=> $b } ( $old_exon->{'exon_end'}-$old_exon->{'exon_start'}+1,
			   $new_exon->{'exon_end'}-$new_exon->{'exon_start'}+1,
			   $old_exon->{'exon_end'}-$new_exon->{'exon_start'}+1,
			   $new_exon->{'exon_end'}-$old_exon->{'exon_start'}+1 );
      if(( $overlap / 
	  ( $new_exon->{'exon_end'} - $nex_exon->{'exon_start'} )) > 0.95 ) {
	
	$old_exon->{'mapped'} = 'direct';
	$new_exon->{'exon_stable'} = $old_exon->{'exon_stable'};
	$old_iter++;
	$new_iter++;
	$direct_maps++;
	next;
      }
    }

    $exon_compare_result = exon_direct_compare( $old_exon, $new_exon );

    if( $exon_compare_result == -1 ) {
      $old_iter++;
      if( $old_iter > scalar( @$sold ) ) {
	last;
      }
      next;

    } elsif( $exon_compare_result == 1 ) {
      $new_iter++;
      if( $new_iter > scalar( @$snew ) ) {
	last;
      }
      next;

    } else {
      $old_iter++;
      $new_iter++;
      if( $old_iter > scalar( @$old ) ||
	  $new_iter > scalar( @$new ) ) {
	last;
      }
	
    }
  }

  print STDERR "Old Exons: ",scalar( @$old ),"\n";
  print STDERR "New Exons: ",scalar( @$new ),"\n";
  print STDERR "Mapped : ",$direct_maps,"\n";
  print STDERR "Direct Mapping finished: ",scalar( localtime() ),"\n";
  return;
}

sub exon_direct_compare {
  my ( $a, $b ) = @_;
  return $a->{'clone_id'} <=> $b->{'clone_id'} ||
    $a->{'clone_version'} <=> $b->{'clone_version'} ||
      $a->{'contig_offset'} <=> $b->{'contig_offset'} ||
	$a->{'exon_start'} <=> $b->{'exon_start'} ||
	  $a->{'exon_strand'} <=> $b->{'exon_strand'} ;
}


  


sub print_arrayref {
  my $aref = shift;
  
  for( my $i=1; $i <= @$aref; $i++ ) {
    print "------- new record -------\n";
    for my $key ( keys %{$aref->[$i]} ) {
      print $key,"\t",$aref->[$i]{$key},"\n";
    }
  }
}
