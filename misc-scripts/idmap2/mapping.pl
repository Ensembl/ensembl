# this script produces idmapping for genes, transcripts, translations and exons
# It needs Database orig with stable_ids set
# it needs Database target where the stable ids are to be set
# it produces log output

# it uses direct sql for speed. All queries are in SQL.pm
# it uses LSF for some jobs I guess ...

use DBI;
use SQL;
use Bio::PrimarySeq;

use strict;
use Bio::EnsEMBL::Pipeline::Runnable::CrossMatch;
use MapGenes;

my $direct_maps = 0;
my $crossmatches = 0;
my $cloneGroups = 0;
my $seqMatches = 0;
my $cmmaps = 0;


my $sdbh = DBI->connect( "DBI:mysql:host=ecs1f;database=arne_ens080", "ensadmin", "ensembl");
my $tdbh = DBI->connect( "DBI:mysql:host=ecs1f;database=arne_ens100", "ensadmin", "ensembl" );

my $starttime = scalar( localtime() );

print STDERR "Start: ",scalar(localtime()),"\n";
my $oExonInfo = &SQL::orig_exon_information( $sdbh );
my $nExonInfo = &SQL::target_exon_information( $tdbh );
print STDERR "Finish: ",scalar(localtime()),"\n";

print STDERR "Count: ",scalar( @$oExonInfo ),"\n";


direct_mapping( $oExonInfo, $nExonInfo );
print "# finished direct mapping\n";

contig_version_update( $oExonInfo, $nExonInfo );
MapGenes::map_genes( $oExonInfo, $nExonInfo );


$sdbh->disconnect();

$tdbh->disconnect();

print STDERR "Start: $starttime\n";
print STDERR "End : ",scalar( localtime()),"\n";

print STDERR "Old Exon Count: ", scalar( @$nExonInfo ), "\n";
print STDERR "New Exon Count: ", scalar( @$nExonInfo ), "\n";
print STDERR "Direct Mappings: ", $direct_maps,"\n";
print STDERR "Crossmatches: $crossmatches\n";
print STDERR "CloneGroups: $cloneGroups\n";
print STDERR "SeqMatches: $seqMatches\n";
print STDERR "Mappings on Crossmatch: $cmmaps\n";
exit;


open( D1, ">dump1.txt" );
open( D2, ">dump2.txt" );


my ( @dumpOld, @dumpNew );

print STDERR "Dumping ...\n";

@dumpOld = sort { ( $a->{'clone_id'} cmp $b->{'clone_id'} ) ||
		( $a->{'exon_start'} <=> $b->{'exon_start'}) } @$oExonInfo;
@dumpNew = sort { ( $a->{'clone_id'} cmp $b->{'clone_id'} ) ||
		( $a->{'exon_start'} <=> $b->{'exon_start'}) } @$nExonInfo;


*STDOUT = *D1;
for my $ex ( @dumpOld ) {
  print_exon( $ex );
}

close D1;

*STDOUT  = *D2;
for my $ex ( @dumpNew ) {
  print_exon( $ex );
}

close D2;


exit;


# extract information (lot of data)
# direct mappings
# contig version update



# similar as direct mapping
# checks if there is an old exon, and a new exon on same clone but
# updated version

sub contig_version_update {
  my ( $old, $new ) = @_;
  
  my ( $pruneOld, $pruneNew, $sortOld, $sortNew  );
  my ( $old_iter, $new_iter );

  @{$pruneOld} = grep { ! exists $_->{'mapped'} } @{$old};
  @{$pruneNew} = grep { ! exists $_->{'mapped'} } @{$new};

  # prune short exons (<10) would be appropriat
  print "Old Exons left to map: ",scalar( @$pruneOld ),"\n";
  print "New Exons left to map: ",scalar( @$pruneNew ),"\n";
 
  @{$sortOld} = sort { $a->{'clone_id'} cmp $b->{'clone_id'} } @$pruneOld;
  @{$sortNew} = sort { $a->{'clone_id'} cmp $b->{'clone_id'} } @$pruneNew;


  $new_iter = 0;
  $old_iter = 0;

  while( 1 ) {

    my $old_exon = $sortOld->[$old_iter];
    my $new_exon = $sortNew->[$new_iter];

    my $cmp = ( $old_exon->{'clone_id'} cmp $new_exon->{'clone_id'} );
    if( $cmp == 0 ) {
      my $clone_id = $old_exon->{'clone_id'};
      my ( $newCloneExons, $oldCloneExons );
      while( defined $old_exon && 
	     $old_exon->{'clone_id'} eq $clone_id ) {
	push( @{$oldCloneExons}, $old_exon );
	$old_iter++;
	if( $old_iter % 10000 == 0 ) {
	  print STDERR "$old_iter exons done.\n";
	}
	$old_exon = $sortOld->[$old_iter];
      }
	  
      while( defined $new_exon && 
	     $new_exon->{'clone_id'} eq $clone_id ) {
	push( @{$newCloneExons}, $new_exon );
	$new_iter++;
	$new_exon = $sortNew->[$new_iter];
      }
      
      clone_map( $oldCloneExons, $newCloneExons );
      # build groups of exons and compare
    } elsif ( $cmp < 0 ) {
      $old_iter++;
      if( $old_iter % 10000 == 0 ) {
	print STDERR "$old_iter exons done.\n";
      }
    } else {
      $new_iter++;
    }
    if( $new_iter >= scalar( @$sortNew ) || 
	$old_iter >= scalar( @$sortOld )) {
      last;
    }
    if( $old_iter % 10000 == 0 ) {
      print STDERR "$old_iter exons done.\n";
    }
  }

}



# map two sets of exons in the same clone from old to new
# by using crossmatch module ...
sub clone_map {
  my ( $oldCloneExons, $newCloneExons ) = @_;
  my ( $oldExHash, $newExHash );
  local *HERR, *HOUT;
 
  open ( NULL, ">/dev/null" );
  
  $cloneGroups++;

  for my $exon ( @$oldCloneExons ) {
    $oldExHash->{$exon->{'exon_id'}} = SQL::exon_sequence( $sdbh, $exon->{'exon_id'} );
    # $exon->{'seq'} = $oldExHash->{$exon->{'exon_id'}};
  }

  for my $exon ( @$newCloneExons ) {
    $newExHash->{$exon->{'exon_id'}} = SQL::exon_sequence( $tdbh, $exon->{'exon_id'} );
    # $exon->{'seq'} = $newExHash->{$exon->{'exon_id'}};
  }

  for my $oEx ( @$oldCloneExons ) {
    for my $nEx ( @$newCloneExons ) {
      if( $oldExHash->{$oEx->{'exon_id'}} eq 
	  $newExHash->{$nEx->{'exon_id'}} ) {
	# we have equal sequence exons here, map them
	$oEx->{'mapped'} = 'same Clone Seq';
	$nEx->{'mapped'} = 'same Clone Seq';
	$nEx->{'exon_stable'} = $oEx->{'exon_stable'};
	print ( "Exon mapped:\t", $oEx->{'exon_id'},"\t",
		$nEx->{'exon_id'},"\t", $oEx->{'exon_stable'},"\n" );
	$seqMatches++;
#	print STDERR "---- DIRECT SEQUENCE MAP ----\n";
      }
    }
  }

  # no crossmatch mappings
  return;

  my ( @oldPrune, @newPrune, @scorelist );
  @oldPrune = grep { ! exists $_->{'mapped'} } @$oldCloneExons;
  @newPrune = grep { ! exists $_->{'mapped'} } @$newCloneExons;

  for my $oEx ( @oldPrune ) {
    for my $nEx ( @newPrune ) {
      my $crossMatcher = new Bio::EnsEMBL::Pipeline::Runnable::CrossMatch
	( -seq1 => new Bio::PrimarySeq( -seq => $oldExHash->{$oEx->{'exon_id'}}, 
					-id => $oEx->{'exon_id'} ),
	  -seq2 => new Bio::PrimarySeq( -seq => $newExHash->{$nEx->{'exon_id'}}, 
					-id => $nEx->{'exon_id'} ),
	  -score => 200
	);

      *HERR = *STDERR;
      *HOUT = *STDOUT;

      *STDOUT = *NULL;
      *STDERR = *NULL;

      $crossMatcher->run();

      *STDOUT = *HOUT;
      *STDERR = *HERR;

      $crossmatches++;
      my $score = 0;
      my $matches = 0;
      for my $fp ( $crossMatcher->output() ) {
	$score += $fp->score();
	$matches = 1;
      }
      if( $matches ) {
	push( @scorelist, { score => $score, 
			    oldEx => $oEx,
			    newEx => $nEx } );
      }
									 
    }
  }

  my @sortedScores;
  @sortedScores = sort { $b->{'score'} <=> $a->{'score'} } @scorelist;
  
  for my $scoreRecord ( @sortedScores ) {
    if( ! exists $scoreRecord->{'oldEx'}{'mapped'} &&
	! exists $scoreRecord->{'newEx'}{'mapped'} ) {
      $scoreRecord->{'oldEx'}{'mapped'} = 'crossmatch Clone Seq';
      $scoreRecord->{'newEx'}{'mapped'} = 'crossmatch Clone Seq';;
      $scoreRecord->{'newEx'}{'exon_stable'} = 
	$scoreRecord->{'oldEx'}{'exon_stable'};
      print ( "Exon mapped:\t", $scoreRecord->{'oldEx'}{'exon_id'},"\t",
	      $scoreRecord->{'newEx'}{'exon_id'},"\t", 
	      $scoreRecord->{'oldEx'}{'exon_stable'},"\n" ); 
      
      $cmmaps++;
#      print STDERR "---- CROSSMATCH MAP ----",$scoreRecord->{'score'},"\n";
    }
  }
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
    if(( $old_exon->{'clone_id'} cmp $new_exon->{'clone_id'}) == 0 &&
	$old_exon->{'clone_version'} == $new_exon->{'clone_version'} &&
	$old_exon->{'exon_strand'} == $new_exon->{'exon_strand'} &&
	$old_exon->{'contig_offset'} == $new_exon->{'contig_offset'} ) {
      $exon_compare_result = 0;
    } else {
      $exon_compare_result = exon_direct_compare( $old_exon, $new_exon );
    }
    
    if(  $exon_compare_result == 0 ) {
      # direct map possible
      
      # find overlap between the two exons. Need >50% of new exon
      # covered

      my ( $overlap ) = 
      sort { $a <=> $b } ( $old_exon->{'exon_end'}-$old_exon->{'exon_start'}+1,
			   $new_exon->{'exon_end'}-$new_exon->{'exon_start'}+1,
			   $old_exon->{'exon_end'}-$new_exon->{'exon_start'}+1,
			   $new_exon->{'exon_end'}-$old_exon->{'exon_start'}+1 );
#      print "---- MAPPING ----\n";

#      print_exon( $old_exon );
#      print_exon( $new_exon );
#      print "OVERLAP: $overlap\n";

      if(( $overlap / 
	  ( $new_exon->{'exon_end'} - $new_exon->{'exon_start'} + 1 )) > 0.5 ) {
	
#	print "MAPPED\n";
	$old_exon->{'mapped'} = 'direct';
	$new_exon->{'mapped'} = 'direct';

	$new_exon->{'exon_stable'} = $old_exon->{'exon_stable'};
	print ( "Exon mapped:\t", $old_exon->{'exon_id'},"\t",
		$new_exon->{'exon_id'},"\t", $old_exon->{'exon_stable'},"\n" );
	$old_iter++;
	$new_iter++;

	if( $old_iter >= scalar( @$old ) ||
	    $new_iter >= scalar( @$new ) ) {
	  last;
	}
	$direct_maps++;
	next;
      }
    }
#    print "NOT MAPPED\n";

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

  print STDERR "Direct Mapping finished: ",scalar( localtime() ),"\n";
  return;
}

sub exon_direct_compare {
  my ( $a, $b ) = @_;
  return $a->{'clone_id'} cmp $b->{'clone_id'} ||
    $a->{'clone_version'} <=> $b->{'clone_version'} ||
      $a->{'contig_offset'} <=> $b->{'contig_offset'} ||
	$a->{'exon_start'} <=> $b->{'exon_start'} ||
	  $a->{'exon_strand'} <=> $b->{'exon_strand'} ;
}


  


sub print_exon {
  my $exon = shift;
  for my $key ( sort keys %{$exon} ) {
    print $key,"\t",$exon->{$key},"\n";
  }
  print "\n";
}
