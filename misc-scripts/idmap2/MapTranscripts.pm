# here we expect a list of exons from old and new database
# each exon contains exon_stable to mark it as mapped 

# each exon contains the internal transcript_id of the transcript it is 
# part of. This module tries to map the transcript_stable ids over to
# the new exons.


package MapTranscripts;
use strict;


# starting with one exon, build_exon_group tries to build an 
# 'envelope' around it. That is a set of transcripts and exons where neither
# of them misses parts. You call build_exon_group as often as needed
# until it doesnt produce new exons any more

my ( %newTranscriptHash, %oldTranscriptHash, 
     %oldExonHash, %newExonHash );

my $mappedTranscripts;


# go through old Exons
# skip if Transcript is mapped
# build envelope if not
sub map_transcripts {
  my ( $oldExonInfo, $newExonInfo )  = @_;

  my %oldGenesMapped; 
  print STDERR ( "Start transcript mapping ",scalar( localtime() ),"\n");
  init_lookup_tables( $oldExonInfo, $newExonInfo );

  for my $oExon ( @$oldExonInfo ) {
    if( ! exists $oExon->{'transcript_mapped'} ) {
      my @map_exons = build_envelope( $oExon, $oldExonInfo, $newExonInfo );
      map_envelope( \@map_exons );
    }
  }
  print STDERR ( "End transcript mapping ",scalar( localtime() ),"\n");
  print STDERR ( "Mapped Transcripts ",$mappedTranscripts,"\n" );
}



sub map_envelope {
  my $exon_stables = shift;
  
  my %transcript_map_hash;
  my $common_exons = 0;
  my ( @oldExList, @newExList );

  for my $exon_stable_id ( @$exon_stables ) {
    if( exists $oldExonHash{ $exon_stable_id } &&
	exists $newExonHash{ $exon_stable_id } ) {
      my ( $oEx, $nEx, $hashkey );
      @oldExList = @{$oldExonHash{ $exon_stable_id }};
      @newExList = @{$newExonHash{ $exon_stable_id }};
      for my $oEx ( @oldExList ) {
	for my $nEx ( @newExList ) {

	  $hashkey = $oEx->{'transcript_id'}."-".
	    $nEx->{'transcript_id'};
	  if( exists $transcript_map_hash{$hashkey} ) {
	    $transcript_map_hash{$hashkey}->[2] += 
	      $nEx->{'exon_end'}-$nEx->{'exon_start'}+1;
	  } else {
	    $transcript_map_hash{$hashkey} = 
	      [ $oEx->{'transcript_id'}, $nEx->{'transcript_id'},
		$nEx->{'exon_end'}-$nEx->{'exon_start'}+1,
		$oEx->{'transcript_stable'}];
	    
	  }
	}
      }
    }
  }

  my @sortedMappings = sort { $b->[2] <=> $a->[2] } values %transcript_map_hash;
  my ( %oldTranscriptsMapped, %newTranscriptsMapped );

  for my $mapRecord ( @sortedMappings ) {
    if( ! exists $oldTranscriptsMapped{$mapRecord->[0]} &&
	! exists $newTranscriptsMapped{$mapRecord->[1]} ) {
      $oldTranscriptsMapped{$mapRecord->[0]} = $mapRecord->[1];
      $newTranscriptsMapped{$mapRecord->[1]} = $mapRecord->[0];
      print ( "Transcript mapped:\t",$mapRecord->[0], "\t",$mapRecord->[1],
	      "\t",$mapRecord->[3],".\n" );
      $mappedTranscripts++;
    }
  }
  
  for my $exon_stable_id ( @$exon_stables ) {
    if( exists $oldExonHash{ $exon_stable_id }) {
      for my $exon ( @{$oldExonHash{ $exon_stable_id }} ) {
	$exon->{'transcript_mapped'} = 1;
      }
    }
    if( exists $newExonHash{ $exon_stable_id }) {
      for my $exon ( @{$newExonHash{ $exon_stable_id }} ) {
	$exon->{'transcript_mapped'} = 1;
      }
    }
  }

}


sub init_lookup_tables {
  my ( $oldExonInfo, $newExonInfo ) = @_;

  # problem: not all exons which are mapped ar emarked as such
  # that happened if they are in more than one transcript
  my %exHash= ();


  for my $exon ( @$newExonInfo ) {
    if( exists $exon->{'mapped'} ) {
      $exHash{$exon->{'exon_id'}." ".$exon->{'sticky_rank'}." mapped"} =
	$exon->{'mapped'};
      $exHash{$exon->{'exon_id'}." ".$exon->{'sticky_rank'}." stable"} =
	$exon->{'exon_stable'};
    }
  }

  for my $exon ( @$newExonInfo ) {
    if( ! exists $exon->{'mapped'} && 
	( exists $exHash{$exon->{'exon_id'}." ".$exon->{'sticky_rank'}." mapped"} )) {
      $exon->{'mapped'}= $exHash{$exon->{'exon_id'}." ".$exon->{'sticky_rank'}." mapped"};
      $exon->{'exon_stable'}= $exHash{$exon->{'exon_id'}." ".$exon->{'sticky_rank'}." stable"};
    }
    push( @{$newTranscriptHash{$exon->{'transcript_id'}}}, $exon );
    push( @{$newExonHash{ $exon->{'exon_stable'} }}, $exon );
  }

  %exHash = ();

  for my $exon ( @$oldExonInfo ) {
    if( exists $exon->{'mapped'} ) {
      $exHash{$exon->{'exon_id'}." ".$exon->{'sticky_rank'}} =
	$exon->{'mapped'};
    }
  }

  for my $exon ( @$oldExonInfo ) {
    if( ! exists $exon->{'mapped'} && 
	( exists $exHash{$exon->{'exon_id'}." ".$exon->{'sticky_rank'}} )) {
      $exon->{'mapped'}= $exHash{$exon->{'exon_id'}." ".$exon->{'sticky_rank'}};
    }
    if( exists $exon->{'mapped'} ) {
      push( @{$oldTranscriptHash{$exon->{'transcript_id'}}}, $exon );
      push( @{$oldExonHash{ $exon->{'exon_stable'} }}, $exon );
    }
  }
  # remove all exons from $oldExonInfo which dont map
  my $iter = 0;
  while( $iter < scalar @$oldExonInfo ) {
    my $exon = $oldExonInfo->[$iter];
    if( ! exists $exon->{'mapped'} ) {
      splice( @$oldExonInfo, $iter ,1 );
    } else {
      $iter++;
    }
  }
}


# given one old exon, build envelope around it
# return set of stable ids which cover the 
# given exon.

# That means: All old and new genes of these exons have all
# their stable_id exons in this set as well.
sub build_envelope {
  my ( $exon, $oExonInfo, $nExonInfo ) = @_;

  my $oTranscriptGroup = {};
  my $nTranscriptGroup = {};
  my $exonSet = {};
  my $newExons = [];

  push( @$newExons, $exon->{'exon_stable'} );

  $newExons = build_exon_group( $oTranscriptGroup, $exonSet, $newExons,
				\%oldTranscriptHash, \%oldExonHash );

  while( scalar( @$newExons )) {
    $newExons = build_exon_group( $nTranscriptGroup, $exonSet, $newExons,
				  \%newTranscriptHash, \%newExonHash );
    if( scalar( @$newExons )) {
      $newExons = build_exon_group( $oTranscriptGroup, $exonSet, $newExons,
				    \%oldTranscriptHash, \%oldExonHash );
    }
  }
  
#  print ( "Exon envelope: ", join( " ", keys %$exonSet),"\n" );
#  print ( "Old Transcripts:", join( " ", keys %$oTranscriptGroup ) ,"\n" );
#  print ( "New Transcripts:", join( " ", keys %$nTranscriptGroup ) ,"\n\n" );
  return keys %$exonSet;
}


# expand newExonSet that it covers all exons of its transcripts
# with respect to one exon/transcript set.
sub build_exon_group {
  my ( $transcriptGroup, $exonSet, $newExonSet, 
       $transcriptLookup, $exonLookup ) = @_;
  # print STDERR ( "Build group ", join( " ", @$newExonSet ),"\n" );
  my (%newTranscripts, @newExons, $exon);

  # get Transcrips from new Exons
  # get Exons from this Transcripts
  
  for my $exon_stable ( @$newExonSet ) {
    for my $exon ( @{$exonLookup->{$exon_stable}} ) {
      my $transcriptId;

      $transcriptId = $exon->{'transcript_id'};
      if( ! defined $transcriptId ) {
	print STDERR "undefined ";
	&::print_exon( $exon );
      }
      if( !exists $transcriptGroup->{$transcriptId} ) {
	$newTranscripts{$transcriptId} = 1;
	$transcriptGroup->{$transcriptId} = 1;
      }
    }
  }
      

  for my $transcriptId ( keys %newTranscripts ) {
    my $exons = $transcriptLookup->{$transcriptId};
    for my $exon ( @$exons ) {
      if( ! exists( $exonSet->{$exon->{'exon_stable'}} ) ) {
	$exonSet->{$exon->{'exon_stable'}} = 1;
	push( @newExons, $exon->{'exon_stable'} );
#	print STDERR "Pushed Exon to Envelope\n";
      }
    }
  }
  # print STDERR ( "Exon envelope: ", join( " ", @newExons),"\n\n" );

  return \@newExons;
}


1;
