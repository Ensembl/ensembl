# here we expect a list of exons from old and new database
# each exon contains exon_stable to mark it as mapped 

# each exon contains the internal gene_id of the gene it is 
# part of. This module tries to map the gene_stable ids over to
# the new exons.


package MapGenes;
use strict;


# starting with one exon, build_exon_group tries to build an 
# 'envelope' around it. That is a set of genes and exons where neither
# of them misses parts. You call build_exon_group as often as needed
# until it doesnt produce new exons any more

my ( %newGeneHash, %oldGeneHash, 
     %oldExonHash, %newExonHash,
     @oldPruneExon, @newPruneExon );
my $mappedGene;

# go through old Exons
# skip if gene is mapped
# build envelope if not
sub map_genes {
  my ( $oldExonInfo, $newExonInfo )  = @_;

  my %oldGenesMapped; 
  print STDERR ( "Start gene mapping ",scalar( localtime() ),"\n");
  init_lookup_tables( $oldExonInfo, $newExonInfo );

  for my $oExon ( @oldPruneExon ) {
    if( ! exists $oExon->{'gene_mapped'} ) {
      my @map_exons = build_envelope( $oExon, \@oldPruneExon, \@newPruneExon );
      map_envelope( \@map_exons );
    }
  }
  print STDERR ( "End gene mapping ",scalar( localtime() ),"\n");
  print STDERR ( "Mapped Genes ",$mappedGene,"\n" );
}



sub map_envelope {
  my $exon_stables = shift;
  
  my %gene_map_hash;
  my $common_exons = 0;

  for my $exon_stable_id ( @$exon_stables ) {
    if( exists $oldExonHash{ $exon_stable_id } &&
	exists $newExonHash{ $exon_stable_id } ) {
      my ( $oEx, $nEx, $hashkey );
      $oEx = $oldExonHash{ $exon_stable_id };
      $nEx = $newExonHash{ $exon_stable_id };
      $common_exons++;

      $hashkey = $oEx->{'gene_id'}."-".$nEx->{'gene_id'};
      if( exists $gene_map_hash{$hashkey} ) {
	$gene_map_hash{$hashkey}->[2] += 
	  $nEx->{'exon_end'}-$nEx->{'exon_start'}+1;
      } else {
	$gene_map_hash{$hashkey} = [ $oEx->{'gene_id'}, $nEx->{'gene_id'},
				     $nEx->{'exon_end'}-$nEx->{'exon_start'}+1,
				     $oEx->{'gene_stable'}];
							       
      }
    }
  }

  my @sortedMappings = sort { $b->[2] <=> $a->[2] } values %gene_map_hash;
  my ( %oldGenesMapped, %newGenesMapped );

  for my $mapRecord ( @sortedMappings ) {
    if( ! exists $oldGenesMapped{$mapRecord->[0]} &&
	! exists $newGenesMapped{$mapRecord->[1]} ) {
      $oldGenesMapped{$mapRecord->[0]} = $mapRecord->[1];
      $newGenesMapped{$mapRecord->[1]} = $mapRecord->[0];
      print ( "Gene mapped:\t",$mapRecord->[0], "\t",$mapRecord->[1],
	      "\t",$mapRecord->[3],".\n" );
      $mappedGene++;
    }
  }

  for my $exon_stable_id ( @$exon_stables ) {
    if( exists $oldExonHash{ $exon_stable_id }) {
      $oldExonHash{ $exon_stable_id }{'gene_mapped'} = 1;
    }
    if( exists $newExonHash{ $exon_stable_id }) {
      $newExonHash{ $exon_stable_id }{'gene_mapped'} = 1;
    }
  }
  
}


sub init_lookup_tables {
  my ( $oldExonInfo, $newExonInfo ) = @_;
  
  # init the lookup tables
  for my $exon ( @$newExonInfo ) {
    if( exists $exon->{'exon_stable'} ) {
      push( @{$newGeneHash{$exon->{'gene_id'}}}, $exon );
      $newExonHash{ $exon->{'exon_stable'} } = $exon;
      push( @newPruneExon, $exon );
    }
  }

  for my $exon ( @$oldExonInfo ) {
    # should exist for all of them ...
    if( exists $newExonHash{ $exon->{'exon_stable'}} ) {
      push( @{$oldGeneHash{$exon->{'gene_id'}}}, $exon );
      $oldExonHash{ $exon->{'exon_stable'} } = $exon;
      push( @oldPruneExon, $exon );
      my $newExon = $newExonHash{ $exon->{'exon_stable'}};
      # print "vvvv Mapping vvvv\n";
      # ::print_exon( $exon );
      # ::print_exon( $newExon );
      # print "^^^^ Mapping ^^^^\n";
    } 
  }
  

#  print "Lookup initialized: ",scalar( localtime() ),"\n";
#  print ( "old exons:",scalar( @oldPruneExon ),"\n" );
#  print ( "new exons:",scalar( @newPruneExon ),"\n" );
}


# given one old exon, build envelope around it
# return set of stable ids which cover the 
# given exon.

# That means: All old and new genes of these exons have all
# their stable_id exons in this set as well.
sub build_envelope {
  my ( $exon, $oExonInfo, $nExonInfo ) = @_;

  my $oGeneGroup = {};
  my $nGeneGroup = {};
  my $exonSet = {};
  my $newExons = [];
  push( @$newExons, $exon->{'exon_stable'} );

  $newExons = build_exon_group( $oGeneGroup, $exonSet, $newExons,
				\%oldGeneHash, \%oldExonHash );

  while( scalar( @$newExons )) {
    $newExons = build_exon_group( $nGeneGroup, $exonSet, $newExons,
				  \%newGeneHash, \%newExonHash );
    if( scalar( @$newExons )) {
      $newExons = build_exon_group( $oGeneGroup, $exonSet, $newExons,
				    \%oldGeneHash, \%oldExonHash );
    }
  }
  
  # print ( "Exon envelope: ", join( " ", keys %$exonSet),"\n" );
  # print ( "Old Genes:", join( " ", keys %$oGeneGroup ) ,"\n" );
  # print ( "New Genes:", join( " ", keys %$nGeneGroup ) ,"\n\n" );
  return keys %$exonSet;
}


# expand newExonSet that it covers all exons of its genes
# wiuth respect to one exon/gene set.
sub build_exon_group {
  my ( $geneGroup, $exonSet, $newExonSet, 
       $geneLookup, $exonLookup ) = @_;
  # print STDERR ( "Build group ", join( " ", @$newExonSet ),"\n" );
  my (%newGenes, @newExons, $exon);

  # get Genes from new Exons
  # get Exons from this genes
  
  for my $exon_stable ( @$newExonSet ) {
    $exon = $exonLookup->{$exon_stable};
    my $geneId;

    $geneId = $exon->{'gene_id'};
    if( ! defined $geneId ) {
      print STDERR "undefined ";
      &::print_exon( $exon );
    }
    if( !exists $geneGroup->{$geneId} ) {
      $newGenes{$geneId} = 1;
      $geneGroup->{$geneId} = 1;
    }
  }
      

  for my $geneId ( keys %newGenes ) {
    my $exons = $geneLookup->{$geneId};
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
