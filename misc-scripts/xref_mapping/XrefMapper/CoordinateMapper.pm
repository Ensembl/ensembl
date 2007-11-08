# $Id$

# This is a set of subroutines used for creating Xrefs based on
# coordinate overlaps.

package XrefMapper::CoordinateMapper;

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Mapper::RangeRegistry;

use Carp;
use IO::File;
use File::Spec::Functions;

use base qw( Exporter );

our @EXPORT = qw( run_coordinatemapping );

our $coding_weight = 2;
our $ens_weight    = 3;

our $transcript_score_threshold = 0.95;
our $gene_score_threshold       = 0.95;

sub run_coordinatemapping {
  my ($mapper, $do_upload) = @_;

  my $xref_db = $mapper->xref();
  my $core_db = $mapper->core();

  my $species = $core_db->species();
  my $species_id =
    XrefMapper::BasicMapper::get_species_id_from_species_name( $xref_db,
                                                             $species );

  # We only do coordinate mapping for mouse and human for now.
  if ( !( $species eq 'mus_musculus' || $species eq 'homo_sapiens' ) ) {
    return;
  }

  my $output_dir = $core_db->dir();

  my $unmapped_reason_filename =
    catfile( $output_dir, 'unmapped_reason_coord.txt' );
  my $unmapped_obj_filename =
    catfile( $output_dir, 'unmapped_object_coord.txt' );

  my $core_db_adaptor =
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                   -host => $core_db->dbc()->host(),
                                   -port => $core_db->dbc()->port(),
                                   -user => $core_db->dbc()->username(),
                                   -pass => $core_db->dbc()->password(),
                                   -dbname => $core_db->dbc()->dbname(),
    );

  my $slice_adaptor = $core_db_adaptor->get_SliceAdaptor();
  my @chromosomes   = @{ $slice_adaptor->fetch_all('Chromosome') };

  my $xref_dbh = $xref_db->dbc()->db_handle();
  my $core_dbh = $core_db->dbc()->db_handle();

  ######################################################################
  # Figure out the last used 'xref_id', 'object_xref_id',              #
  # 'unmapped_object_id', and 'unmapped_reason_id' from the Core       #
  # database.                                                          #
  ######################################################################

  my $xref_id =
    $core_dbh->selectall_arrayref('SELECT MAX(xref_id) FROM xref')
    ->[0][0];
  my $object_xref_id = $core_dbh->selectall_arrayref(
                 'SELECT MAX(object_xref_id) FROM object_xref')->[0][0];
  my $unmapped_object_id = $core_dbh->selectall_arrayref(
         'SELECT MAX(unmapped_object_id) FROM unmapped_object')->[0][0];
  my $unmapped_reason_id = $core_dbh->selectall_arrayref(
         'SELECT MAX(unmapped_reason_id) FROM unmapped_reason')->[0][0];

  log_progress( "Last used xref_id            is %d\n", $xref_id );
  log_progress( "Last used object_xref_id     is %d\n",
                $object_xref_id );
  log_progress( "Last used unmapped_object_id is %d\n",
                $unmapped_object_id );
  log_progress( "Last used unmapped_reason_id is %d\n",
                $unmapped_reason_id );

  ######################################################################
  # Get an 'analysis_id', or discover that we need to add our analysis #
  # to the 'analyis' table later.                                      #
  ######################################################################

  my $analysis_params =
    sprintf( "weights(coding,ensembl)="
               . "%.2f,%.2f;"
               . "thresholds(transcript,gene)="
               . "%.2f,%.2f",
             $coding_weight, $ens_weight, $transcript_score_threshold,
             $gene_score_threshold );

  my $analysis_sql = qq(
    SELECT  analysis_id
    FROM    analysis
    WHERE   logic_name = 'XrefCoordinateMapping'
    AND     parameters = ?
  );

  my $analysis_sth = $core_dbh->prepare($analysis_sql);
  $analysis_sth->execute($analysis_params);

  my $analysis_id = $analysis_sth->fetchall_arrayref()->[0][0];
  if ( !defined($analysis_id) ) {
    log_progress("Can not find analysis ID for this analysis:\n");
    log_progress("  logic_name = 'XrefCoordinateMapping'\n");
    log_progress( "  parameters = '%s'\n", $analysis_params );
    log_progress("A new analysis will be added\n");
  } else {
    log_progress( "Analysis ID                  is %d\n",
                  $analysis_id );
  }

  ######################################################################
  # Read and store available Xrefs from the Xref database.             #
  ######################################################################

  my %unmapped;
  my %mapped;

  my $xref_sql = qq(
    SELECT  coord_xref_id, source_id, accession
    FROM    coordinate_xref
    WHERE   species_id = ?
  );

  my $xref_sth = $xref_dbh->prepare($xref_sql);
  $xref_sth->execute($species_id);

  while ( my $xref = $xref_sth->fetchrow_hashref() ) {
    $unmapped{ $xref->{'coord_xref_id'} } = {
      'external_db_id' =>
        $XrefMapper::BasicMapper::source_to_external_db_id{ $xref->{
          'source_id'} }
        || 11000,    # FIXME ('external_db' needs to be updated)
      'accession' => $xref->{'accession'},
      'reason'    => 'No overlap',
      'reason_full' =>
        'No coordinate overlap with any Ensembl gene or transcript' };
  }
  $xref_sth->finish();

  ######################################################################
  # Do coordinate matching.                                            #
  ######################################################################

  my $sql = qq(
    SELECT    coord_xref_id, accession,
              txStart, txEnd,
              cdsStart, cdsEnd,
              exonStarts, exonEnds
    FROM      coordinate_xref
    WHERE     species_id = ?
    AND       chromosome = ? AND strand = ?
    AND       ((txStart >= ? AND txStart <= ?)    -- txStart in region
    OR         (txEnd   >= ? AND txEnd   <= ?)    -- txEnd in region
    OR         (txStart <= ? AND txEnd   >= ?))   -- region is contained
    ORDER BY  accession
  );

  foreach my $chromosome (@chromosomes) {
    my $chr_name = $chromosome->seq_region_name();

    log_progress( "Processing chromsome '%s'\n", $chr_name );

    my @genes = @{ $chromosome->get_all_Genes( undef, undef, 1 ) };

    log_progress( "There are %4d genes on chromosome '%s'\n",
                  scalar(@genes), $chr_name );

    while ( my $gene = shift(@genes) ) {
      my @transcripts = @{ $gene->get_all_Transcripts() };

      my %gene_result;

      foreach my $transcript ( sort { $a->start() <=> $b->start() }
                               @transcripts )
      {
        my @exons = @{ $transcript->get_all_Exons() };

        my %transcript_result;

        my $rr1 = Bio::EnsEMBL::Mapper::RangeRegistry->new();

        my $coding_transcript;
        if ( defined( $transcript->translation() ) ) {
          $coding_transcript = 1;
        } else {
          $coding_transcript = 0;
        }

        foreach my $exon (@exons) {
          $rr1->check_and_register( 'exon', $exon->start(),
                                    $exon->end() );

          if (    $coding_transcript
               && defined( $exon->coding_region_start($transcript) )
               && defined( $exon->coding_region_end($transcript) ) )
          {
            $rr1->check_and_register(
                                'coding',
                                $exon->coding_region_start($transcript),
                                $exon->coding_region_end($transcript) );
          }
        }

        my $sth = $xref_dbh->prepare_cached($sql);
        $sth->execute( $species_id,        $chr_name,
                       $gene->strand(),    $transcript->start(),
                       $transcript->end(), $transcript->start(),
                       $transcript->end(), $transcript->start(),
                       $transcript->end() );

        my ( $coord_xref_id, $accession, $txStart, $txEnd, $cdsStart,
             $cdsEnd, $exonStarts, $exonEnds );

        $sth->bind_columns(
                        \( $coord_xref_id, $accession, $txStart, $txEnd,
                           $cdsStart, $cdsEnd, $exonStarts, $exonEnds
                        ) );

        my $exonCount = 0;

        while ( $sth->fetch() ) {
          my @exonStarts = split( /,\s*/, $exonStarts );
          my @exonEnds   = split( /,\s*/, $exonEnds );
          my $exonCount = scalar(@exonStarts);

          my $rr2 = Bio::EnsEMBL::Mapper::RangeRegistry->new();

          my $exon_match   = 0;
          my $coding_match = 0;

          my $coding_count = 0;

          for ( my $i = 0 ; $i < $exonCount ; ++$i ) {
            my $overlap =
              $rr1->overlap_size( 'exon', $exonStarts[$i],
                                  $exonEnds[$i] );

            $exon_match +=
              $overlap/( $exonEnds[$i] - $exonStarts[$i] + 1 );

            $rr2->check_and_register( 'exon', $exonStarts[$i],
                                      $exonEnds[$i] );

            if ( !defined($cdsStart) || !defined($cdsEnd) ) {
              # Non-coding transcript.
            } else {
              my $codingStart = (   $exonStarts[$i] > $cdsStart
                                  ? $exonStarts[$i]
                                  : $cdsStart );
              my $codingEnd =
                ( $exonEnds[$i] < $cdsEnd ? $exonEnds[$i] : $cdsEnd );

              if ( $codingStart < $codingEnd ) {
                my $coding_overlap =
                  $rr1->overlap_size( 'coding', $codingStart,
                                      $codingEnd );

                $coding_match +=
                  $coding_overlap/( $codingEnd - $codingStart + 1 );

                $rr2->check_and_register( 'coding', $codingStart,
                                          $codingEnd );

                ++$coding_count;
              }
            }
          } ## end for ( my $i = 0 ; $i < ...

          my $rexon_match   = 0;
          my $rcoding_match = 0;

          my $rcoding_count = 0;

          foreach my $exon (@exons) {
            my $overlap =
              $rr2->overlap_size( 'exon', $exon->start(),
                                  $exon->end() );

            $rexon_match +=
              $overlap/( $exon->end() - $exon->start() + 1 );

            if (    $coding_transcript
                 && defined( $exon->coding_region_start($transcript) )
                 && defined( $exon->coding_region_end($transcript) ) )
            {
              my $coding_overlap =
                $rr2->overlap_size( 'coding',
                                    $exon->coding_region_start(
                                                           $transcript),
                                    $exon->coding_region_end(
                                                            $transcript)
                );

              $rcoding_match +=
                $coding_overlap/
                ( $exon->coding_region_end($transcript) -
                  $exon->coding_region_start($transcript) +
                  1 );

              ++$rcoding_count;
            }
          } ## end foreach my $exon (@exons)

          my $ucsc_exon_hit   = $exon_match/$exonCount;
          my $ucsc_coding_hit = $coding_match/( $coding_count || 1 );
          my $ens_exon_hit    = $rexon_match/scalar(@exons);
          my $ens_coding_hit  = $rcoding_match/( $rcoding_count || 1 );

          ##############################################################
          # Calculate the score.                                       #
          ##############################################################

          my $score = (
                ( $exon_match/$ens_weight + $ens_weight*$rexon_match )/
                  $coding_weight + $coding_weight*(
                  $coding_match/$ens_weight + $ens_weight*$rcoding_match
                  )
            )/( ( $exonCount/$ens_weight + $ens_weight*scalar(@exons) )/
                  $coding_weight + $coding_weight*(
                  $coding_count/$ens_weight + $ens_weight*$rcoding_count
                  ) );

          if ( !defined( $transcript_result{$coord_xref_id} )
               || $transcript_result{$coord_xref_id} < $score )
          {
            $transcript_result{$coord_xref_id} = $score;
          }

        } ## end while ( $sth->fetch() )
        $sth->finish();

        while ( my ( $coord_xref_id, $score ) =
                each(%transcript_result) )
        {
          if ( !defined( $gene_result{$coord_xref_id} )
               || $gene_result{$coord_xref_id} < $score )
          {
            $gene_result{$coord_xref_id} = $score;
          }
        }

        ################################################################
        # Apply transcript threshold.                                  #
        ################################################################

        foreach my $coord_xref_id ( keys(%transcript_result) ) {
          my $score = $transcript_result{$coord_xref_id};

          if ( $score > $transcript_score_threshold ) {
            if ( exists( $unmapped{$coord_xref_id} ) ) {
              $mapped{$coord_xref_id} = $unmapped{$coord_xref_id};
              delete( $unmapped{$coord_xref_id} );
              $mapped{$coord_xref_id}{'reason'}      = undef;
              $mapped{$coord_xref_id}{'reason_full'} = undef;
            }

            push( @{ $mapped{$coord_xref_id}{'mapped_to'} }, {
                    'ensembl_id'          => $transcript->dbID(),
                    'ensembl_object_type' => 'Transcript'
                  } );

          } elsif ( exists( $unmapped{$coord_xref_id} ) ) {
            $unmapped{$coord_xref_id}{'reason'} =
              'Did not meet threshold';
            $unmapped{$coord_xref_id}{'reason_full'} =
              sprintf( "Overlap score for transcript "
                         . "lower than threshold (%.2f)",
                       $transcript_score_threshold );
          }
        } ## end foreach my $coord_xref_id (...

        ################################################################
        # For this transcript, pick out the best matche(s) from the    #
        # ones passing the transcript threshold.                       #
        ################################################################

        # TODO

      } ## end foreach my $transcript ( sort...

      ##################################################################
      # Apply gene threshold.                                          #
      ##################################################################

      foreach my $coord_xref_id ( keys(%gene_result) ) {
        my $score = $gene_result{$coord_xref_id};

        if ( $score > $gene_score_threshold ) {
          if ( exists( $unmapped{$coord_xref_id} ) ) {
            $mapped{$coord_xref_id} = $unmapped{$coord_xref_id};
            delete( $unmapped{$coord_xref_id} );
            $mapped{$coord_xref_id}{'reason'}      = undef;
            $mapped{$coord_xref_id}{'reason_full'} = undef;
          }

          push( @{ $mapped{$coord_xref_id}{'mapped_to'} }, {
                  'ensembl_id'          => $gene->dbID(),
                  'ensembl_object_type' => 'Gene'
                } );

        } elsif ( exists( $unmapped{$coord_xref_id} ) ) {
          $unmapped{$coord_xref_id}{'reason'} =
            'Did not meet threshold';
          $unmapped{$coord_xref_id}{'reason_full'} =
            sprintf( "Overlap score for gene "
                       . "lower than threshold (%.2f)",
                     $gene_score_threshold );
        }
      } ## end foreach my $coord_xref_id (...

      ##################################################################
      # For this gene, pick out the best matche(s) from the ones       #
      # passing the gene threshold.                                    #
      ##################################################################

      # TODO

    } ## end while ( my $gene = shift(...
    last;
  } ## end foreach my $chromosome (@chromosomes)

  # Make all dumps.  Order is important.
  dump_xref( $output_dir, $xref_id, \%mapped, \%unmapped );
  dump_object_xref( $output_dir, $object_xref_id, \%mapped );
  dump_unmapped_reason( $output_dir, $unmapped_reason_id, \%unmapped );
  dump_unmapped_object( $output_dir,  $unmapped_object_id,
                        $analysis_id, \%unmapped );
}

sub dump_xref {
  my ( $output_dir, $xref_id, $mapped, $unmapped ) = @_;

  ######################################################################
  # Dump for 'xref'.                                                   #
  ######################################################################

  my $filename = catfile( $output_dir, 'xref_coord.txt' );

  my $fh = IO::File->new( '>' . $filename )
    or croak( sprintf( "Can not open '%s' for writing", $filename ) );

  log_progress( "Dumping for 'xref' to '%s'\n", $filename );

  foreach my $xref ( values( %{$unmapped} ), values( %{$mapped} ) ) {
    # Assign 'xref_id' to this Xref.
    $xref->{'xref_id'} = ++$xref_id;

    $fh->printf("%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n",
                $xref->{'xref_id'},
                $xref->{'external_db_id'},
                $xref->{'accession'},
                $xref->{'accession'},
                '0',
                '\N',
                'MISC',                               # FIXME
                '\N'                                  # FIXME (possibly)
    );
  }
  $fh->close();

  log_progress("Dumping for 'xref' done\n");

} ## end sub dump_xref

sub dump_object_xref {
  my ( $output_dir, $object_xref_id, $mapped ) = @_;

  ######################################################################
  # Dump for 'object_xref'.                                            #
  ######################################################################

  my $filename = catfile( $output_dir, 'object_xref_coord.txt' );

  my $fh = IO::File->new( '>' . $filename )
    or croak( sprintf( "Can not open '%s' for writing", $filename ) );

  log_progress( "Dumping for 'object_xref' to '%s'\n", $filename );

  foreach my $xref ( values( %{$mapped} ) ) {
    foreach my $object_xref ( @{ $xref->{'mapped_to'} } ) {
      # Assign 'object_xref_id' to this Object Xref.
      $object_xref->{'object_xref_id'} = ++$object_xref_id;

      $fh->printf( "%d\t%d\t%s\t%d\t%s\n",
                   $object_xref->{'object_xref_id'},
                   $object_xref->{'ensembl_id'},
                   $object_xref->{'ensembl_object_type'},
                   $xref->{'xref_id'},
                   '\N' );
    }
  }
  $fh->close();

  log_progress("Dumping for 'object_xref' done\n");

} ## end sub dump_objexref

sub dump_unmapped_reason {
  my ( $output_dir, $unmapped_reason_id, $unmapped ) = @_;

  ######################################################################
  # Dump for 'unmapped_reason'.                                        #
  ######################################################################

  # Create a list of the unique reasons.
  my %reasons;

  foreach my $xref ( values( %{$unmapped} ) ) {
    if ( !exists( $reasons{ $xref->{'reason_full'} } ) ) {
      $reasons{ $xref->{'reason_full'} } = {
                                        'summary' => $xref->{'reason'},
                                        'full' => $xref->{'reason_full'}
      };
    }
  }

  my $filename = catfile( $output_dir, 'unmapped_reason_coord.txt' );

  my $fh = IO::File->new( '>' . $filename )
    or croak( sprintf( "Can not open '%s' for writing", $filename ) );

  log_progress( "Dumping for 'unmapped_reason' to '%s'\n", $filename );

  foreach my $reason ( values(%reasons) ) {
    # Assign 'unmapped_reason_id' to this reason.
    $reason->{'unmapped_reason_id'} = ++$unmapped_reason_id;

    $fh->printf( "%d\t%s\t%s\n", $reason->{'unmapped_reason_id'},
                 $reason->{'summary'}, $reason->{'full'} );

  }
  $fh->close();

  log_progress("Dumping for 'unmapped_reason' done\n");

  # Assign reasons to the unmapped Xrefs from %reasons.
  foreach my $xref ( values( %{$unmapped} ) ) {
    $xref->{'reason'}      = $reasons{ $xref->{'reason_full'} };
    $xref->{'reason_full'} = undef;
  }

} ## end sub dump_unmapped_reason

sub dump_unmapped_object {
  my ( $output_dir, $unmapped_object_id, $analysis_id, $unmapped ) = @_;

  ######################################################################
  # Dump for 'unmapped_object'.                                        #
  ######################################################################

  my $filename = catfile( $output_dir, 'unmapped_object_coord.txt' );

  my $fh = IO::File->new( '>' . $filename )
    or croak( sprintf( "Can not open '%s' for writing", $filename ) );

  log_progress( "Dumping for 'unmapped_object' to '%s'\n", $filename );

  foreach my $xref ( values( %{$unmapped} ) ) {
    # Assign 'unmapped_object_id' to this Xref.
    $xref->{'unmapped_object_id'} = ++$unmapped_object_id;

    $fh->printf(
      "%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n",
      $xref->{'unmapped_object_id'},
      'xref',
      $analysis_id || -1,    # -1 means create analysis for
                             # this when uploading
      $xref->{'external_db_id'},
      $xref->{'accession'},
      $xref->{'reason'}->{'unmapped_reason_id'},
      '\N',
      '\N',
      '\N',
      '\N',
      '\N' );
  }
  $fh->close();

  log_progress("Dumping for 'unmapped_object' done\n");
} ## end sub dump_unmapped_object

sub log_progress {
  my ( $fmt, @params ) = @_;
  printf( STDERR "COORD==> %s", sprintf( $fmt, @params ) );
}

1;
