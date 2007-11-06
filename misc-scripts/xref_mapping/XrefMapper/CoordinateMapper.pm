# $Id$

# This is a set of subroutines used for creating Xrefs based on
# coordinate overlaps.

package XrefMapper::CoordinateMapper;

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Mapper::RangeRegistry;

use base qw( Exporter );

our @EXPORT = qw( run_coordinatemapping );

sub run_coordinatemapping {
  my $self = shift;

  #local $| = 1;

  my $xref_db = $self->xref();
  my $core_db = $self->core();

  my $species = $core_db->species();
  my $species_id =
    XrefMapper::BasicMapper::get_species_id_from_species_name( $xref_db,
                                                             $species );

  # We only do coordinate mapping for mouse and human for now.
  if ( !( $species eq 'mus_musculus' || $species eq 'homo_sapiens' ) ) {
    return;
  }

  my $xref_dbh = $xref_db->dbc()->db_handle();

  my $core_db_adaptor =
    Bio::EnsEMBL::DBSQL::DBAdaptor->new( -dbconn => $core_db->dbc() );
  my $slice_adaptor = $core_db_adaptor->get_SliceAdaptor();
  my @chromosomes   = @{ $slice_adaptor->fetch_all('Chromosome') };

  my $sql = qq(
    SELECT    accession,
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

      foreach my $transcript ( sort { $a->start() <=> $b->start() }
                               @transcripts )
      {
        my @exons = @{ $transcript->get_all_Exons() };

        my $rr1 = Bio::EnsEMBL::Mapper::RangeRegistry->new();

        log_progress( '=' x 80, "\n" );
        log_progress( "%15s %12s [%d,%d] [%d,%d] chr%s\n",
                $transcript->stable_id(), (
                  defined( $transcript->external_name() )
                  ? $transcript->external_name()
                  : ''
                ),
                $transcript->start(),
                $transcript->end(),
                $transcript->coding_region_start() || 0,
                $transcript->coding_region_end()   || 0,
                $chr_name );
        print("\n");

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

        my ( $accession, $txStart,    $txEnd, $cdsStart,
             $cdsEnd,    $exonStarts, $exonEnds );

        $sth->bind_columns(
                            \( $accession, $txStart, $txEnd, $cdsStart,
                               $cdsEnd, $exonStarts, $exonEnds
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

          my $coding_weight = 2;
          my $ens_weight    = 3;

          my $score = (
                ( $exon_match/$ens_weight + $ens_weight*$rexon_match )/
                  $coding_weight + $coding_weight*(
                  $coding_match/$ens_weight + $ens_weight*$rcoding_match
                  )
            )/( ( $exonCount/$ens_weight + $ens_weight*scalar(@exons) )/
                  $coding_weight + $coding_weight*(
                  $coding_count/$ens_weight + $ens_weight*$rcoding_count
                  ) );

          log_progress( "%10s\t"
                    . "%3d%% (%4.1f/%2d):"
                    . "%3d%% (%4.1f/%2d)|"
                    . "%3d%% (%4.1f/%2d):"
                    . "%3d%% (%4.1f/%2d) "
                    . "%3d%%\n",
                  $accession,
                  int( 100*$ucsc_exon_hit + 0.5 ),
                  $exon_match,
                  $exonCount,
                  int( 100*$ucsc_coding_hit + 0.5 ),
                  $coding_match,
                  $coding_count,
                  int( 100*$ens_exon_hit + 0.5 ),
                  $rexon_match,
                  scalar(@exons),
                  int( 100*$ens_coding_hit + 0.5 ),
                  $rcoding_match,
                  $rcoding_count,
                  int( 100*$score + 0.5 ) );

        } ## end while ( $sth->fetch() )
        $sth->finish();

      } ## end foreach my $transcript ( sort...
    } ## end while ( my $gene = shift(...

  } ## end foreach my $chromosome (@chromosomes)
} ## end sub run_coordinatemapping

sub log_progress {
  my ( $fmt, @params ) = @_;
  printf( STDERR "COORD==> %s", sprintf( $fmt, @params ) );
}

1;
