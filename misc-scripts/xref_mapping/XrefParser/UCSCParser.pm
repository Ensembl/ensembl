# $Id$

package XrefParser::UCSCParser;

use strict;
use warnings;

use base qw( XrefParser::CoordinateParser );

sub run {
  my $self = shift;
  my ( $source_id, $species_id, $data_file, $release_file ) = @_;

  # Get the $source_id for the "UCSC" source.
  $source_id = $self->get_source_id_for_source_name('UCSC');

  my $data_io = $self->get_filehandle($data_file);

  while ( defined( my $line = $data_io->getline() ) ) {
    chomp($line);

    # Each line will have the following tab-delimited fields:
    # 0.    name        (UCSC stable ID)
    # 1.    chrom       (chromosome name, a la UCSC)
    # 2.    strand      (plus or minus)
    # 3.    txStart     (transcript start)
    # 4.    txEnd       (transcript end)
    # 5.    cdsStart    (CDS start)
    # 6.    cdsEnd      (CDS end)
    # 7.    exonCount   (number of exons in transcript)
    # 8.    exonStarts  (comma-separated list of exon start positions)
    # 9.    exonEnds    (comma-separated list of exon end positions)
    # 10.   proteinID   (cross reference to a protein ID, e.g. UniProt)
    # 11.   alignID     (not sure what this is right now)

    my ( $name,     $chrom,  $strand,     $txStart, $txEnd,
         $cdsStart, $cdsEnd, $exonStarts, $exonEnds
    ) = ( split( /\t/, $line ) )[ 0 .. 6, 8, 9 ];

    # UCSC uses slightly different chromosome names, at least for
    # human and mouse, so chop off the 'chr' in the beginning.  We do
    # not yet translate the names of the special chromosomes, e.g.
    # "chr6_cox_hap1" (UCSC) into "c6_COX" (Ensembl).
    $chrom =~ s/^chr//;

    # They also use '+' and '-' for the strand, instead of -1, 0, or 1.
    if    ( $strand eq '+' ) { $strand = 1 }
    elsif ( $strand eq '-' ) { $strand = -1 }
    else                     { strand  = 0 }

    # ... and non-coding transcripts have cdsStart == cdsEnd.  We would
    # like these to be stored as NULLs.
    if ( $cdsStart == $cdsEnd ) {
      undef($cdsStart);
      undef($cdsEnd);
    }

    # ... and they use the same kind of "inbetween" coordinates as e.g.
    # exonerate, so increment all start coordinates by one.
    $txStart += 1;
    $exonStarts =
      join( ',', map( { ++$_ } split( /,/, $exonStarts ) ) );
    if ( defined($cdsStart) ) { $cdsStart += 1 }

    # Cut off the last comma from $exonEnds, if it exists.  This is done
    # for $exonStarts already (above).
    if ( substr( $exonEnds, -1, 1 ) eq ',' ) { chop($exonEnds) }

    my %xref = ( 'accession'  => $name,
                 'chromosome' => $chrom,
                 'strand'     => $strand,
                 'txStart'    => $txStart,
                 'txEnd'      => $txEnd,
                 'cdsStart'   => $cdsStart,
                 'cdsEnd'     => $cdsEnd,
                 'exonStarts' => $exonStarts,
                 'exonEnds'   => $exonEnds );

    $self->add_xref( $source_id, $species_id, \%xref );
  } ## end while ( defined( my $line...
  $data_io->close();

  return 0;
} ## end sub run

1;
