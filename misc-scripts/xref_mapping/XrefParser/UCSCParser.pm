=head1 LICENSE

See the NOTICE file distributed with this work for additional information
regarding copyright ownership.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

XrefParser::UCSCParser

=head1 DESCRIPTION

A parser class to parse UCSC data for human and mouse.

-data_uri = ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/knownGene.txt.gz
-file_format = TSV
-columns = [
    ensembl_id
    chromosome
    strand
    tx_start
    tx_end
    cds_start
    cds_end
    nb_exons
    exon_starts
    exon_ends
    uniprot_accession
    ucsc_accession
  ]

Only columns listed in @required_columns are mandatory.

=head1 SYNOPSIS

  my $parser = XrefParser::UCSCParser->new($db->dbh);
  $parser->run({
    source_id  => 1,
    species_id => 9606,
    files      => ['UCSC_human/knownGene.txt.gz'],
  });

=cut

package XrefParser::UCSCParser;

use strict;
use warnings;

use Carp;
use Text::CSV;

use parent qw( XrefParser::CoordinateParser );


=head2 run
  Description: Runs the UCSCParser
  Return type: N/A
  Caller     : internal
=cut

sub run {
  my ( $self, $ref_arg ) = @_;

  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose} // 0;
  my $dbi          = $ref_arg->{dbi} // $self->{dbi};

  if ( (!defined $source_id) || (!defined $species_id) || (!defined $files) ) {
    confess 'Need to pass source_id, species_id and files as pairs';
  }

  my $file = @{$files}[0];

  my $count = 0;

  my $file_io = $self->get_filehandle($file);
  if ( !defined $file_io ) {
    confess "Can't open UCSC file $file\n";
  }

  my $input_file = Text::CSV->new({
    sep_char       => "\t",
    empty_is_undef => 1,
    strict         => 1,
  }) || confess "Cannot use file $file: " . Text::CSV->error_diag();

  while ( my $data = $input_file->getline( $file_io ) ) {
    my (undef, $chromosome, $strand, $tx_start, $tx_end, $cds_start, $cds_end,
        undef, $exon_starts, $exon_ends, undef, $accession) = @{ $data };

    # UCSC uses slightly different chromosome names, at least for
    # human and mouse, so chop off the 'chr' in the beginning.  We do
    # not yet translate the names of the special chromosomes, e.g.
    # "chr6_cox_hap1" (UCSC) into "c6_COX" (Ensembl).
    $chromosome =~ s{ \A chr }{}msx;

    # They also use '+' and '-' for the strand, instead of -1, 0, or 1.
    if ( $strand eq q{+} ) {
      $strand = 1;
    }
    elsif ( $strand eq q{-} ) {
      $strand = -1;
    }
    else {
      $strand = 0;
    }

    # ... and non-coding transcripts have cds_start == cds_end.  We would
    # like these to be stored as NULLs.
    if ( $cds_start == $cds_end ) {
      undef $cds_start;
      undef $cds_end;
    }

    # $exon_starts and $exon_ends usually (always?) have trailing commas,
    # remove them.
    $exon_starts =~ s{ , \z }{}msx;
    $exon_ends   =~ s{ , \z }{}msx;

    # ... and they use the same kind of "inbetween" coordinates as e.g.
    # exonerate, so increment all start coordinates by one.
    $tx_start += 1;
    if ( defined $cds_start ) {
      $cds_start += 1;
    }
    # The string exon_starts is a comma-separated list of start coordinates
    # for subsequent exons and we must increment each one. Split the string
    # on commas, use map() to apply the "+1" transformation to every
    # element of the resulting array, then join the result into a new
    # comma-separated list.
    $exon_starts =
      join q{,}, map { $_ + 1 } split qr{ , }msx, $exon_starts;

    $self->add_xref( $source_id, $species_id, {
      accession  => $accession,
      chromosome => $chromosome,
      strand     => $strand,
      txStart    => $tx_start,
      txEnd      => $tx_end,
      cdsStart   => $cds_start,
      cdsEnd     => $cds_end,
      exonStarts => $exon_starts,
      exonEnds   => $exon_ends,
      dbi        => $dbi,
    });

    $count += 1;

  }

  $input_file->eof || confess "Error parsing file $file: " . $input_file->error_diag();
  $file_io->close();

  if ($verbose) {
    print "Loaded a total of $count UCSC xrefs\n";
  }

  return 0;
}


1;
