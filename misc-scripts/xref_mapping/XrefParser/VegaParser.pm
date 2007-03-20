# $Id$

package XrefParser::VegaParser;

use warnings;
use strict;

use base qw( XrefParser::BaseParser );

# Parses the Vega CDNA and Peptide Fasta file format:
#
# >OTTMUST00000004500 cdna:tot chromosome:VEGA:1:60690948:60709172:1 Gene:OTTMUSG00000002254
# GTGACTTCAGTTCACACCACACTCTGCCTTGCTCACAGAGGAGGGGCTGCAGCCCTGGCC
# CTCATCAGAACAATGACACTCAGGCTGCTGTTCTTGGCTCTCAACTTCTTCTCAGTTCAA
# GTAACAGAAAACAAGATTTTGGTAAAGCAGTCGCCCCTGCTTGTGGTAGATAGCAACGAG
#
# >OTTMUSP00000002157 pep:known chromosome:VEGA:1:60690904:60717905:1 Gene:OTTMUSG00000002254 Transcript:OTTMUST00000004499
# MTLRLLFLALNFFSVQVTENKILVKQSPLLVVDSNEVSLSCRYSYNLLAKEFRASLYKGV
# NSDVEVCVGNGNFTYQPQFRSNAEFNCDGDFDNETVTFRLWNLHVNHTDIYFCKIEFMYP
# PPYLDNERSNGTIIHIKEKHLCHTQSSPKLFWALVVVAGVLFCYGLLVTVALCVIWTNSR
# RNRLLQSDYMNMTPRRPGLTRKPYQPYAPARDFAAYRP

sub run
{
    my $self = shift;
    my ( $source_id, $species_id, $file_name ) = @_;

    my $file_io = $self->get_filehandle($file_name);

    if ( !defined $file_io ) {
        return 1;    # Failed.
    }

    my @xrefs;
    while ( defined( my $line = $file_io->getline() ) ) {
        chomp $line;

        if ( substr( $line, 0, 1 ) eq '>' ) {
            # New sequence header.

            substr( $line, 0, 1, '' );    # Remove initial '>'

            my ( $vega_id, $vega_alphabet ) =
              ( $line =~ /^(\S+)\s([^:]+):/ );

            my %xref = (
                'ACCESSION'   => $vega_id,
                'LABEL'       => $vega_id,
                'DESCRIPTION' => $line,
                'SEQUENCE'    => '',
                'SOURCE_ID'   => $source_id,
                'SPECIES_ID'  => $species_id,
                'SEQUENCE_TYPE' =>
                  ( $vega_alphabet eq 'pep' ? 'peptide' : 'dna' ),
                'STATUS' => 'experimental'
            );

            push @xrefs, \%xref;

        } else {
            $xrefs[-1]->{'SEQUENCE'} .= $line;
        }
    }

    print scalar(@xrefs) . " Vega Fasta Xrefs successfully parsed\n";

    $self->upload_xref_object_graphs( \@xrefs );

    print "Done\n";

    return 0;    # Successful.
}

1;
