# $Id$

package XrefParser::VegaParser;

use warnings;
use strict;

use base qw( XrefParser::BaseParser );

# Parses the Vega Fasta file format:
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
            my (
                $vega_protein_id, $vega_type,
                $vega_position,   $vega_gene_id,
                $vega_transcript_id
            ) = split / /, $line;

            substr( $vega_protein_id, 0, 1, '' ); # Remove initial '>',
            substr( $vega_gene_id,    0, 5, '' ); # initial 'Gene:', and
            substr( $vega_transcript_id, 0, 11, '' );  #  'Transcript:'.

            my ( $vega_alphabet, $vega_status ) =
              ( $vega_type =~ /(.*):(.*)/ );

            my %xref = (
                'ACCESSION' => $vega_transcript_id,
                'LABEL'     => $vega_transcript_id,
                'DESCRIPTION' =>
                  sprintf( "%s %s", $vega_type, $vega_position ),
                'SEQUENCE'   => '',
                'SOURCE_ID'  => $source_id,
                'SPECIES_ID' => $species_id,
                'SEQUENCE_TYPE' =>
                  ( $vega_alphabet eq 'pep' ? 'peptide' : 'dna' ),
                'STATUS' => $vega_status,
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
