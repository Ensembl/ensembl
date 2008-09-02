package XrefParser::DirectParser;

use strict;

use base qw( XrefParser::BaseParser );

# This parser will read Direct Xrefs from a simple tab-delimited file.
# The columns of the file should be the following:
#
# 1)    Accession ID
# 2)    Ensembl ID
# 3)    Type (one of 'transcript', 'gene', or 'translation',
#       will default to 'gene')
# 4)    Display label (will default to the Accession ID in column 1)
# 5)    Description (will default to the empty string)
# 6)    Version (will default to '1')
#
# Columns 1 and 2 are obligatory.

sub new {
    my $proto = shift;

    my $class = ref $proto || $proto;
    my $self = bless {}, $class;

    return $self;
}

sub run {
    my $self = shift;

    my ( $source_id, $species_id, $filename ) = @_;

    my $file_io = $self->get_filehandle($filename);
    if ( !defined($file_io) ) {
        return 1;
    }

    my $parsed_count = 0;

    printf( STDERR "source = %d\t species = %d\n",
            $source_id, $species_id );

    while ( defined( my $line = $file_io->getline() ) ) {
        chomp $line;

        my ( $accession, $ensembl_id, $type, $label, $description,
             $version )
          = split( /\t/, $line );

        if ( !defined($accession) || !defined($ensembl_id) ) {
            printf( "Line %d contains  has less than two columns.\n",
                    1 + $parsed_count );
            print("The parsing failed\n");
            return 1;
        }

        $type        ||= 'gene';
        $label       ||= $accession;
        $description ||= '';
        $version     ||= '1';

        ++$parsed_count;

        my $xref_id =
          XrefParser::BaseParser->get_xref( $accession, $source_id, $species_id );

        if ( !defined($xref_id) || $xref_id eq '' ) {
            $xref_id =
              XrefParser::BaseParser->add_xref(
                                   $accession,   $version,   $label,
                                   $description, $source_id, $species_id
              );
        }
        XrefParser::BaseParser->add_direct_xref( $xref_id, $ensembl_id,
                                                 $type, $accession );
    } ## end while ( defined( my $line...

    printf( "%d direct xrefs succesfully parsed\n", $parsed_count );

    $file_io->close();

    print "Done\n";

    return 0;
} ## end sub run

1;
