package XrefParser::HPAParser;

use strict;

use base qw( XrefParser::BaseParser);

# This parser will read direct xrefs from a simple comma-delimited file downloaded from the Human Protein Atlas (HPA) database.
# The database contains two types of antibody, their own HPA antibodies and Collaborator antibody (CAB) commercial antibodies. 
# The columns of the file should be the following:
#
# 1)    Antibody
# 2)    Antibody ID
# 3)    Ensembl Peptide ID
# 4)	Link (URL)


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

    $file_io->getline();
    
    while ( defined( my $line = $file_io->getline() ) ) {
 
	$line =~ s/\s*$//;
		

        my ( $antibody, $antibody_id, $ensembl_peptide_id, $link) = split( /,/, $line );

        if ( !defined($antibody) || !defined($ensembl_peptide_id) ) {
            printf( "Line %d contains  has less than two columns.\n",
                    1 + $parsed_count );
            print ("The parsing failed\n");
            return 1;
        }
	
        my $label       = $antibody_id;
	my $type        = 'translation';
        my $description = '';
        my $version     = '1';

        ++$parsed_count;

        my $xref_id = XrefParser::BaseParser->get_xref( $antibody, $source_id );

        if ( !defined($xref_id) || $xref_id eq '' ) {
            $xref_id = XrefParser::BaseParser->add_xref($antibody, $version, $label, $description, $source_id, $species_id);
        }
	
	
	XrefParser::BaseParser->add_direct_xref( $xref_id, $ensembl_peptide_id, $type, $antibody);
	
    } ## end while ( defined( my $line...

    printf( "%d direct xrefs succesfully parsed\n", $parsed_count );

    $file_io->close();

    print "Done\n";

    return 0;
} ## end sub run

1;
