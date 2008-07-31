package XrefParser::DBASSParser;

use strict;

use base qw( XrefParser::BaseParser);
use DBI;
use Carp;

# This parser will read direct xrefs from a simple comma-delimited file downloaded from the DBASS database.
# The columns of the file should be the following:
#
# 1)    DBASS Gene ID
# 2)    DBASS Gene Name
# 3)    Ensembl Gene ID


my $dbi;
my $xref_id;
my $source_id;

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
		

        my ( $dbass_gene_id, $dbass_gene_name, $ensembl_id) = split( /,/, $line );

        if ( !defined($dbass_gene_id) || !defined($ensembl_id) ) {
            printf( "Line %d contains  has less than two columns.\n",
                    1 + $parsed_count );
            print ("The parsing failed\n");
            return 1;
        }
	
	
	my $first_gene_name;
	my $second_gene_name;
	
	
	if ($dbass_gene_name =~ /.\/./){
		($first_gene_name, $second_gene_name) = split( /\//, $dbass_gene_name );
	}
	
	if ($dbass_gene_name =~ /(.*)\((.*)\)/){
		$first_gene_name = $1;
	 	$second_gene_name = $2;
		print $first_gene_name, "\n", $second_gene_name, "\n";
	}
       
       
       
        my $label       = $first_gene_name;
	my $type        = 'gene';
        my $description = '';
        my $version     = '1';
	my $synonym	= $second_gene_name;

        ++$parsed_count;

        my $xref_id = XrefParser::BaseParser->get_xref( $dbass_gene_id, $source_id );

        if ( !defined($xref_id) || $xref_id eq '' ) {
            $xref_id = XrefParser::BaseParser->add_xref($dbass_gene_id, $version, $label, $description, $source_id, $species_id);
        }
	
	XrefParser::BaseParser->add_direct_xref( $xref_id, $ensembl_id, $type, '');
	
	if (defined ($synonym)) {
		XrefParser::DBASSParser->add_synonym($xref_id, $synonym);
	}
	elsif ($synonym =~ /^\s/){
		print "There is white space! \n";	
	}
	else {
		next;
	}
	
	
    } ## end while ( defined( my $line...

    printf( "%d direct xrefs succesfully parsed\n", $parsed_count );

    $file_io->close();

    print "Done\n";

    return 0;
} ## end sub run


1;
