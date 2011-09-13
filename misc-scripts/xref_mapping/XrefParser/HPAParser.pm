package XrefParser::HPAParser;

use strict;
use warnings;
use Carp;
use base qw( XrefParser::BaseParser);

# This parser will read direct xrefs from a simple comma-delimited file downloaded from the Human Protein Atlas (HPA) database.
# The database contains two types of antibody, their own HPA antibodies and Collaborator antibody (CAB) commercial antibodies. 
# The columns of the file should be the following:
#
# 1)    Antibody
# 2)    Antibody ID
# 3)    Ensembl Peptide ID
# 4)	Link (URL)


sub run {
  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id, files and rel_file as pairs";
  }
  $verbose |=0;

  my $filename = @{$files}[0];

  my $file_io = $self->get_filehandle($filename);
  if ( !defined($file_io) ) {
    return 1;
  }

  my $parsed_count = 0;

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
	
    my $label       = $antibody;
    my $type        = 'translation';
    my $description = '';
    my $version     = '1';

    ++$parsed_count;

    my $xref_id = XrefParser::BaseParser->get_xref( $antibody_id, $source_id, $species_id );

    if ( !defined($xref_id) || $xref_id eq '' ) {
      $xref_id = XrefParser::BaseParser->add_xref($antibody_id, $version, $label, $description, $source_id, $species_id, "DIRECT");
    }
	
	
    XrefParser::BaseParser->add_direct_xref( $xref_id, $ensembl_peptide_id, $type, '');
	
  } ## end while ( defined( my $line...

  printf( "%d direct xrefs succesfully parsed\n", $parsed_count ) if($verbose);

  $file_io->close();

  return 0;
} ## end sub run

1;
