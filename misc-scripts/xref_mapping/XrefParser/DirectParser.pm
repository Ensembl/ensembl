=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package XrefParser::DirectParser;

use strict;
use warnings;
use Carp;

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

sub run {
  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $filename     = $ref_arg->{file};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $filename) ){
    croak "Need to pass source_id, species_id and file as pairs";
  }
  $verbose |=0;


  my $file_io = $self->get_filehandle($filename);
  if ( !defined($file_io) ) {
    return 1;
  }

  my $parsed_count = 0;

  printf( STDERR "source = %d\t species = %d\n",
	  $source_id, $species_id );

  while ( defined( my $line = $file_io->getline() ) ) {
    chomp $line;

    my ( $accession, $ensembl_id, $type, $label, $description, $version )
      = split( /\t/, $line );

    if ( !defined($accession) || !defined($ensembl_id) ) {
      print {*STDERR} "Line $parsed_count contains  has less than two columns.\n";
      print {*STDERR} ("The parsing failed\n");
      return 1;
    }

    $type        ||= 'gene';
    $label       ||= $accession;
    $description ||= '';
    $version     ||= '1';

    ++$parsed_count;

    my $xref_id =
      $self->get_xref( $accession, $source_id, $species_id );

    if ( !defined($xref_id) || $xref_id eq '' ) {
      $xref_id =
	$self->add_xref({ acc        => $accession,
			  version    => $version,
			  label      => $label,
			  desc	     => $description,
			  source_id  => $source_id,
			  species_id => $species_id,
			  info_type  => "DIRECT"} );
    }
    $self->add_direct_xref( $xref_id, $ensembl_id,
					     $type, $accession );
  } ## end while ( defined( my $line...

  printf( "%d direct xrefs succesfully parsed\n", $parsed_count ) if($verbose);

  $file_io->close();

  print "Done\n" if($verbose);;

  return 0;
} ## end sub run

1;
