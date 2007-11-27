# $Id$

package XrefParser::FlybaseParser;

use strict;
use warnings;

use Data::Dumper;

use base qw( XrefParser::BaseParser );

our %object_types = ( mRNA       => 1,
                      miRNA      => 1,
                      ncRNA      => 1,
                      pseudogene => 1,
                      rRNA       => 1,
                      snRNA      => 1,
                      snoRNA     => 1,
                      tRNA       => 1,
                      gene       => 1 );

sub run {
  my $self = shift;
  my ( $source_id, $species_id, $data_file, $release_file ) = @_;

  my $data_io = $self->get_filehandle($data_file);

  while ( defined( my $line = $data_io->getline() ) ) {
    chomp($line);

    # Skip comment lines at the start of the file.
    if ( $line =~ /^#/ ) { next }

    # Split each line into fields.
    my @fields = split( /\t/, $line );

    # Only pick out the interesting lines.
    if ( !exists( $object_types{ $fields[2] } ) ) { next }

    # Go though each attribute (from the 9th field), split them up into
    # key-value pairs and store them.
    my %attributes;
    foreach my $attribute ( split( /;/, $fields[8] ) ) {
      my ( $key, $value ) = split( /=/, $attribute );
      if ( $key ne '' && $value ne '' ) {
        $attributes{$key} = $value;
      }
    }

    # For the 'Dbxref' and 'Ontology_term' attributes, split them up on
    # commas, divide into key-value pairs, and store them.
    foreach my $attribute_key ( 'Dbxref', 'Ontology_term' ) {
      if ( exists( $attributes{$attribute_key} ) ) {
        my %tmphash;
        foreach
          my $subattribute ( split( /,/, $attributes{$attribute_key} ) )
        {
          my ( $key, $value ) = split( /:/, $subattribute );
          push( @{ $tmphash{$key} }, $value );
        }

        # Replace the attribute entry with the hash.
        $attributes{$attribute_key} = \%tmphash;
      }
    }

    print( Dumper( \%attributes ) );
  }
  $data_io->close();
} ## end sub run

1;
