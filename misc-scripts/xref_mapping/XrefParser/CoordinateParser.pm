# $Id$

# This sub-class of XrefParser::BaseParser serves as the parent class
# for parsers of Xref source data that we use coordinate overlap to
# determine the Xrefs to.

package XrefParser::CoordinateParser;

use strict;
use warnings;

use Carp;

use base qw( XrefParser::BaseParser );

our $add_xref_sth;
our $add_xref_sql = q(
  INSERT INTO coordinate_xref
  ( source_id,  species_id,
    accession,
    chromosome, strand,
    txStart,    txEnd,
    cdsStart,   cdsEnd,
    exonStarts, exonEnds )
  VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
);

sub add_xref {
  my $self = shift;
  my ( $source_id, $species_id, $xref ) = @_;

  if ( !defined($add_xref_sth) ) {
    $add_xref_sth = $self->dbi()->prepare_cached($add_xref_sql);
  }

  for my $required_key ( 'accession', 'chromosome',
                         'strand',    'txStart',
                         'txEnd',     'exonStarts',
                         'exonEnds' )
  {
    if ( !defined( $xref->{$required_key} ) ) {
      croak(
          sprintf( "Missing required key '%s' for Xref", $required_key )
      );
    }
  }

  $add_xref_sth->execute( $source_id,           $species_id,
                          $xref->{'accession'}, $xref->{'chromosome'},
                          $xref->{'strand'},    $xref->{'txStart'},
                          $xref->{'txEnd'},     $xref->{'cdsStart'},
                          $xref->{'cdsEnd'},    $xref->{'exonStarts'},
                          $xref->{'exonEnds'}
  ) or croak( $add_xref_sth->errstr() );
}

1;
