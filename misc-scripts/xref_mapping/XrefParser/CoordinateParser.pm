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
    my $dbh = $self->dbi();
    $add_xref_sth = $dbh->prepare_cached($add_xref_sql);
    if ( !defined($add_xref_sth) ) {
      croak( $dbh->errstr() );
    }
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

  $add_xref_sth->bind_param( 1,  $source_id,            SQL_INTEGER );
  $add_xref_sth->bind_param( 2,  $species_id,           SQL_INTEGER );
  $add_xref_sth->bind_param( 3,  $xref->{'accession'},  SQL_VARCHAR );
  $add_xref_sth->bind_param( 4,  $xref->{'chromosome'}, SQL_VARCHAR );
  $add_xref_sth->bind_param( 5,  $xref->{'strand'},     SQL_INTEGER );
  $add_xref_sth->bind_param( 6,  $xref->{'txStart'},    SQL_INTEGER );
  $add_xref_sth->bind_param( 7,  $xref->{'txEnd'},      SQL_INTEGER );
  $add_xref_sth->bind_param( 8,  $xref->{'cdsStart'},   SQL_INTEGER );
  $add_xref_sth->bind_param( 9,  $xref->{'cdsEnd'},     SQL_INTEGER );
  $add_xref_sth->bind_param( 10, $xref->{'exonStarts'}, SQL_VARCHAR );
  $add_xref_sth->bind_param( 11, $xref->{'exonEnds'},   SQL_VARCHAR );

  $add_xref_sth->execute();
    or croak( $add_xref_sth->errstr() );
} ## end sub add_xref

1;
