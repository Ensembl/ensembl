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

# $Id$

# This sub-class of XrefParser::BaseParser serves as the parent class
# for parsers of Xref source data that we use coordinate overlap to
# determine the Xrefs to.

package XrefParser::CoordinateParser;

use strict;
use warnings;

use Carp;
use DBI qw( :sql_types );

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

  $add_xref_sth->execute() or croak( $add_xref_sth->errstr() );
} ## end sub add_xref

1;
