# $Id$

package Bio::EnsEMBL::Collection::Exon;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument  ('rearrange');
use Bio::EnsEMBL::Utils::Exception ('throw');

use base( 'Bio::EnsEMBL::Collection',
          'Bio::EnsEMBL::DBSQL::ExonAdaptor' );

sub _create_feature {
  my ( $this, $feature_type, $args ) = @_;

  my $feature = $this->SUPER::_create_feature( $feature_type, $args );

  if ( !$this->_lightweight() ) {
    my ( $analysis, $stable_id, $version ) =
      rearrange( [ 'STABLE_ID', 'VERSION', ], %{$args} );

    push( @{$feature}, $stable_id, $version );
  }

  return $feature;
}

sub _create_feature_fast {
  my ( $this, $feature_type, $args ) = @_;

  my $feature =
    $this->SUPER::_create_feature_fast( $feature_type, $args );

  return $feature;
}

sub _tables {
  my ($this) = @_;

  my @tables = $this->SUPER::_tables();

  if ( $this->_lightweight() ) {
    return ( $tables[0] );
  }

  return @tables;
}

sub _columns {
  my ($this) = @_;

  my @columns = $this->SUPER::_columns();

  if ( $this->_lightweight() ) {
    @columns[ 5 .. $#columns ] = map( 1, 5 .. $#columns );
  }

  return @columns;
}

sub _default_where_clause {
  my ($this) = @_;

  if ( $this->_lightweight() ) {
    return '';
  }

  return $this->SUPER::_default_where_clause();
}

1;
