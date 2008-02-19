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

  my ( $analysis, $stable_id, $version ) =
    rearrange( [ 'STABLE_ID', 'VERSION', ], %{$args} );

  push( @{$feature}, $stable_id, $version );

  return $feature;
}

sub _create_feature_fast {
  my ( $this, $feature_type, $args ) = @_;

  my $feature =
    $this->SUPER::_create_feature_fast( $feature_type, $args );

  return $feature;
}

1;
