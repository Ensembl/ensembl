# $Id$

package Bio::EnsEMBL::Collection;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument  ('rearrange');
use Bio::EnsEMBL::Utils::Exception ('throw');

sub _remap {
  my ( $features, $mapper, $slice ) = @_;

  return $features;
}

sub _create_feature {
  my ( $this, $feature_type, $args ) = @_;

  my ( $dbid, $start, $end, $strand, $slice ) =
    rearrange( [ 'DBID', 'START', 'END', 'STRAND', 'SLICE' ],
               @{$args} );

  return [ $dbid, $slice->get_seq_region_id(), $start, $end, $strand ];
}

sub _create_feature_fast {
  my ( $this, $feature_type, $args ) = @_;

  return [ $args->{'dbID'},  $args->{'slice'}->get_seq_region_id(),
           $args->{'start'}, $args->{'end'},
           $args->{'strand'} ];
}

1;
