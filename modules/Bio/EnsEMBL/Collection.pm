# $Id$

package Bio::EnsEMBL::Collection;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument  ('rearrange');
use Bio::EnsEMBL::Utils::Exception ('throw');

use constant { FEATURE_DBID        => 0,
               FEATURE_SEQREGIONID => 1,
               FEATURE_START       => 2,
               FEATURE_END         => 3,
               FEATURE_STRAND      => 4 };

sub _remap {
  my ( $this, $features, $mapper, $slice ) = @_;

  if ( scalar( @{$features} ) > 0 ) {
    if ( $slice->get_seq_region_id() !=
         $features->[0][FEATURE_SEQREGIONID] )
    {
      throw(   '_remap() in Bio::EnsEMBL::Collection '
             . 'was unexpectedly expected to be intelligent' );
    }
  }

  return $features;
}

sub _create_feature {
  my ( $this, $feature_type, $args ) = @_;

  my ( $dbid, $start, $end, $strand, $slice ) =
    rearrange( [ 'DBID', 'START', 'END', 'STRAND', 'SLICE' ],
               %{$args} );
  my $feature =
    [ $dbid, $slice->get_seq_region_id(), $start, $end, $strand ];

  return $feature;
}

sub _create_feature_fast {
  my ( $this, $feature_type, $args ) = @_;

  if ( !defined( $this->{'recent_slice'} )
       || $this->{'recent_slice'} ne $args->{'slice'} )
  {
    $this->{'recent_slice'} = $args->{'slice'};
    $this->{'recent_slice_seq_region_id'} =
      $args->{'slice'}->get_seq_region_id();
  }

  my $feature = [ $args->{'dbID'},
                  $this->{'recent_slice_seq_region_id'},
                  $args->{'start'},
                  $args->{'end'},
                  $args->{'strand'} ];

  return $feature;
}

1;
