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

# ----------------------------------------------------------------------
# The public interface added by Feature Collection

sub fetch_bins_by_Slice {
  my $this = shift;
  my ( $slice, $method, $nbins, $logic_name ) =
    rearrange( [ 'SLICE', 'METHOD', 'NBINS', 'LOGIC_NAME' ], @_ );

  my $old_value = $this->_lightweight();
  $this->_lightweight(1);

  my $bins =
    $this->_bin_features(
        -slice  => $slice,
        -nbins  => $nbins,
        -method => $method,
        -features =>
          $this->fetch_by_Slice_constraint( $slice, undef, $logic_name )
    );

  $this->_lightweight($old_value);

  return $bins;
}

# ----------------------------------------------------------------------
# Specialization of methods in Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor

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

# ----------------------------------------------------------------------
# Private methods

sub _lightweight {
  my ( $this, $value ) = @_;

  if ( defined($value) ) {
    $self->{'lightweight'} = ( $value != 0 );
  }

  return $self->{'lightweight'} || 0;
}

our %VALID_BINNING_METHODS = (
               'count'            => 0,
               'density'          => 0,    # Same as 'count'.
               'indices'          => 1,
               'index'            => 1,    # Same as 'indices'.
               'entries'          => 2,
               'entry'            => 2,    # Same as 'entries'.
               'fractional_count' => 3,
               'weight'           => 3,    # Same as 'fractional_count'.
               'coverage'         => 4 );

sub _bin_features {
  my $this = shift;
  my ( $slice, $nbins, $method_name, $features ) =
    rearrange( [ 'SLICE', 'NBINS', 'METHOD', 'FEATURES' ], @_ );

  if ( !defined($features) || !@{$features} ) { return [] }

  if ( !defined($nbins) ) {
    throw('Missing NBINS argument');
  } elsif ( $nbins <= 0 ) {
    throw('Negative or zero NBINS argument');
  }

  $method_name ||= 'count';
  if ( !exists( $VALID_BINNING_METHODS{$method_name} ) ) {
    throw(
           sprintf(
              "Invalid binning method '%s', valid methods are:\n\t%s\n",
              $method_name,
              join( "\n\t", sort( keys(%VALID_BINNING_METHODS) ) ) ) );
  }
  my $method = $VALID_BINNING_METHODS{$method_name};

  my $slice_start = $slice->start();

  my $bin_length = ( $slice->end() - $slice_start + 1 )/$nbins;

  my @bins = map( $_ = undef, 0 .. $nbins - 1 );

  my $feature_index = 0;
  my @bin_masks;

  foreach my $feature ( @{$features} ) {
    my $start_bin =
      int( ( $feature->[FEATURE_START] - $slice_start )/$bin_length );
    my $end_bin =
      int( ( $feature->[FEATURE_END] - $slice_start )/$bin_length );

    if ( $end_bin >= $nbins ) {
      # This might happen for the very last entry.
      #       $end_bin = $nbins - 1;
    }

    if ( $method == 0 ) {

      # For 'count' and 'density'.

      for ( my $bin_index = $start_bin ;
            $bin_index <= $end_bin ;
            ++$bin_index )
      {
        ++$bins[$bin_index];
      }

    } elsif ( $method == 1 ) {

      # For 'indices' and 'index'

      for ( my $bin_index = $start_bin ;
            $bin_index <= $end_bin ;
            ++$bin_index )
      {
        push( @{ $bins[$bin_index] }, $feature_index );
      }

      ++$feature_index;

    } elsif ( $method == 2 ) {

      # For 'entries' and 'entry'.

      for ( my $bin_index = $start_bin ;
            $bin_index <= $end_bin ;
            ++$bin_index )
      {
        push( @{ $bins[$bin_index] }, $feature );
      }

    } elsif ( $method == 3 ) {

      # For 'fractional_count' and 'weight'.

      if ( $start_bin == $end_bin ) {
        ++$bins[$start_bin];
      } else {

        my $feature_length =
          $feature->[FEATURE_END] - $feature->[FEATURE_START] + 1;

        # The first bin...
        $bins[$start_bin] +=
          ( ( $start_bin + 1 )*$bin_length -
            ( $feature->[FEATURE_START] - $slice_start ) )/
          $feature_length;

        # The intermediate bins (if there are any)...
        for ( my $bin_index = $start_bin + 1 ;
              $bin_index <= $end_bin - 1 ;
              ++$bin_index )
        {
          $bins[$bin_index] += $bin_length/$feature_length;
        }

        # The last bin...
        $bins[$end_bin] +=
          ( ( $feature->[FEATURE_END] - $slice_start ) -
            $end_bin*$bin_length +
            1 )/$feature_length;

      } ## end else [ if ( $start_bin == $end_bin)

    } elsif ( $method == 4 ) {

      # For 'coverage'.

      my $feature_start = $feature->[FEATURE_START] - $slice_start;
      my $feature_end   = $feature->[FEATURE_END] - $slice_start;

      if ( !defined( $bin_masks[$start_bin] )
           || ( defined( $bin_masks[$start_bin] )
                && $bin_masks[$start_bin] != 1 ) )
      {
        # Mask the $start_bin from the start of the feature to the end
        # of the bin, or to the end of the feature (whichever occurs
        # first).
        my $bin_start = int( $start_bin*$bin_length );
        my $bin_end = int( ( $start_bin + 1 )*$bin_length - 1 );
        for ( my $pos = $feature_start ;
              $pos <= $bin_end && $pos <= $feature_end ;
              ++$pos )
        {
          $bin_masks[$start_bin][ $pos - $bin_start ] = 1;
        }
      }

      for ( my $bin_index = $start_bin + 1 ;
            $bin_index <= $end_bin - 1 ;
            ++$bin_index )
      {
        # Mark the middle bins between $start_bin and $end_bin as fully
        # masked out.
        $bin_masks[$bin_index] = 1;
      }

      if ( $end_bin != $start_bin ) {

        if ( !defined( $bin_masks[$end_bin] )
             || ( defined( $bin_masks[$end_bin] )
                  && $bin_masks[$end_bin] != 1 ) )
        {
          # Mask the $end_bin from the start of the bin to the end of
          # the feature, or to the end of the bin (whichever occurs
          # first).
          my $bin_start = int( $end_bin*$bin_length );
          my $bin_end = int( ( $end_bin + 1 )*$bin_length - 1 );
          for ( my $pos = $bin_start ;
                $pos <= $feature_end && $pos <= $bin_end ;
                ++$pos )
          {
            $bin_masks[$end_bin][ $pos - $bin_start ] = 1;
          }
        }

      }

    } ## end elsif ( $method == 4 )

  } ## end foreach my $feature ( @{$features...

  if ( $method == 4 ) {

    # For the 'coverage' method: Finish up by going through @bin_masks
    # and sum up the arrays.

    for ( my $bin_index = 0 ; $bin_index < $nbins ; ++$bin_index ) {
      if ( defined( $bin_masks[$bin_index] ) ) {
        if ( !ref( $bin_masks[$bin_index] ) ) {
          $bins[$bin_index] = 1;
        } else {
          $bins[$bin_index] =
            scalar( grep ( defined($_), @{ $bin_masks[$bin_index] } ) )/
            $bin_length;
        }
      }
    }

  }

  return \@bins;
} ## end sub _bin_features

1;
