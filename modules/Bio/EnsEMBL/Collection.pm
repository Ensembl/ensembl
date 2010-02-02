=head1 LICENSE

  Copyright (c) 1999-2010 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Collection - Abstract base class for feature collection
classes.

=head1 SYNOPSIS

  # Create an exon collection adaptor.
  my $collection =
    $registry->get_adaptor( 'Human', 'Core', 'ExonCollection' );

  # Create a slice.
  my $slice =
    $slice_adaptor->fetch_by_region( 'Chromosome', '2', 1, 300_000 );

  # Fetch the exons, represented as simple arrays of attributes, on
  # the slice.
  my @exons = @{ $collection->fetch_all_by_Slice($slice) };


  # Fetch the exons on the slice and count them into 50 bins (each
  # bin contains the number of exons overlapping that bin).
  my @bins = @{
    $collection->fetch_bins_by_Slice( -slice  => $slice,
                                      -nbins  => 50,
                                      -method => 'density'
    ) };

  # Fetch the exons on the slice and calculate the coverage of 50
  # bins (each bin contains the coverage of that bin by any exon as
  # a number between 0 and 1).
  my @bins = @{
    $collection->fetch_bins_by_Slice( -slice  => $slice,
                                      -nbins  => 50,
                                      -method => 'coverage'
    ) };

=head1 DESCRIPTION

This is the abstract base class for feature collections.

In short, this class replaces the object creation routines in
Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor with routines that creates
flat arrays.  It also provides a way of binning features on a slice in
various different ways (see below).

A feature collection provides a compact representation of features of
a particular type on a slice.  Each feature is represented by a short
array of attribute data.  This data is divided into two halfs:

=over 4

=item 1.

Basic feature representation.

=item 2.

Extended feature representation.

=back

=head2 Basic feature representation

The basic feature representation is common to all features of any type
and consists of a minimal set of values.  Each feature is represented by
an array that contains at least the following data (in this order)

=over 4

=item 1.

Ensembl internal database ID.

=item 2.

Ensembl internal sequence region ID.

=item 3.

Feature start position.

=item 4.

Feature end position.

=item 5.

Feature strand.

=back

The module defines a number of constants that may be used in place
of the index numbers 0 to 4: FEATURE_DBID, FEATURE_SEQREGIONID,
FEATURE_START, FEATURE_END, FEATURE_STRAND.  For a feature $f,
$f->[Bio::EnsEMBL::Collection::FEATURE_END] will thus be the end
position for the feature that the array @{$f} represents.

The position of the feature is in the same coordinate system as the
slice used to fetch the features.

=head2 Extended feature representation

A sub-class of this abstract base class will specify further
data to be added to the features in order to account for the
particular feature type.  A feature from a gene feature collection
(Bio::EnsEMBL::Collection::Gene) might, for example, contain the Ensembl
Stable ID of the gene.

The extended feature representation is defined by the method
_create_feature() or _create_feature_fast() which is implemented by the
sub-class.  The extended representation is also documented in the POD of
each sub-class.

=head2 Light-weight features

A light-weight feature is a feature whose representation does not
contain the extended feature representation.

Light-weight features may be fetched by setting the lightweight
attribute using the lightweight() method before calling e.g.
fetch_all_by_Slice().

=cut

package Bio::EnsEMBL::Collection;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument  ('rearrange');
use Bio::EnsEMBL::Utils::Exception ('throw');

use constant BASIC_SLOTS => 5;

use constant { FEATURE_DBID        => 0,
               FEATURE_SEQREGIONID => 1,
               FEATURE_START       => 2,
               FEATURE_END         => 3,
               FEATURE_STRAND      => 4 };

our %VALID_BINNING_METHODS = (
               'count'            => 0,
               'density'          => 0,    # Same as 'count'.
               'fractional_count' => 3,
               'weight'           => 3,    # Same as 'fractional_count'.
               'coverage'         => 4 );

# ======================================================================
# The public interface added by Feature Collection

sub lightweight {
  my ( $this, $value ) = @_;

  if ( defined($value) ) {
    $this->{'lightweight'} = ( $value != 0 );
  }

  return $this->{'lightweight'} || 0;
}

sub fetch_bins_by_Slice {
  my ( $this, @args ) = @_;

  my ( $slice, $method, $nbins, $logic_name, $lightweight ) =
    rearrange( [ 'SLICE', 'METHOD', 'NBINS', 'LOGIC_NAME' ], @args );

  # Temporarily set the colleciton to be lightweight.
  my $old_value = $this->lightweight();
  $this->lightweight(1);

  my $bins =
    $this->_bin_features(
           -slice    => $slice,
           -nbins    => $nbins,
           -method   => $method,
           -features => $this->fetch_all_by_Slice( $slice, $logic_name )
    );

  # Reset the lightweightness to whatever it was before.
  $this->lightweight($old_value);

  return $bins;
}

# ======================================================================
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

# ======================================================================
# Private methods

sub _bin_features {
  my ( $this, @args ) = @_;

  my ( $slice, $nbins, $method_name, $features ) =
    rearrange( [ 'SLICE', 'NBINS', 'METHOD', 'FEATURES' ], @args );

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

  my @bins;
  if ( $method == 0 ||    # 'count' or 'density'
       $method == 3 ||    # 'fractional_count' or 'weight'
       $method == 4       # 'coverage'
    )
  {
    # For binning methods where each bin contain numerical values.
    @bins = map( 0, 1 .. $nbins );
  } else {
    # For binning methods where each bin does not contain numerical
    # values.
    @bins = map( undef, 1 .. $nbins );
  }

  my $feature_index = 0;
  my @bin_masks;

  foreach my $feature ( @{$features} ) {
    my $start_bin =
      int( ( $feature->[FEATURE_START] - $slice_start )/$bin_length );
    my $end_bin =
      int( ( $feature->[FEATURE_END] - $slice_start )/$bin_length );

    if ( $end_bin >= $nbins ) {
      # This might happen for the very last entry.
      $end_bin = $nbins - 1;
    }

    if ( $method == 0 ) {
      # ----------------------------------------------------------------
      # For 'count' and 'density'.

      for ( my $bin_index = $start_bin;
            $bin_index <= $end_bin;
            ++$bin_index )
      {
        ++$bins[$bin_index];
      }

    }

    # Method 1 & 2 were deleted.

    elsif ( $method == 3 ) {
      # ----------------------------------------------------------------
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
        for ( my $bin_index = $start_bin + 1;
              $bin_index <= $end_bin - 1;
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
      # ----------------------------------------------------------------
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
        for ( my $pos = $feature_start;
              $pos <= $bin_end && $pos <= $feature_end;
              ++$pos )
        {
          $bin_masks[$start_bin][ $pos - $bin_start ] = 1;
        }
      }

      for ( my $bin_index = $start_bin + 1;
            $bin_index <= $end_bin - 1;
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
          for ( my $pos = $bin_start;
                $pos <= $feature_end && $pos <= $bin_end;
                ++$pos )
          {
            $bin_masks[$end_bin][ $pos - $bin_start ] = 1;
          }
        }

      }

    } ## end elsif ( $method == 4 )

  } ## end foreach my $feature ( @{$features...

  if ( $method == 4 ) {

    # ------------------------------------------------------------------
    # For the 'coverage' method: Finish up by going through @bin_masks
    # and sum up the arrays.

    for ( my $bin_index = 0; $bin_index < $nbins; ++$bin_index ) {
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
