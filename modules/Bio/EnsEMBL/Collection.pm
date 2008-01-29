# $Id$

# Ensembl module for Bio::EnsEMBL::Collection
#
# You may distribute this module under the same terms as perl itself.

=head1 NAME

Bio::EnsEMBL::Collection - Abstract base class for feature collection
classes.

=head1 SYNOPSIS

  use Bio::EnsEMBL::Collection::RepeatFeature;

  # Pick a slice.
  my $slice =
    $slice_adaptor->fetch_by_region( 'Chromosome', '2', 1, 1e9 );

  # Create a feature collection on the slice, restricting it to
  # 'Dust' features.
  my $collection =
    Bio::EnsEMBL::Collection::RepeatFeature->new(-slice => $slice,
                                                 -light => 0,
                                                 -analysis => 'Dust'
    );

  # Populate the feature collection from the slice.
  $collection->populate();

  # Populate the feature collection from the slice.  Sort the
  # entries on feature start position.
  $collection->populate( -sorted => 1 );

  # Populate the feature collection from the slice.  Sort the
  # entries on feature length.
  my $entry_start_idx = Bio::EnsEMBL::Collection::ENTRY_SEQREGIONSTART;
  my $entry_end_idx   = Bio::EnsEMBL::Collection::ENTRY_SEQREGIONEND;
  $collection->populate(
    -sorted  => 1,
    -sortsub => sub {
      $a->[$entry_end_idx] - $a->[$entry_start_idx]
      <=>
      $b->[$entry_end_idx] - $b->[$entry_start_idx];
    } );

  # Retrieve the entries from the collection as an array of arrays.
  my @entries = @{ $collection->entries() };

  # Retrieve a binned representation of the entries using 100 bins.
  my @binned_entries = @{
    $collection->get_bins( -nbins  => 100,
                           -method => 'entries'
    ) };

  # Retrieve a binned representation of the entries using 100 bins,
  # where each entry is represented by its index in the feature
  # collection array (@entries above).
  my @binned_entry_indicies = @{
    $collection->get_bins( -nbins  => 100,
                           -method => 'indices'
    ) };

  # Retrieve only the bin counts/densities.
  my @bin_counts = @{ $collection->get_bins( -nbins => 100 ) };

=head1 DESCRIPTION

This is the abstract base class for feature collections.

A feature collection provides a compact representation of features of a
particular type on a slice.  Each entry in a collection is a short array
of data representing a feature.  This data is divided into two halfs:

=over 4

=item 1.

Basic feature representation.

=item 2.

Extended feature representation.

=back

=head2 Basic feature representation

The basic feature representation is common to all entries in any type
of feature collection and consists of a minimal set of values.  Each
collection entry is an array that contains at least the following data
(in this order)

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

The module defines a number of constants that may be used
as symbolic constants in place of the index numbers 0 to
4: ENTRY_DBID, ENTRY_SEQREGIONID, ENTRY_SEQREGIONSTART,
ENTRY_SEQREGIONEND, ENTRY_SEQREGIONSTRAND.  For an entry $entry,
$entry->[Bio::EnsEMBL::Collection::ENTRY_SEQREGIONEND] will thus be the
end position for the feature that the entry represents.

The position of the feature is in the same coordinate system as the
slice associated with the collection object.

=head2 Extended feature representation

A sub-class of this abstract base class will specify further
data to be added to the entries in order to account for the
particular feature type.  An entry from a gene feature collection
(Bio::EnsEMBL::Collection::Gene) might, for example, contain the Ensembl
Stable ID of the gene.

The extended feature representation is defined by the method
_extra_columns() which is implemented by the sub-class.

=head2 Light-weight collections/entries

A light-weight collection is a feature collection whose collection
entries are light-weight, whose entries does not contain the extended
feature representation.

A collection which is light-weight may be created by using the argument
'-light=>1' in the constructor, new().

=head2 Binning methods

This module allows for various ways of binning the result of the
populate() method by using the get_bins() method.

Each binning method bins the collection entries and gives an array
with the specified length (number of bins).  An entry, which basically
consists of a start and a end position, is allocated to one or several
bins depending on its span and the size of the individual bins.

 Features:     |----|     |----------------|  |--|         |------|
               |-------------|                |-----|             |--|
               |-------------------|                 |----|       |--|

 Finer bins:   3 3 3 2 2 2 3 3 2 2 2 1 1 1 1 0 2 2 1 1 1 1 1 1 1 1 2 2
 Coarser bins: 3  3  2  2  3  2  2  1  1  1  0  2  1  1  1  1  1  3  2

The example above shows the arrays that might be returned from
get_bins() when using the 'count' binning method with 28 and 19 bins
respectively.

Here follows a brief description of each implemented binning method.

=over 4

=item  'count' and 'density'

Returns an array of bins, each bin containing the number of entries
allocated to (i.e. overlapping) that bin.  The 'density' binning method
is equivalent to 'count' and this is also the default binning method.

=item 'indices' and 'index'

Returns an array of bins, each bin containing an array of indices into
the collection entry array (available from the entries() method) for the
entries allocated to that bin.  The 'index' binning method is equivalent
to 'indicies'.

=item 'entries' and 'entry'

Returns an array of bins, each bin containing an array of entries
(references into the array of entries retrieved from the entries()
method) allocated to that bin.  The 'entry' binning method is equivalent
to 'entries'.

=item 'fractional_count' and 'fcount'

Returns an array of bins, each bin containing the sum of the fractions
of features overlapping that bin.  A feature fully inside a bin will
contribute 1 to the sum while a feature spanning exactly three bins
(from the very start of the first to the very end of the third) will
contribute 1/3 to the sum of each bin.

=item 'coverage'

Returns an array of bins, each bin containing the fraction of the bin
that is coverage by any feature.

=back

=head1 CONTACT

This modules is part of the Ensembl project.  See
http://www.ensembl.org/ for further information.

Questions may be posted to the ensembl-dev mailing list:
ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Collection;

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );

use base qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );

# Here are some class constants and class variables.

# Symbolic constants that acts as indices into the individual feature
# collection entries.  These must be in sync with the columns returned
# by the _columns() method.
use constant { ENTRY_DBID            => 0,
               ENTRY_SEQREGIONID     => 1,
               ENTRY_SEQREGIONSTART  => 2,
               ENTRY_SEQREGIONEND    => 3,
               ENTRY_SEQREGIONSTRAND => 4 };

# We'll keep the current slice and the projection segments of all slices
# on the various coordinate systems as class variables rather than
# as object attributes.  This way, the Ensembl drawing code can have
# several collection objects attached to one web view (e.g. one per
# track in ContigView), all sharing the same projection segments.  These
# structures may be emptied using the flush() method.

our $SLICE;
our %SEGMENTS;
our %SEQ_REG_MAP;

our %VALID_BINNING_METHODS = (
               'count'            => 0,
               'density'          => 0,    # Same as 'count'.
               'indices'          => 1,
               'index'            => 1,    # Same as 'indices'.
               'entries'          => 2,
               'entry'            => 2,    # Same as 'entries'.
               'fractional_count' => 3,
               'fcount'           => 3,    # Same as 'fractional_count'.
               'coverage'         => 4 );

=head1 METHODS (constructor)

=head2 new

  Arg [SLICE]   : Bio::EnsEMBL::Slice
                  The slice to be associated with this feature
                  collection.

  Arg [LIGHT]   : Boolean (optional, default false)
                  If true, the collection will be 'light-weight',
                  i.e. no type-specific data will be stored in its
                  entries when populate() is called.

  Arg [ANALYSIS]: String (optional, default undef)
                  Restrict the feature collection to a specific
                  analysis logic name, e.g. 'Dust' in a feature
                  collection of repeat features or 'ncRNA' in a gene
                  feature collection.

  Example       : my $collection =
                    Bio::EnsEMBL::Collection::<feature_type>->new(
                                                   -slice => $slice,
                                                   -light => 1 );

  Description   : When called for a sub-class, creates a feature
                  collection object and associates a slice with it.

  Return type   : Bio::EnsEMBL::Collection::<feature_type>

  Exceptions    : Throws if no slice is specified.
                  Warns if trying to restrict by analysis for a
                  feature type that does not have an analysis
                  associated with it, e.g. exon.

  Caller        : General (through a sub-class).

  Status        : At Risk (under development)

=cut

sub new {
  my $proto = shift;
  my ( $slice, $light, $analysis_logic_name ) =
    rearrange( [ 'SLICE', 'LIGHT', 'ANALYSIS' ], @_ );

  if ( !defined($slice) ) {
    throw(   'Unable to create a feature collection '
           . 'without a slice object, I am.' );
  }

  my $this = $proto->SUPER::new( $slice->adaptor()->db() );

  if ( defined($analysis_logic_name) ) {
    if ( !$this->_has_analysis() ) {
      warning(   'This feature type does not have '
               . 'an analysis to restrict by.' );
    } else {
      my $table_alias = $this->_feature_table()->[1];
      $this->__attrib( 'restrict_by',
                       sprintf( 'AND %s.analysis_id = a.analysis_id '
                                  . 'AND a.logic_name = %s',
                                $table_alias,
                                $this->dbc()->db_handle()
                                  ->quote($analysis_logic_name) ) );
    }
  }

  my $sql = qq(
    SELECT  cs.name, mc.max_length, m.meta_value
    FROM    coord_system cs,
            meta_coord mc
    LEFT JOIN meta m ON m.meta_key  = ?
    WHERE   mc.table_name           = ?
    AND     mc.coord_system_id      = cs.coord_system_id
  );

  my $sth           = $this->prepare($sql);
  my $feature_table = $this->_feature_table()->[0];
  $sth->execute( $feature_table . 'build.level', $feature_table );

  my ( $cs_name, $max_length, $build_level );
  $sth->bind_columns( \( $cs_name, $max_length, $build_level ) );

  my %coordinate_systems;

  while ( $sth->fetch() ) {
    $coordinate_systems{$cs_name} = {
            'name'               => $cs_name,
            'max_feature_length' => $max_length,
            'build_level'        => $build_level    # Not currently used
    };
  }
  $sth->finish();

  $this->__attrib( 'coordinate_systems', \%coordinate_systems );
  $this->__attrib( 'is_lightweight',     $light );

  $this->slice($slice);

  return $this;
} ## end sub new

#-----------------------------------------------------------------------

=head1 METHODS (public)

=head2 slice

  Arg [1]       : Bio::EnsEMBL::Slice (optional)
                  The new slice to be associated with this
                  collection.

  Example       : my $slice = $collection->slice();

  Description   : Getter/setter for the main slice associated with
                  this collection.  Associating a new slice with
                  the collection will empty all entries from the
                  collection.

  Return type   : Bio::EnsEMBL::Slice

  Exceptions    : None

  Caller        : General

  Status        : At Risk (under development)

=cut

sub slice {
  my ( $this, $slice ) = @_;

  if ( defined($slice) ) {
    my $slice_name = $slice->name();

    foreach
      my $cs_name ( keys( %{ $this->__attrib('coordinate_systems') } ) )
    {
      if ( exists( $SEGMENTS{$slice_name}{$cs_name} ) ) { next }

      my @segments = @{ $slice->project($cs_name) };

      foreach my $segment (@segments) {
        # This is a simple map between seq_region_id and a
        # ProjectionSegment, used in the mapping done by the
        # _objs_from_sth() method.
        $SEQ_REG_MAP{$slice_name}
          { $segment->to_Slice()->get_seq_region_id() } = $segment;
      }

      @segments = ( undef, [@segments] );
      $segments[0] = shift( @{ $segments[1] } );
      if ( scalar( @{ $segments[1] } ) > 0 ) {
        $segments[2] = pop( @{ $segments[1] } );
      }
      if ( scalar( @{ $segments[1] } ) == 0 ) {
        if ( exists( $segments[2] ) ) {
          $segments[1] = $segments[2];
        }
        pop(@segments);
      }

      # This list is a list of ProjectionSegment objects and of arrays
      # of ProjectionSegment objects.  For segments in arrays, no
      # constraint on seq_region_start or seq_region_end is needed.
      $SEGMENTS{$slice_name}{$cs_name} = \@segments;

    } ## end foreach my $cs_name ( keys(...

    $this->entries( [] );
    $this->__attrib( 'is_populated', 0 );

    $SLICE = $slice;

  } ## end if ( defined($slice) )

  return $SLICE;
} ## end sub slice

=head2 entries

  Arg [1]       : List reference
                  The new list of feature collection entries,
                  typically produced within the populate() method.

  Example       : my @entries = @{ $collection->entries() };

  Description   : Getter/setter for the list of collection entries.

  Return type   : List reference to list of list references

  Exceptions    : None

  Caller        : General

  Status        : At Risk (under development)

=cut

sub entries {
  my ( $this, $entries ) = @_;
  return $this->__attrib( 'entries', $entries );
}

=head2 is_populated

  Args          : None

  Example       : if ( !$collection->is_populated() ) {
                    $collection->populate();
                  }

  Description   : Returns a true value if the collection has been
                  populated, false of not.

  Return type   : Boolean

  Exceptions    : None

  Caller        : General

  Status        : At Risk (under development)

=cut

# Returns true if the collection has been populated, false of not.
sub is_populated {
  my ($this) = @_;
  return $this->__attrib('is_populated');
}

=head2 populate

  Arg [SORTED]  : Boolean (optional, default false)
                  If true, will sort the entries after popluating
                  the collection.

      [SORTSUB] : Subroutine reference (optional, default undef)
                  If defined, will be used to sort the entries after
                  popluating the collection.

  Example       : $collection->populate( -sorted => 1 );

  Description   : Populates the collection with a compact
                  representation of the features overlapping the
                  current slice, and optionally sorts the results.

  Return type   : None

  Exceptions    : None

  Caller        : General

  Status        : At Risk (under development)

=cut

sub populate {
  my $this = shift;
  my ( $sorted, $sort_like_this ) =
    rearrange( [ 'SORTED', 'SORTSUB' ], @_ );

  if ( $this->is_populated() ) { return }

  my @entries;

  my $slice_name = $this->slice()->name();

  foreach
    my $cs_name ( keys( %{ $this->__attrib('coordinate_systems') } ) )
  {
    foreach my $segment ( @{ $SEGMENTS{$slice_name}{$cs_name} } ) {
      if ( defined($segment)
           && ( ref($segment) ne 'ARRAY'
                || scalar( @{$segment} ) > 0 ) )
      {
        push( @entries,
              @{$this->generic_fetch( $this->__constraint($segment) ) }
        );
      }
    }
  }

  if ($sorted) {
    if ( !defined($sort_like_this) ) {
      $sort_like_this = sub {
        $a->[ENTRY_SEQREGIONSTART] <=> $b->[ENTRY_SEQREGIONSTART]
          || $a->[ENTRY_SEQREGIONEND] <=> $b->[ENTRY_SEQREGIONEND];
      };
    }

    @entries = sort( $sort_like_this @entries );
  }

  $this->entries( \@entries );
  $this->__attrib( 'is_populated', 1 );

} ## end sub populate

=head2 is_lightweight

  Arg [1]       : Boolean (optional)

  Example       : if ( !$collection->is_lightweight() ) { ... }

  Description   : Returns true if the collection was created as
                  a light-weight feature collection.  If the
                  collection is light-weight, its entries does not
                  contain any feature-specific data (e.g. transcript
                  or gene stable IDs for a transcript feature
                  collection).

  Return type   : Boolean

  Exceptions    : None

  Caller        : General

  Status        : At Risk (under development)

=cut

sub is_lightweight {
  my ($this) = @_;
  return $this->__attrib('is_lightweight');
}

=head2 count

  Args          : None

  Example       : my $count = $collection->count();

  Description   : Returns the number of features that overlaps the
                  current slice.

  Return type   : Integer

  Exceptions    : Throws if the collection has not yet been
                  populated.

  Caller        : General

  Status        : At Risk (under development)

=cut

sub count {
  my ($this) = @_;

  if ( !$this->is_populated() ) {
    throw(   'Can not count the entries in a feature collection '
           . 'without first having called populate().' );
  }

  return scalar( @{ $this->entries() } );
}

=head2 flush

  Args          : None

  Example       : $collection->flush();

  Description   : Flushes (deletes) all cached data not associated
                  with the current slice.

  Return type   : None

  Exceptions    : None

  Caller        : General

  Status        : At Risk (under development)

=cut

sub flush {
  my ($this) = @_;

  my $slice_name = $this->slice()->name();

  foreach my $key ( keys(%SEGMENTS) ) {
    if ( $key ne $slice_name ) {
      delete( $SEGMENTS{$key} );
    }
  }

  foreach my $key ( keys(%SEQ_REG_MAP) ) {
    if ( $key ne $slice_name ) {
      delete( $SEQ_REG_MAP{$key} );
    }
  }

}

=head2 get_bins

  Arg [NBINS]   : Integer
                  The number of bins to use.

  Arg [METHOD]  : String (optional, default 'count')
                  The binning method to use.  The valid methods are
                  described above, in the section called 'Binning
                  methods'.

  Example       : my @bins = @{
                    $collection->get_bins( -nbins  => 640,
                                           -method => 'count'
                    ) };

  Description   : Performs binning of the collection entries.

  Return type   : List reference

  Exceptions    : Throws if the population has not been populated
                  using the populate() method, or if the provided
                  binning method does not exist.

                  Throws if the NBINS argument is missing or if any
                  of the arguments are out of bounds.

  Caller        : General

  Status        : At Risk (under development)

=cut

sub get_bins {
  my $this = shift;
  my ( $nbins, $method_name ) = rearrange( [ 'NBINS', 'METHOD' ], @_ );

  if ( !$this->is_populated() ) {
    throw(   'Can not bin a feature collection '
           . 'without first having called populate()' );
  }

  if ( !defined($nbins) ) {
    throw('Missing NBINS argument');
  } elsif ( $nbins <= 0 ) {
    throw('Negative or zero NBINS argument');
  }

  $method_name ||= 'count';
  if ( !exists( $VALID_BINNING_METHODS{$method_name} ) ) {
    throw(
           sprintf("Invalid binning method '%s', valid methods are: %s",
                   $method_name,
                   join( ', ', sort( keys(%VALID_BINNING_METHODS) ) ) )
    );
  }
  my $method = $VALID_BINNING_METHODS{$method_name};

  my $slice       = $this->slice();
  my $slice_start = $slice->start();

  my $bin_length = ( $slice->end() - $slice_start + 1 )/$nbins;

  my @bins = map( $_ = undef, 0 .. $nbins - 1 );

  my $entry_index = 0;
  my @bin_masks;

  foreach my $entry ( @{ $this->entries() } ) {
    my $start_bin = int(
        ( $entry->[ENTRY_SEQREGIONSTART] - $slice_start )/$bin_length );
    my $end_bin = int(
          ( $entry->[ENTRY_SEQREGIONEND] - $slice_start )/$bin_length );

    if ( $end_bin >= $nbins ) {
      # This might happen for the very last entry.
      $end_bin = $nbins - 1;
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
        push( @{ $bins[$bin_index] }, $entry_index );
      }

      ++$entry_index;

    } elsif ( $method == 2 ) {

      # For 'entries' and 'entry'.

      for ( my $bin_index = $start_bin ;
            $bin_index <= $end_bin ;
            ++$bin_index )
      {
        push( @{ $bins[$bin_index] }, $entry );
      }

    } elsif ( $method == 3 ) {

      # For 'fractional_count' and 'fcount'.

      if ( $start_bin == $end_bin ) {
        ++$bins[$start_bin];
      } else {

        my $feature_length =
          $entry->[ENTRY_SEQREGIONEND] -
          $entry->[ENTRY_SEQREGIONSTART] + 1;

        # The first bin...
        $bins[$start_bin] +=
          ( ( $start_bin + 1 )*$bin_length -
            ( $entry->[ENTRY_SEQREGIONSTART] - $slice_start ) )/
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
          ( ( $entry->[ENTRY_SEQREGIONEND] - $slice_start ) -
            $end_bin*$bin_length +
            1 )/$feature_length;

      } ## end else [ if ( $start_bin == $end_bin)

    } elsif ( $method == 4 ) {

      # For 'coverage'.

      my $feature_start = $entry->[ENTRY_SEQREGIONSTART] - $slice_start;
      my $feature_end   = $entry->[ENTRY_SEQREGIONEND] - $slice_start;

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

  } ## end foreach my $entry ( @{ $this...

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
} ## end sub get_bins

#-----------------------------------------------------------------------

=head1 METHODS (private)

These are methods that should only ever be used by the class itself.

=head2 __attrib

  Args [1]      : String
                  The name of the attribute to get/set.

  Args [2]      : Any (optional)
                  The new value of the attribute.

  Description   : Simple generic getter/setter method for the
                  attributes of the class.

  Return type   : Varies

  Caller        : Various methods in this class.

=cut

sub __attrib {
  my ( $this, $attribute, $value ) = @_;

  if ( defined($value) ) {
    $this->{'attributes'}{$attribute} = $value;
  }

  return $this->{'attributes'}{$attribute};
}

=head2 __constraint

  Arg [1]       : Bio::EnsEMBL::ProjectionSegment or a reference to
                  list thereof.

  Description   : Produces the constraint for
                  Bio::EnsEMBL::DBSQL::BaseAdaptor::generic_fetch()

  Return type   : String

  Caller        : populate()

=cut

sub __constraint {
  my ( $this, $arg ) = @_;

  my $constraint;
  my $table_alias = $this->_feature_table()->[1];

  if ( ref($arg) ne 'ARRAY' ) {
    my $slice = $arg->to_Slice();

    my $max_feature_length =
      $this->__attrib('coordinate_systems')
      ->{ $slice->coord_system_name() }{'max_feature_length'};

    my $constraint_fmt = q(
        %1$s.seq_region_id     = %2$10d
    AND %1$s.seq_region_start <= %4$10d
    AND %1$s.seq_region_end   >= %3$10d
  );

    $constraint =
      sprintf( $constraint_fmt,
               $table_alias, $slice->get_seq_region_id(),
               $slice->start(), $slice->end() );

  } else {
    my @seq_region_ids =
      sort( { $a <=> $b }
            map( $_->to_Slice()->get_seq_region_id(), @{$arg} ) );

    my $constraint_fmt = '%s.seq_region_id IN (%s)';

    $constraint = sprintf( $constraint_fmt,
                           $table_alias, join( ',', @seq_region_ids ) );
  }

  return $constraint;
} ## end sub __constraint

#-----------------------------------------------------------------------

=head1 METHODS (protected)

These are methods that may be overridden (specialized) by sub-classes.

=head2 _extra_tables

  Args          : None

  Description   : The method provides a list of tables and
                  table aliases, in addition to the ones
                  returned by _tables(), that are used to
                  create collections of particular feature
                  types, e.g. the tables transcript_stable_id
                  [tsi], gene [g], and gene_stable_id [gsi] for
                  Bio::EnsEMBL::Collection::Transcript.

  Return type   : List of list references

  Exceptions    : None

  Caller        : _tables()

  Status        : At Risk (under development)

=cut

sub _extra_tables { return () }

=head2 _extra_columns

  Args          : None

  Description   : The method provides a list of columns, in
                  addition to the ones returned by _columns(),
                  that constitutes the contents of an entry of
                  a particular feature type, e.g. t.biotype,
                  t.status, tsi.stable_id, and gsi.stable_id for
                  Bio::EnsEMBL::Collection::Transcript.

  Return type   : List of strings

  Exceptions    : None

  Caller        : _columns()

  Status        : At Risk (under development)

=cut

sub _extra_columns { return () }

=head2 _extra_where_clause

  Args          : None

  Description   : The method should supply the necessary SQL for
                  joining the tables from _extra_tables().

  Return type   : String

  Exceptions    : None

  Caller        : _default_where_clause()

  Status        : At Risk (under development)

=cut

sub _extra_where_clause { return undef }

=head2 _dbID_column

  Args          : None

  Description   : The method provides the name of the column in the
                  primary feature table that holds the dbID for the
                  collection elements, i.e. 'transcript_id' for
                  Bio::EnsEMBL::Collection::Transcript.

                  If this method is not specialized by a sub-class,
                  it is assumed that the name of the dbID column is
                  the name of the primary feature table suffixed by
                  the string '_id'.

  Return type   : String

  Exceptions    : None

  Caller        : _columns()

  Status        : At Risk (under development)

=cut

sub _dbID_column {
  my ($this) = @_;

  if ( !defined( $this->__attrib('dbID_column') ) ) {
    $this->__attrib( 'dbID_column',
                     $this->_feature_table()->[0] . '_id' );
  }

  return $this->__attrib('dbID_column');
}

=head2 _has_analysis

  Args          : None

  Description   : Some feature types may have an analysis_id column
                  in their main feature table.  For these feature
                  types, a restriction on analysis.logic_name may
                  be enforced by using the ANALYSIS argument of the
                  constructor.  This method should return a false
                  value if there is no such an analysis_id column in
                  the feature table.

  Return type   : Boolean

  Exceptions    : None

  Caller        : new()

  Status        : At Risk (under development)

=cut

sub _has_analysis { return 1 }

#-----------------------------------------------------------------------

=head1 METHODS (abstract protected)

These are methods that needs to be implemented by a sub-class since this
base class can not implement them.

=head2 _feature_table

  Args          : None

  Description   : The method should return the primary
                  feature table and table alias used by the
                  collection, i.e. the transcript [t] table for
                  Bio::EnsEMBL::Collection::Transcript.

  Return type   : List reference to a list of two strings (the
                  feature table name and its alias).

  Exceptions    : Throws if called on the base class.

  Caller        : new(), _dbID_column(), __constraint()

  Status        : At Risk (under development)

=cut

sub _feature_table {
  throw("Called abstract method '_feature_table()'");
}

#-----------------------------------------------------------------------

=head1 Implemented abstract protected methods from base class
Bio::EnsEMBL::BaseAdaptor

=head2 _tables

  Args          : None

  Description   : Returns the primary feature table and, if the
                  collection is not light-weight, the additional
                  tables that are needed to construct the collection
                  entries.

  Return type   : List of list references

  Exceptions    : None

  Caller        : Bio::EnsEMBL::DBSQL::BaseAdaptor::generic_fetch()

  Status        : At Risk (under development)

=cut

sub _tables {
  my ($this) = @_;

  my @tables = ( $this->_feature_table() );

  if ( !$this->is_lightweight() ) {
    push( @tables, $this->_extra_tables() );
  }

  my $restrict_by = $this->__attrib('restrict_by');
  if ( defined($restrict_by) ) {
    push( @tables, [ 'analysis', 'a' ] );
  }

  return @tables;
}

=head2 _columns

  Args          : None

  Description   : Returns the columns that makes up a collection
                  entry.

  Return type   : List of strings

  Exceptions    : None

  Caller        : Bio::EnsEMBL::DBSQL::BaseAdaptor::generic_fetch()

  Status        : At Risk (under development)

=cut

sub _columns {
  my ($this) = @_;

  my $table_alias = $this->_feature_table()->[1];

  my @columns = ( $table_alias . '.' . $this->_dbID_column(),
                  $table_alias . '.seq_region_id',
                  $table_alias . '.seq_region_start',
                  $table_alias . '.seq_region_end',
                  $table_alias . '.seq_region_strand' );

  if ( !$this->is_lightweight() ) {
    push( @columns, $this->_extra_columns() );
  }

  return @columns;
}

=head2 _default_where_clause

  Args          : None

  Description   : For a light-weight collection, unrestricted by
                  any analysis, returns an empty string (no joins),
                  otherwise returns the WHERE clause that joins the
                  tables returned by _tables() and/or that restricts
                  the query by analysis.logic_name.

  Return type   : String

  Exceptions    : None

  Caller        : Bio::EnsEMBL::DBSQL::BaseAdaptor::generic_fetch()

  Status        : At Risk (under development)

=cut

sub _default_where_clause {
  my ($this) = @_;

  my $where_clause = '';

  if ( !$this->is_lightweight() ) {
    my $extra_where = $this->_extra_where_clause();
    if ( defined($extra_where) ) {
      $where_clause = $extra_where;
    }
  }

  my $restrict_by = $this->__attrib('restrict_by');
  if ( defined($restrict_by) ) {
    $where_clause .= ' ' . $restrict_by;
  }

  return $where_clause;
}

=head2 _objs_from_sth

  Arg [1]       : DBI statment handle

  Description   : Given a DBI statment handle, reads the entries,
                  maps them to the appropriate coordinate system,
                  and returns a list of them.

  Return type   : List reference to list of list references

  Exceptions    : None

  Caller        : Bio::EnsEMBL::DBSQL::BaseAdaptor::generic_fetch()

  Status        : At Risk (under development)

=cut

sub _objs_from_sth {
  my ( $this, $sth ) = @_;

  my @features;

  my ( $segment, $segment_slice, $segment_slice_start,
       $segment_slice_strand, $segment_offset );

  my $slice_name  = $this->slice()->name();
  my $slice_start = $this->slice()->start();

  while ( my $entry = $sth->fetchrow_arrayref() ) {
    if ( !defined($segment)
         || $segment !=
         $SEQ_REG_MAP{$slice_name}{ $entry->[ENTRY_SEQREGIONID] } )
    {
      $segment =
        $SEQ_REG_MAP{$slice_name}{ $entry->[ENTRY_SEQREGIONID] };
      $segment_slice        = $segment->to_Slice();
      $segment_slice_start  = $segment_slice->start();
      $segment_slice_strand = $segment_slice->strand();

      if ( $segment_slice_strand == -1 ) {
        $segment_offset = $segment->from_end();
      } else {
        $segment_offset = $segment->from_start();
      }
    }

    my $start = $slice_start + $segment_offset;
    my $end   = $start;

    if ( $segment_slice_strand == -1 ) {
      $start -= $entry->[ENTRY_SEQREGIONEND] - $segment_slice_start;
      $end   -= $entry->[ENTRY_SEQREGIONSTART] - $segment_slice_start;
    } else {    # Assumes '0' is really the positive strand.
      $start += $entry->[ENTRY_SEQREGIONSTART] - $segment_slice_start;
      $end   += $entry->[ENTRY_SEQREGIONEND] - $segment_slice_start;
    }

    $entry->[ENTRY_SEQREGIONSTART] = $start - 1;
    $entry->[ENTRY_SEQREGIONEND]   = $end - 1;

    push( @features, [ @{$entry} ] );
  } ## end while ( my $entry = $sth->fetchrow_arrayref...
  $sth->finish();

  return \@features;
} ## end sub _objs_from_sth

1;
