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

  # Create a feature collection on the slice.
  my $collection =
    Bio::EnsEMBL::Collection::RepeatFeature->new( -slice => $slice,
                                                  -light => 0 );

  # Populate the feature collection from the slice.
  $collection->populate();

  # Populate the feature collection from the slice.  Sort the
  # entries on feature start position.
  $collection->populate( -sorted => 1 );

  # Populate the feature collection from the slice.  Sort the
  # entries on feature length and don't bother about feature
  # type-specific data.
  $collection->populate(
    -light   => 1,
    -sorted  => 1,
    -sortsub => sub {
      $a->[3] - $a->[2] <=> $b->[3] - $b->[2];
    } );

  # Retrieve the entries from the collection as an array of arrays.
  my @entries = @{ $collection->entries() };

  # Retrieve a binned representation of the entries using 100 bins.
  my @binned_entries = @{ $collection->get_bin_entries(100) };

  # Retrieve a binned representation of the entries using 100 bins,
  # where each entry is represented by its index in the feature
  # collection array (@entries above).
  my @binned_entry_indicies =
    @{ $collection->get_bin_entry_indices(100) };

  # Retrieve only the bin counts.
  my @bin_counts = @{ $collection->get_bin_counts(100) };

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

A sub-class of this abstract base class will specify further data to be
added to the entries in order to account of the particular feature type.
An entry from a gene feature collection (Bio::EnsEMBL::Collection::Gene)
might, for example, contain the Ensembl Stable ID of the gene.

=head2 Light-weight collections/entries

A light-weight collection is a feature collection whose collection
entries are light-weight, whose entries does not contain the extended
feature representation.

A collection which is light-weight by default may be created by using
the argument '-light=>1' in the constructor, new(), but one may also use
the same argument with populate().

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

# Symbolic constants that acts as indices into the individual feature
# collection entries.  These must be in sync with the columns returned
# by the _columns() method.
use constant { ENTRY_DBID            => 0,
               ENTRY_SEQREGIONID     => 1,
               ENTRY_SEQREGIONSTART  => 2,
               ENTRY_SEQREGIONEND    => 3,
               ENTRY_SEQREGIONSTRAND => 4 };

=head1 METHODS (constructor)

=head2 new

  Arg [SLICE]   : Bio::EnsEMBL::Slice
                  The slice to be associated with this feature
                  collection.

  Arg [LIGHT]   : Boolean (optional, default false)
                  If true, the collection will be 'light-weight',
                  i.e. no type-specific data will be stored in its
                  entries when populate() is called.

  Example       : my $collection =
                    Bio::EnsEMBL::Collection::<feature_type>->new(
                                                       -slice => $slice,
                                                       -light => 1 );

  Description   : When called for a sub-class, creates a feature
                  collection object and associates a slice with it.

  Return type   : Bio::EnsEMBL::Collection::<feature_type>

  Exceptions    : Throws if no slice is specified.

  Caller        : General (through a sub-class).

  Status        : At Risk (under development)

=cut

sub new {
  my $proto = shift;
  my ( $slice, $light ) = rearrange( [ 'SLICE', 'LIGHT' ], @_ );

  if ( !defined($slice) ) {
    throw(   'Unable to create a feature collection '
           . 'without a slice object, I am.' );
  }

  my $this = $proto->SUPER::new( $slice->adaptor()->db() );

  my $sql = qq(
    SELECT  cs.name, mc.max_length, m.meta_value
    FROM    coord_system cs,
            meta_coord mc
    LEFT JOIN meta m ON m.meta_key = ?
    WHERE   mc.table_name = ?
    AND     mc.coord_system_id = cs.coord_system_id
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

  $this->lightweight($light);
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
    my %seqreg_map;
    my @all_segments;

    foreach
      my $cs_name ( keys( %{ $this->__attrib('coordinate_systems') } ) )
    {
      my @segments = @{ $slice->project($cs_name) };

      foreach my $segment (@segments) {
        $seqreg_map{ $segment->to_Slice()->get_seq_region_id() } =
          $segment;
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

      push( @all_segments, @segments );
    }

    # The 'segments' list is a list of ProjectionSegment objects and of
    # arrays of ProjectionSegment objects.  For segments in arrays, no
    # constraint on seq_region_start or seq_region_end is needed.
    $this->__attrib( 'segments', \@all_segments );

    # This is a simple map between seq_region_id and a
    # ProjectionSegment, used in the mapping done by the
    # _objs_from_sth() method.
    $this->__attrib( 'seqreg_map', \%seqreg_map );

    $this->entries( [] );
    $this->__attrib( 'is_populated', 0 );

    $this->__attrib( 'slice', $slice );
  } ## end if ( defined($slice) )

  return $this->__attrib('slice');
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

      [LIGHT]   : Boolean (optional, default undef)
                  If true, will populate the collection with
                  light-weight entries.  If false, will populate the
                  collection with entries that are not light-weight.
                  If unset, will use the 'lightweight' boolean set
                  when calling new().

  Example       : $collection->populate( -sorted => 1, -light => 0 );

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
  my ( $sorted, $sort_like_this, $light ) =
    rearrange( [ 'SORTED', 'SORTSUB', 'LIGHT' ], @_ );

  if ( $this->is_populated() ) { return }

  # Save the old lightweight() value so that we can restore it if the
  # -light argument was used.

  my $old_light = $this->lightweight();
  if ( defined($light) ) {
    $this->lightweight($light);
  }

  my @entries;

  foreach my $segment ( @{ $this->__attrib('segments') } ) {
    if ( defined($segment)
         && ( ref($segment) ne 'ARRAY'
              || scalar( @{$segment} ) > 0 ) )
    {
      push( @entries,
            @{ $this->generic_fetch( $this->__constraint($segment) ) }
      );
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

  # Restore the lightweight() value if -light was used.
  if ( defined($light) ) {
    $this->lightweight($old_light);
  }
} ## end sub populate

=head2 lightweight

  Arg [1]       : Boolean (optional)

  Example       : if ( !$collection->lightweight() ) { ... }

  Description   : Getter/setter for the 'lightweight' boolean.
                  If the collection is light-weight, its entries
                  does not contain any feature-specific data (e.g.
                  transcript or gene stable IDs for a transcript
                  feature collection).

  Return type   : Boolean

  Exceptions    : None

  Caller        : General

  Status        : At Risk (under development)

=cut

sub lightweight {
  my ( $this, $light ) = @_;
  return $this->__attrib( 'lightweight', $light );
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

  Caller        : General

  Status        : At Risk (under development)

=cut

our %valid_binning_methods = ( 'count'    => 0,
                               'density'  => 0,    # Same as 'count'.
                               'indices'  => 1,
                               'index'    => 1,    # Same as 'indices'.
                               'entries'  => 2,
                               'entry'    => 2,    # Same as 'entries'.
                               'coverage' => 3 );

sub get_bins {
  my $this = shift;
  my ( $nbins, $method_name ) = rearrange( [ 'NBINS', 'METHOD' ], @_ );

  if ( !$this->is_populated() ) {
    throw(   'Can not bin a feature collection '
           . 'without first having called populate()' );
  }

  $method_name ||= 'count';
  if ( !exists( $valid_binning_methods{$method_name} ) ) {
    throw(
           sprintf(
                "Invalid binning method '%s', valid methods are: %s",
                $method_name, join( ', ', keys(%valid_binning_methods) )
           ) );
  }
  my $method = $valid_binning_methods{$method_name};

  my $slice       = $this->slice();
  my $slice_start = $slice->start();

  my $bin_length = ( $slice->end() - $slice_start + 1 )/$nbins;

  my @bins = map( $_ = undef, 0 .. $nbins - 1 );
  my $entry_index = 0;

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

      # For 'coverage'.

      throw("The 'coverage' binning method is not yet implemented");

      # FIXME

    }

  } ## end foreach my $entry ( @{ $this...

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

  Description   : The provides the name of the column in the
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

  if ( !$this->lightweight() ) {
    push( @tables, $this->_extra_tables() );
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

  if ( !$this->lightweight() ) {
    push( @columns, $this->_extra_columns() );
  }

  return @columns;
}

=head2 _default_where_clause

  Args          : None

  Description   : For a light-weight collection, returns an
                  empty string (no joins), otherwise returns the
                  WHERE clause that joins the tables returned by
                  _tables().

  Return type   : String

  Exceptions    : None

  Caller        : Bio::EnsEMBL::DBSQL::BaseAdaptor::generic_fetch()

  Status        : At Risk (under development)

=cut

sub _default_where_clause {
  my ($this) = @_;

  if ( !$this->lightweight() ) {
    my $extra_where = $this->_extra_where_clause();
    if ( defined($extra_where) ) {
      return $extra_where;
    }
  }

  return '';
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

  my $seqreg_map = $this->__attrib('seqreg_map');

  my @features;

  my ( $segment, $segment_slice, $segment_slice_start,
       $segment_slice_strand, $segment_offset );

  my $slice_start = $this->slice()->start();

  while ( my $entry = $sth->fetchrow_arrayref() ) {
    if ( !defined($segment)
         || $segment != $seqreg_map->{ $entry->[ENTRY_SEQREGIONID] } )
    {
      $segment       = $seqreg_map->{ $entry->[ENTRY_SEQREGIONID] };
      $segment_slice = $segment->to_Slice();
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
