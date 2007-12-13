# $Id$

package Bio::EnsEMBL::Collection;

use strict;
use warnings;

use Carp;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use base qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );

# Symbolic constants that acts as indices into the individual feature
# collection entries.  These must be in sync with the columns returned
# by the _columns() method.
use constant { ENTRY_DBID            => 0,
               ENTRY_SEQREGIONID     => 1,
               ENTRY_SEQREGIONSTART  => 2,
               ENTRY_SEQREGIONEND    => 3,
               ENTRY_SEQREGIONSTRAND => 4 };

#-----------------------------------------------------------------------
# Constructor
#-----------------------------------------------------------------------

sub new {
  my ( $proto, $slice, $light ) = @_;

  my $this = $proto->SUPER::new( $slice->adaptor()->db() );

  $this->__attrib( 'is_light', $light );

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
                                    'build_level'        => $build_level
    };
  }
  $sth->finish();

  $this->__attrib( 'coordinate_systems', \%coordinate_systems );

  $this->slice($slice);

  return $this;
} ## end sub new

#-----------------------------------------------------------------------
# PUBLIC methods
#-----------------------------------------------------------------------

# Returns true if the collection is lightweight, i.e. if its entries
# does not contain any feature-specific data (e.g. transcript and gene
# stable IDs for a transcript feature collection).
sub is_light {
  my ($this) = @_;
  return $this->__attrib('is_light');
}

# Getter/setter method for the main slice associated with this
# collection.  Associating a new slice with the collection will empty
# all entries from the collection.
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

    $this->collection( [] );
    $this->__attrib( 'is_populated', 0 );

    $this->__attrib( 'slice', $slice );
  } ## end if ( defined($slice) )

  return $this->__attrib('slice');
} ## end sub slice

# Getter/setter for the collection array.
sub collection {
  my ( $this, $collection ) = @_;
  return $this->__attrib( 'collection', $collection );
}

# Returns true if the collection has been populated, false of not.
sub is_populated {
  my ($this) = @_;
  return $this->__attrib('is_populated');
}

# Populates the collection with a compact representation of the features
# overlapping the current slice.
sub populate {
  my ($this) = @_;

  if ( $this->is_populated() ) { return }

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

  $this->__attrib(
        'collection', [
          sort( {
              $a->[ENTRY_SEQREGIONSTART] <=> $b->[ENTRY_SEQREGIONSTART]
                || $a->[ENTRY_SEQREGIONEND] <=> $b->[ENTRY_SEQREGIONEND]
            } @entries ) ] );

  $this->__attrib( 'is_populated', 1 );
} ## end sub populate

# Counts the number of features that overlaps the current slice.
sub count {
  my ($this) = @_;

  $this->populate();

  return scalar( @{ $this->collection() } );
}

#
# Binning methods
#
# Each binning method bins the collection entries and returns an array
# with the specified length (number of bins).  An entry, which basically
# consists of a start and a end position, is allocated to one or several
# bins depending on its span and the size of the individual bins.
#
# Features:     |----|     |----------------|  |--|         |------|
#               |-------------|                |-----|             |--|
#               |-------------------|                 |----|       |--|
#
# Fine bins:    3 3 3 2 2 2 3 3 2 2 2 1 1 1 1 0 2 2 1 1 1 1 1 1 1 1 2 2
# Coarse bins:  3  3  2  2  3  2  2  1  1  1  0  2  1  1  1  1  1  3  2
#
# The example above shows the arrays returned from get_bin_counts(28)
# and get_bin_counts(19).

# Returns an array of bins, each bin containing the number of entries
# allocated to that bin.  This is equivalent to the density of the
# entries.
sub get_bin_counts {
  my ( $this, $nbins ) = @_;

  return $this->__bin(
    $nbins,
    sub {
      my ( $bins, $bin_index, $entry, $entry_index ) = @_;
      ++$bins->[$bin_index];
    } );
}

# Returns an array of bins, each bin containing an array of indices into
# the collection array for the entries allocated to that bin.
sub get_bin_entry_indices {
  my ( $this, $nbins ) = @_;

  return $this->__bin(
    $nbins,
    sub {
      my ( $bins, $bin_index, $entry, $entry_index ) = @_;
      push( @{ $bins->[$bin_index] }, $entry_index );
    } );
}

# Returns an array of bins, each bin containing an array of entries
# (array references into the collection array) allocated to that bin.
sub get_bin_entries {
  my ( $this, $nbins ) = @_;

  return $this->__bin(
    $nbins,
    sub {
      my ( $bins, $bin_index, $entry, $entry_index ) = @_;
      push( @{ $bins->[$bin_index] }, $entry );
    } );
}

#-----------------------------------------------------------------------
# PRIVATE methods
#-----------------------------------------------------------------------

sub __attrib {
  my ( $this, $attribute, $value ) = @_;

  if ( defined($value) ) {
    $this->{'attributes'}{$attribute} = $value;
  }

  return $this->{'attributes'}{$attribute};
}

sub __bin {
  my ( $this, $nbins, $bin_entry_sub ) = @_;

  $this->populate();

  my $slice       = $this->slice();
  my $slice_start = $slice->start();

  my $bin_length = $slice->length()/$nbins;

  my @bins;
  my $entry_index = 0;

  foreach my $entry ( @{ $this->collection() } ) {
    my $start_bin = int(
        ( $entry->[ENTRY_SEQREGIONSTART] - $slice_start )/$bin_length );
    my $end_bin = int(
          ( $entry->[ENTRY_SEQREGIONEND] - $slice_start )/$bin_length );

    if ( $end_bin >= $nbins ) { $end_bin = $nbins - 1 }

    for ( my $bin_index = $start_bin ;
          $bin_index <= $end_bin ;
          ++$bin_index )
    {
      $bin_entry_sub->( \@bins, $bin_index, $entry, $entry_index );
    }

    ++$entry_index;
  }

  return \@bins;
} ## end sub __bin

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
# PROTECTED methods
#-----------------------------------------------------------------------

# _extra_tables() provides a list of tables and table aliases, in
# addition to the ones returned by _tables(), that are used to
# create collections of particular feature types, e.g. the tables
# transcript_stable_id [tsi], gene [g], and gene_stable_id [gsi] for
# Bio::EnsEMBL::Collection::Transcript.
sub _extra_tables { return () }

# _extra_columns() provides a list of columns, in addition to the ones
# returned by _columns(), that constitutes the contents of an entry of a
# particular feature type, e.g. t.biotype, t.status, tsi.stable_id, and
# gsi.stable_id for Bio::EnsEMBL::Collection::Transcript.
sub _extra_columns { return () }

# _extra_where_clause() should join the tables from _extra_tables().
sub _extra_where_clause { return undef }

# _dbID_column() provides the name of the column in the primary
# feature table that holds the dbID for the collection elements, i.e.
# 'transcript_id' for Bio::EnsEMBL::Collection::Transcript.
#
# If this method is not specialized by a subclass, it is assumed that
# the name of the dbID column is the name of the primary feature table
# suffixed by '_id'.
sub _dbID_column {
  my ($this) = @_;

  if ( !defined( $this->__attrib('dbID_column') ) ) {
    $this->__attrib( 'dbID_column',
                     $this->_feature_table()->[0] . '_id' );
  }

  return $this->__attrib('dbID_column');
}

#-----------------------------------------------------------------------
# ABSTRACT PROTECTED methods
#-----------------------------------------------------------------------

# _feature_table() should return the primary feature table and table
# alias used by the collection elements, i.e. the transcript [t] table
# for Bio::EnsEMBL::Collection::Transcript.
sub _feature_table {
  croak("Called abstract method '_feature_table()'");
}

#-----------------------------------------------------------------------
# Implemented abstract protected methods from base class
# Bio::EnsEMBL::BaseAdaptor
#-----------------------------------------------------------------------

sub _tables {
  my ($this) = @_;

  my @tables = ( $this->_feature_table() );

  if ( !$this->is_light() ) {
    push( @tables, $this->_extra_tables() );
  }

  return @tables;
}

sub _columns {
  my ($this) = @_;

  my $table_alias = $this->_feature_table()->[1];

  my @columns = ( $table_alias . '.' . $this->_dbID_column(),
                  $table_alias . '.seq_region_id',
                  $table_alias . '.seq_region_start',
                  $table_alias . '.seq_region_end',
                  $table_alias . '.seq_region_strand' );

  if ( !$this->is_light() ) {
    push( @columns, $this->_extra_columns() );
  }

  return @columns;
}

sub _default_where_clause {
  my ($this) = @_;

  if ( !$this->is_light() ) {
    my $extra_where = $this->_extra_where_clause();
    if ( defined($extra_where) ) {
      return $extra_where;
    }
  }

  return '';
}

sub _objs_from_sth {
  my ( $this, $sth ) = @_;

  my $seqreg_map = $this->__attrib('seqreg_map');

  my @features;

  my ( $segment, $segment_slice, $segment_slice_start,
       $segment_slice_strand, $segment_offset );

  my $slice_start = $this->slice()->start();

  while ( my $entry = $sth->fetchrow_arrayref() ) {
    if ( !defined($segment)
         || $segment != $seqreg_map->{ $entry->[1] } )
    {
      $segment              = $seqreg_map->{ $entry->[1] };
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
