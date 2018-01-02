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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor

=head1 SYNOPSIS

  my $dfa = $database_adaptor->get_DensityFeatureAdaptor();

  my $interpolate   = 1;
  my $blocks_wanted = 50;

  @dense_feats = @{
    $dfa->fetch_all_by_Slice( $slice, 'SNPDensity', $blocks_wanted,
      $interpolate );
    }

=head1 DESCRIPTION

Density Feature Adaptor - An adaptor responsible for the creation of density
features from the database.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor;
use vars qw(@ISA);
use strict;


use POSIX;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Utils::Cache;
use Bio::EnsEMBL::DensityFeature;
use Bio::EnsEMBL::DensityFeatureSet;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);

our $DENSITY_FEATURE_CACHE_SIZE = 20;


=head2 new

  Arg [1]    : list of args @args
               Superclass constructor arguments
  Example    : none
  Description: Constructor which just initializes internal cache structures
  Returntype : Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor
  Exceptions : none
  Caller     : implementing subclass constructors
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  #initialize an LRU cache
  my %cache;
  tie(%cache, 'Bio::EnsEMBL::Utils::Cache', $DENSITY_FEATURE_CACHE_SIZE);
  $self->{'_density_feature_cache'} = \%cache;

  return $self;
}

=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice - The slice representing the region
               to retrieve density features from.
  Arg [2]    : string $logic_name - The logic name of the density features to
               retrieve.
  Arg [3]    : int $num_blocks (optional; default = 50) - The number of
               features that are desired. The ratio between the size of these
               features and the size of the features in the database will be
               used to determine which database features will be used.
  Arg [4]    : boolean $interpolate (optional; default = 0) - A flag indicating
               whether the features in the database should be interpolated to
               fit them to the requested number of features.  If true the
               features will be interpolated to provide $num_blocks features.
               This will not guarantee that exactly $num_blocks features are
               returned due to rounding etc. but it will be close.
  Arg [5]    : float $max_ratio - The maximum ratio between the size of the
               requested features (as determined by $num_blocks) and the actual
               size of the features in the database.  If this value is exceeded
               then an empty list will be returned.  This can be used to
               prevent excessive interpolation of the database values.
  Example    : #interpolate:
               $feats = $dfa->fetch_all_by_Slice($slice,'SNPDensity', 10, 1);
               #do not interpoloate, get what is in the database:
               $feats = $dfa->fetch_all_by_Slice($slice,'SNPDensity', 50);
               #interpolate, but not too much
               $feats = $dfa->fetch_all_by_Slice($slice,'SNPDensity',50,1,5.0);
  Description: Retrieves a set of density features which overlap the region
               of this slice. Density features are a discrete representation
               of a continuous value along a sequence, such as a density or
               percent coverage.  Density Features may be stored in chunks of
               different sizes in the database, and interpolated to fit the
               desired size for the requested region.  For example the database
               may store a single density value for each 1KB and also for each
               1MB.  When fetching for an entire chromosome the 1MB density
               chunks will be used if the requested number of blocks is not
               very high.
               Note that features which have been interpolated are not stored
               in the database and as such will have no dbID or adaptor set.
  Returntype : Bio::EnsEMBL::DensityFeature
  Exceptions : warning on invalid $num_blocks argument
               warning if no features with logic_name $logic_name exist
               warning if density_type table has invalid block_size value
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Slice {
  my ($self, $slice, $logic_name, $num_blocks, $interpolate, $max_ratio) = @_;

  if(defined($num_blocks) && $num_blocks < 1) {
    warning("Invalid number of density blocks [$num_blocks] requested.\n" .
	   "Returning empty list.");
    return [];
  }

  $num_blocks ||= 50;
  my $length = $slice->length();

  my $wanted_block_size = POSIX::ceil($length/$num_blocks);

  #
  # get out all of the density types and choose the one with the
  # block size closest to our desired block size
  #

  my $dta = $self->db()->get_DensityTypeAdaptor();

  my @dtypes = @{$dta->fetch_all_by_logic_name($logic_name)};
  if( ! @dtypes ){
    my @all_dtypes =  @{ $dta->fetch_all() };
    @all_dtypes or warning( "No DensityTypes in $dta" ) && return [];
    my $valid_list = join( ", ", map{$_->analysis->logic_name} @all_dtypes );
    warning( "No DensityTypes for logic name $logic_name. ".
	     "Select from $valid_list" );
    return [];
  }

  my $best_ratio   = undef;
  my $density_type = undef;
  my $best_ratio_large = undef;
  my $density_type_large = undef;
  my %dt_ratio_hash;

  foreach my $dt (@dtypes) {

    my $ratio;
    if( $dt->block_size() > 0 ) {
      $ratio = $wanted_block_size/$dt->block_size();
    } else {
      # This is only valid if the requested seq_region is the one the 
      # features are stored on. Please use sensibly. Or find better implementation.

      my $block_size = $slice->seq_region_length() / $dt->region_features();
      $ratio = $wanted_block_size / $block_size;
    }
    
    # we prefer to use a block size that's smaller than the required one
    # (better results on interpolation).
    # give larger bits a disadvantage and make them comparable
    if( $ratio < 1 ) {
      $ratio = 5/$ratio;
    }

    $dt_ratio_hash{ $ratio } = $dt;
  }

  $best_ratio = (sort {$a<=>$b} (keys %dt_ratio_hash))[0];
  
  #the ratio was not good enough, or this logic name was not in the
  #density type table, return an empty list
  if(!defined($best_ratio) ||
     (defined($max_ratio) && $best_ratio > $max_ratio)) {
    return [];
  }

  $density_type = $dt_ratio_hash{$best_ratio};

  my $constraint = "df.density_type_id = " . $density_type->dbID();

  my @features =
    @{$self->fetch_all_by_Slice_constraint($slice,$constraint)};

  return \@features if(!$interpolate);

  #interpolate the features into new features of a different size
  my @out;
  #sort the features on start position
  @features = sort( { $a->start() <=> $b->start() } @features );

  #resize the features that were returned
  my $start = 1;
  my $end   = $start+$wanted_block_size-1;

  my $value_type = $density_type->value_type();

  # create a new density type object for the interpolated features that
  # is not stored in the database
  $density_type = Bio::EnsEMBL::DensityType->new
    (-VALUE_TYPE => $value_type,
     -ANALYSIS   => $density_type->analysis(),
     -BLOCK_SIZE => $wanted_block_size);

  while($start < $length) {
#    $end = $length if($end > $length);

    my $density_value = 0.0;
    my ($f, $fstart, $fend, $portion);
    my @dvalues;

    #construct a new feature using all of the old density features that
    #we overlapped
    while(($f = shift(@features)) && $end > $f->{'start'}) {

      #what portion of this feature are we using to construct our new block?
      $fstart = ($f->{'start'} < $start) ? $start : $f->{'start'};
      $fend   = ($f->{'end'}   > $end  ) ? $end   : $f->{'end'};
      $fend   = $length if($fend > $length);

      if($value_type eq 'sum') {

        $portion = ($fend-$fstart+1)/$f->length();

        #take a percentage of density value, depending on how much of the
        #feature we overlapped
        $density_value += $portion * $f->{'density_value'};

      } elsif($value_type eq 'ratio') {

        #maintain a running total of the length used from each feature
        #and its value
        push(@dvalues,[$fend-$fstart+1,$f->{'density_value'}]);

      } else {
        throw("Unknown density value type [$value_type].");
      }

      #do not want to look at next feature if we only used part of this one:
      last if($fend < $f->{'end'});
    }

    #if we did not completely overlap the last feature, put it back on so
    #it can be partially used by the next block
    if(defined($f) && (!defined($fend) || $fend < $f->{'end'})) {
      unshift(@features, $f);
    }

    if($value_type eq 'ratio') {
      #take a weighted average of the all the density values of the features
      #used to construct this one
      my $total_len = $end - $start + 1;
      if($total_len > 0) {
        foreach my $pair (@dvalues) {
          my ($dlen, $dval) = @$pair;
          $density_value += $dval * ($dlen/$total_len);
        }
      }
    }

    # interpolated features are not stored in the db so do not set the dbID
    # or the adaptor
    push @out, Bio::EnsEMBL::DensityFeature->new
      (-seq_region    => $slice,
       -start         => $start,
       -end           => $end,
       -density_type  => $density_type,
       -density_value => $density_value);

    $start = $end + 1;
    $end  += $wanted_block_size;
  }

  return \@out;
}


sub fetch_all {
  my $self = shift;
  my $logic_name = shift;

  my ($sth, $seq_region_id, $start, $end, $density_value, $density_type_id);
  my ($slice, $density_type, @out);
  my $sa = $self->db()->get_SliceAdaptor();
  my $dta = $self->db()->get_DensityTypeAdaptor();

  if ($logic_name) {
    $sth = $self->prepare("SELECT df.seq_region_id, df.seq_region_start, df.seq_region_end, df.density_value, df.density_type_id "
                         . "FROM density_feature df, density_type dt, analysis a "
                         . "WHERE df.density_type_id = dt.density_type_id "
                         . "AND dt.analysis_id = a.analysis_id "
                         . "AND logic_name = ?" );
    $sth->bind_param(1, $logic_name, SQL_VARCHAR);
    $sth->execute();
  } else {
    $sth = $self->prepare("SELECT df.seq_region_id, df.seq_region_start, df.seq_region_end, df.density_value, df.density_type_id "
                         . "FROM density_feature df, density_type dt "
                         . "WHERE df.density_type_id = dt.density_type_id" );
    $sth->execute();
  }


  $sth->bind_columns(\($seq_region_id, $start, $end, $density_value, $density_type_id));


  while ($sth->fetch()) {
    $slice = $sa->fetch_by_seq_region_id($seq_region_id);
    $density_type = $dta->fetch_by_dbID($density_type_id);
    push (@out, Bio::EnsEMBL::DensityFeature->new
      (-seq_region    => $slice,
       -start         => $start,
       -end           => $end,
       -density_type  => $density_type,
       -density_value => $density_value));
  }
  return \@out;
}


=head2 fetch_all_by_Slice_constraint

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Arg [2]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Example    : $fs = $a->fetch_all_by_Slice_constraint($slc, 'density_type_id = 88');
  Description: Returns a listref of features created from the database which 
               are on the Slice defined by $slice and fulfill the SQL 
               constraint defined by $constraint. 
               Note that this is a re-implementation of a method with the same name
               in the BaseFeatureAdaptor.
               This was necessary to remove the use of symliked sequences for density features
  Returntype : listref of Bio::EnsEMBL::SeqFeatures in Slice coordinates
  Exceptions : thrown if $slice is not defined
  Caller     : Bio::EnsEMBL::Slice
  Status     : Stable

=cut


sub fetch_all_by_Slice_constraint {
  my $self = shift;
  my $slice = shift;
  my $constraint = shift;

  my @result;

  #Pull out here as we need to remember them & reset accordingly
  my $bind_params = $self->bind_param_generic_fetch();

  if ( !ref($slice)
       || !(    $slice->isa('Bio::EnsEMBL::Slice')
             or $slice->isa('Bio::EnsEMBL::LRGSlice') ) )
  {
    throw("Bio::EnsEMBL::Slice argument expected.");
  }

  $constraint ||= '';

  # If the logic name was invalid, undef was returned
  if ( !defined($constraint) ) { return [] }

  my $key;
  my $cache;

  # Will only use feature_cache if hasn't been no_cache attribute set
  if (
    !( defined( $self->db()->no_cache() ) && $self->db()->no_cache() ) )
  {

    #strain test and add to constraint if so to stop caching.
    if ( $slice->isa('Bio::EnsEMBL::StrainSlice') ) {
      my $string =
        $self->dbc()->db_handle()->quote( $slice->strain_name() );

      if ( $constraint ne "" ) {
        $constraint .= " AND $string = $string ";
      } else {
        $constraint .= " $string = $string ";
      }
    }

    # Check the cache and return the cached results if we have already
    # done this query.  The cache key is the made up from the slice
    # name, the constraint, and the bound parameters (if there are any).
    $key = uc( join( ':', $slice->name(), $constraint ) );

    if ( defined($bind_params) ) {
      $key .= ':'
        . join( ':', map { $_->[0] . '/' . $_->[1] } @{$bind_params} );
    }

    $cache = $self->_slice_feature_cache();
    if ( exists( $cache->{$key} ) ) {
      # Clear the bound parameters and return the cached data.
      $self->{'_bind_param_generic_fetch'} = ();
      return $cache->{$key};
    }
  } ## end if ( !( defined( $self...)))

  $self->_bind_param_generic_fetch($bind_params);
  my $features = $self->_slice_fetch( $slice, $constraint );

  if ( defined($key) ) {
    $cache->{$key} = $features;
  }

  return $features;
}


sub _tables {
  my $self = shift;

  return (['density_feature', 'df']);
}


sub _columns {
  my $self = shift;

  return qw( df.density_feature_id
             df.seq_region_id df.seq_region_start df.seq_region_end
             df.density_value df.density_type_id );
}



sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;

  #
  # This code is ugly because an attempt has been made to remove as many
  # function calls as possible for speed purposes.  Thus many caches,
  # and a fair bit of gymnastics have been used.
  #

  my $sa  = $self->db()->get_SliceAdaptor();
  my $dta = $self->db()->get_DensityTypeAdaptor();

  my @features;
  my %dtype_hash;
  my %slice_hash;
  my %sr_name_hash;
  my %sr_cs_hash;

  my(
    $density_feature_id, $seq_region_id, $seq_region_start,
    $seq_region_end, $density_value, $density_type_id );

  $sth->bind_columns(\(
                      $density_feature_id, $seq_region_id, $seq_region_start,
                      $seq_region_end, $density_value, $density_type_id));

  my $dest_slice_start;
  my $dest_slice_end;
  my $dest_slice_strand;
  my $dest_slice_length;
  my $dest_slice_cs;
  my $dest_slice_sr_name;
  my $dest_slice_sr_id;
  my $asma;

  if ($dest_slice) {
    $dest_slice_start   = $dest_slice->start();
    $dest_slice_end     = $dest_slice->end();
    $dest_slice_strand  = $dest_slice->strand();
    $dest_slice_length  = $dest_slice->length();
    $dest_slice_cs      = $dest_slice->coord_system();
    $dest_slice_sr_name = $dest_slice->seq_region_name();
    $dest_slice_sr_id   = $dest_slice->get_seq_region_id();
    $asma               = $self->db->get_AssemblyMapperAdaptor();
  }

  FEATURE: while($sth->fetch()) {

    #get the density type object
    my $density_type = $dtype_hash{$density_type_id} ||= $dta->fetch_by_dbID($density_type_id);

    #need to get the internal_seq_region, if present
    $seq_region_id = $self->get_seq_region_id_internal($seq_region_id);
    my $slice = $slice_hash{"ID:".$seq_region_id};

    if (!$slice) {
      $slice                            = $sa->fetch_by_seq_region_id($seq_region_id);
      $slice_hash{"ID:".$seq_region_id} = $slice;
      $sr_name_hash{$seq_region_id}     = $slice->seq_region_name();
      $sr_cs_hash{$seq_region_id}       = $slice->coord_system();
    }

    #obtain a mapper if none was defined, but a dest_seq_region was
    if(!$mapper && $dest_slice && !$dest_slice_cs->equals($slice->coord_system)) {
      $mapper = $asma->fetch_by_CoordSystems($dest_slice_cs, $slice->coord_system);
    }

    my $sr_name = $sr_name_hash{$seq_region_id};
    my $sr_cs   = $sr_cs_hash{$seq_region_id};

    #
    # remap the feature coordinates to another coord system
    # if a mapper was provided
    #

    if ($mapper) {

      my $len = $seq_region_end - $seq_region_start + 1;

      if (defined $dest_slice && $mapper->isa('Bio::EnsEMBL::ChainedAssemblyMapper') ) {
        ($seq_region_id, $seq_region_start, $seq_region_end) =
         $mapper->map( $sr_name, $seq_region_start, $seq_region_end, 1, $sr_cs, 0, $dest_slice);

      } else {
        ($seq_region_id, $seq_region_start, $seq_region_end) =
         $mapper->map($sr_name, $seq_region_start, $seq_region_end, 1, $sr_cs);
      }

      #skip features that map to gaps or coord system boundaries
      next FEATURE if (!defined($seq_region_id));

      if($density_type->value_type() eq 'sum') {
        #adjust density value so it reflects length of feature actually used
        my $newlen = $seq_region_end - $seq_region_start +1;
        $density_value *= $newlen/$len if($newlen != $len);
      }

      #get a slice in the coord system we just mapped to
      $slice = $slice_hash{"ID:".$seq_region_id} ||= $sa->fetch_by_seq_region_id($seq_region_id);
    }

    #
    # If a destination slice was provided convert the coords.
    #
    if (defined($dest_slice)) {
      my $seq_region_len = $dest_slice->seq_region_length();

      if ( $dest_slice_strand == 1 ) {
        $seq_region_start = $seq_region_start - $dest_slice_start + 1;
        $seq_region_end   = $seq_region_end - $dest_slice_start + 1;

        if ( $dest_slice->is_circular ) {
        # Handle circular chromosomes.

          if ( $seq_region_start > $seq_region_end ) {
            # Looking at a feature overlapping the chromosome origin.

            if ( $seq_region_end > $dest_slice_start ) {
              # Looking at the region in the beginning of the chromosome
              $seq_region_start -= $seq_region_len;
            }
            if ( $seq_region_end < 0 ) {
              $seq_region_end += $seq_region_len;
            }
          } else {
            if ($dest_slice_start > $dest_slice_end && $seq_region_end < 0) {
              # Looking at the region overlapping the chromosome
              # origin and a feature which is at the beginning of the
              # chromosome.
              $seq_region_start += $seq_region_len;
              $seq_region_end   += $seq_region_len;
            }
          }
        }
      } else {

        my $start = $dest_slice_end - $seq_region_end + 1;
        my $end = $dest_slice_end - $seq_region_start + 1;

        if ($dest_slice->is_circular()) {

          if ($dest_slice_start > $dest_slice_end) {
            # slice spans origin or replication

            if ($seq_region_start >= $dest_slice_start) {
              $end += $seq_region_len;
              $start += $seq_region_len if $seq_region_end > $dest_slice_start;

            } elsif ($seq_region_start <= $dest_slice_end) {
              # do nothing
            } elsif ($seq_region_end >= $dest_slice_start) {
              $start += $seq_region_len;
              $end += $seq_region_len;

            } elsif ($seq_region_end <= $dest_slice_end) {
              $end += $seq_region_len if $end < 0;

            } elsif ($seq_region_start > $seq_region_end) {
              $end += $seq_region_len;
            }

          } else {

            if ($seq_region_start <= $dest_slice_end and $seq_region_end >= $dest_slice_start) {
              # do nothing
            } elsif ($seq_region_start > $seq_region_end) {
              if ($seq_region_start <= $dest_slice_end) {
                $start -= $seq_region_len;
              } elsif ($seq_region_end >= $dest_slice_start) {
                $end += $seq_region_len;
              }
            }
          }
        }

        $seq_region_start = $start;
        $seq_region_end = $end;
      }

      # Throw away features off the end of the requested slice or on
      # different seq_region.
      if ($seq_region_end < 1
          || $seq_region_start > $dest_slice_length
          || ($dest_slice_sr_id != $seq_region_id)) {
        next FEATURE;
      }
      $slice = $dest_slice;
    }

    push(
      @features,
      $self->_create_feature(
      'Bio::EnsEMBL::DensityFeature', {
       -dbID       => $density_feature_id,
       -adaptor    => $self,
       -start      => $seq_region_start,
       -end        => $seq_region_end,
       -seq_region => $slice,
       -density_value => $density_value,
       -density_type  => $density_type
     } ) );

  }

  return \@features;
}


=head2 list_dbIDs

  Arg [1]    : none
  Example    : @feature_ids = @{$density_feature_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all density features in the
               current db
  Arg[1]     : <optional> int. not set to 0 for the ids to be sorted by the seq_region.
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub list_dbIDs {
   my ($self, $ordered) = @_;

   return $self->_list_dbIDs("density_feature",undef, $ordered);
}


=head2 store

  Arg [1]    : list of Bio::EnsEMBL::DensityFeatures @df
               the simple features to store in the database
  Example    : $density_feature_adaptor->store(1234, @density_feats);
  Description: Stores a list of density feature objects in the database
  Returntype : none
  Exceptions : thrown if @df is not defined, if any of the features do not
               have an attached slice.
               or if any elements of @df are not Bio::EnsEMBL::SeqFeatures 
  Caller     : general
  Status     : Stable

=cut

sub store{
  my ($self,@df) = @_;

  if( scalar(@df) == 0 ) {
    throw("Must call store with list of DensityFeatures");
  }
#mysql> desc density_feature;
#+--------------------+---------+------+-----+---------+----------------+
#| Field              | Type    | Null | Key | Default | Extra          |
#+--------------------+---------+------+-----+---------+----------------+
#| density_feature_id | int(11) |      | PRI | NULL    | auto_increment |
#| density_type_id    | int(11) |      | MUL | 0       |                |
#| seq_region_id      | int(11) |      |     | 0       |                |
#| seq_region_start   | int(11) |      |     | 0       |                |
#| seq_region_end     | int(11) |      |     | 0       |                |
#| density_value      | float   |      |     | 0       |                |
#+--------------------+---------+------+-----+---------+----------------+

  my $sth = $self->prepare
    ("INSERT INTO density_feature (seq_region_id, seq_region_start, " .
                                  "seq_region_end, density_type_id, " .
                                  "density_value) " .
     "VALUES (?,?,?,?,?)");

  my $db = $self->db();
  my $analysis_adaptor = $db->get_AnalysisAdaptor();

 FEATURE: foreach my $df ( @df ) {

    if( !ref $df || !$df->isa("Bio::EnsEMBL::DensityFeature") ) {
      throw("DensityFeature must be an Ensembl DensityFeature, " .
            "not a [".ref($df)."]");
    }
    
    # we dont store 0 value density features
    next if( $df->density_value == 0 );
    if($df->is_stored($db)) {
      warning("DensityFeature [".$df->dbID."] is already stored" .
              " in this database.");
      next FEATURE;
    }

    if(!defined($df->density_type)) {
      throw("A density type must be attached to the features to be stored.");
    }

    #store the density_type if it has not been stored yet

    if(!$df->density_type->is_stored($db)) {
      my $dta = $db->get_DensityTypeAdaptor();
      $dta->store($df->density_type());
    }

    my $original = $df;
    my $seq_region_id;
    ($df, $seq_region_id) = $self->_pre_store($df);

    $sth->bind_param(1,$seq_region_id,SQL_INTEGER);
    $sth->bind_param(2,$df->start,SQL_INTEGER);
    $sth->bind_param(3,$df->end,SQL_INTEGER);
    $sth->bind_param(4,$df->density_type->dbID,SQL_INTEGER);
    $sth->bind_param(5,$df->density_value,SQL_FLOAT);
    $sth->execute();

    $original->dbID($self->last_insert_id('density_feature_id', undef, 'density_feature'));
    $original->adaptor($self);
  }
}

=head2 fetch_Featureset_by_Slice 

  Arg [1-5]  : see
               Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor::fetch_all_by_Slice()
               for argument documentation
  Example    : $featureset = $dfa->fetch_FeatureSet_by_Slice($slice,'SNPDensity', 10, 1);
  Description: wrapper around
               Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor::fetch_all_by_Slice()
               which returns a Bio::EnsEMBL::DensityFeatureSet and also caches
               results
  Returntype : Bio::EnsEMBL::DensityFeatureSet
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_Featureset_by_Slice {
  my ($self, $slice, $logic_name, $num_blocks, $interpolate, $max_ratio) = @_;

  my $key = join(":", $slice->name,
                      $logic_name,
                      $num_blocks || 50,
                      $interpolate || 0,
                      $max_ratio);
                      
  unless ($self->{'_density_feature_cache'}->{$key}) {
    my $dfeats = $self->fetch_all_by_Slice($slice, $logic_name, $num_blocks,
                                           $interpolate, $max_ratio);
    $self->{'_density_feature_cache'}->{$key} =
      new Bio::EnsEMBL::DensityFeatureSet($dfeats);
  }
  return $self->{'_density_feature_cache'}->{$key};
}


1;
