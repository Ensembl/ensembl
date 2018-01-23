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

Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor - An Abstract Base class for all
FeatureAdaptors

=head1 SYNOPSIS

Abstract class - should not be instantiated.  Implementation of
abstract methods must be performed by subclasses.

=head1 DESCRIPTION

This is a base adaptor for feature adaptors. This base class is simply a way
of eliminating code duplication through the implementation of methods
common to all feature adaptors.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use vars qw(@ISA @EXPORT);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Cache;
use Bio::EnsEMBL::Utils::Exception qw(warning throw deprecate stack_trace_dump);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Iterator;
use Bio::EnsEMBL::Utils::Scalar qw/assert_ref/;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});

our $SLICE_FEATURE_CACHE_SIZE    = 4;
our $MAX_SPLIT_QUERY_SEQ_REGIONS = 3;
our $SILENCE_CACHE_WARNINGS = 0;

=head2 new

  Arg [1]    : list of args @args
               Superclass constructor arguments
  Example    : none
  Description: Constructor which warns if caching has been switched off
  Returntype : Bio::EnsEMBL::BaseFeatureAdaptor
  Exceptions : none
  Caller     : implementing subclass constructors
  Status     : Stable

=cut

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  if ( defined $self->db->no_cache() && $self->db->no_cache() && ! $SILENCE_CACHE_WARNINGS) {
    warning(  "You are using the API without caching most recent features. "
            . "Performance might be affected." );
  }
  return $self;
}

=head2 start_equals_end

  Arg [1]    : (optional) boolean $newval
  Example    : $bfa->start_equals_end(1);
  Description: Getter/Setter for the start_equals_end flag.  If set
               to true sub _slice_fetch will use a simplified sql to retrieve 1bp slices.
  Returntype : boolean
  Exceptions : none
  Caller     : EnsemblGenomes variation DB build
  Status     : Stable
  
=cut

sub start_equals_end {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'start_equals_end'} = $value;
  }
  return $self->{'start_equals_end'};
}


=head2 clear_cache

  Args      : None
  Example   : my $sa =
                $registry->get_adaptor( 'Mus musculus', 'Core',
                                        'Slice' );
              my $ga =
                $registry->get_adaptor( 'Mus musculus', 'Core',
                                        'Gene' );

              my $slice =
                $sa->fetch_by_region( 'Chromosome', '1', 1e8,
                                      1.05e8 );

              my $genes = $ga->fetch_all_by_Slice($slice);

              $ga->clear_cache();

  Description   : Empties the feature cache associated with this
                  feature adaptor.
  Return type   : None
  Exceptions    : None
  Caller        : General
  Status        : At risk (under development)

=cut

sub clear_cache {
  my ($self) = @_;
  $self->_clear_slice_feature_cache();
  if(!$self->_no_id_cache()) {
    $self->_id_cache()->clear_cache();
  }
  return;
}

sub _clear_slice_feature_cache {
  my ($self) = @_;
  %{$self->{_slice_feature_cache}} = ();
  return;
}

=head2 _slice_feature_cache
 
  Description	: Returns the feature cache if we are allowed to cache and
                will build it if we need to. We will never return a reference
                to the hash to avoid unintentional auto-vivfying caching
  Returntype 	: Bio::EnsEMBL::Utils::Cache
  Exceptions 	: None
  Caller     	: Internal

=cut

sub _slice_feature_cache {
  my ($self) = @_;
  return if $self->db()->no_cache();
  if(! exists $self->{_slice_feature_cache}) {
    tie my %cache, 'Bio::EnsEMBL::Utils::Cache', $SLICE_FEATURE_CACHE_SIZE;
    $self->{_slice_feature_cache} = \%cache;
  }
  return $self->{_slice_feature_cache};
}

=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Arg [2]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Example    : $fts = $a->fetch_all_by_Slice($slice, 'Swall');
  Description: Returns a listref of features created from the database 
               which are on the Slice defined by $slice. If $logic_name is 
               defined only features with an analysis of type $logic_name 
               will be returned. 
               NOTE: only features that are entirely on the slice's seq_region
               will be returned (i.e. if they hang off the start/end of a
               seq_region they will be discarded). Features can extend over the
               slice boundaries though (in cases where you have a slice that
               doesn't span the whole seq_region).
  Returntype : listref of Bio::EnsEMBL::SeqFeatures in Slice coordinates
  Exceptions : none
  Caller     : Bio::EnsEMBL::Slice
  Status     : Stable

=cut

sub fetch_all_by_Slice {
  my ($self, $slice, $logic_name) = @_;
  #fetch by constraint with empty constraint
  return $self->fetch_all_by_Slice_constraint($slice, '', $logic_name);
}



=head2 fetch_Iterator_by_Slice_method

  Arg [1]    : CODE ref of Slice fetch method
  Arg [2]    : ARRAY ref of parameters for Slice fetch method
  Arg [3]    : Optional int: Slice index in parameters array
  Arg [4]    : Optional int: Slice chunk size. Default=500000
  Example    : my $slice_iter = $feature_adaptor->fetch_Iterator_by_Slice_method
                               	      ($feature_adaptor->can('fetch_all_by_Slice_Arrays'),
	                                   \@fetch_method_params,
	                                   0,#Slice idx
	                                  );

               while(my $feature = $slice_iter->next && defined $feature){
                 #Do something here
               }

  Description: Creates an Iterator which chunks the query Slice to facilitate
               large Slice queries which would have previously run out of memory
  Returntype : Bio::EnsEMBL::Utils::Iterator
  Exceptions : Throws if mandatory params not valid
  Caller     : general
  Status     : at risk

=cut

#Does not support Collections. See Funcgen ResultFeatureAdaptor::fetch_collection_Iterator_by_Slice_method

sub fetch_Iterator_by_Slice_method{
  my ($self, $slice_method_ref, $params_ref, $slice_idx, $chunk_size) = @_;

  if(! ( defined $slice_method_ref &&
		 ref($slice_method_ref) eq 'CODE')
	){
	throw('Must pass a valid Slice fetch method CODE ref');
  }

  if (! ($params_ref && 
		 ref($params_ref) eq 'ARRAY')) {
	#Don't need to check size here so long as we have valid Slice
	throw('You must pass a method params ARRAYREF');
  }
  
  $slice_idx    = 0 if(! defined $slice_idx);
  my $slice     = $params_ref->[$slice_idx];
  $chunk_size ||= 1000000;
		
  my @feat_cache;
  my $finished     = 0;
  my $start        = 1;	#local coord for sub slice
  my $end          = $slice->length;
  my $num_overlaps = 0;
  
  my $coderef = 
	sub {
	  
	  while (scalar(@feat_cache) == 0 &&
			 ! $finished) {
		
		my $new_end = ($start + $chunk_size - 1);
		
		if ($new_end >= $end) {
		  # this is our last chunk
		  $new_end = $end;
		  $finished = 1;  
		}
		
		#Chunk by sub slicing
		my $sub_slice             = $slice->sub_Slice($start, $new_end);
		$params_ref->[$slice_idx] = $sub_slice;
		@feat_cache = @{ $slice_method_ref->($self, @$params_ref)};
		
		#Remove & count overlapping features
		splice(@feat_cache, 0, $num_overlaps) if($num_overlaps);
		my $i;
		
		if (scalar(@feat_cache) > 0) {
		  
		  my $feat_end  = $feat_cache[$#feat_cache]->seq_region_end;
		  my $slice_end = $sub_slice->end;
		  $num_overlaps = 0;
		  for ($i = $#feat_cache; $i >=0; $i--) {
                        $feat_end  = $feat_cache[$i]->seq_region_end;
			
			if ($feat_end > $slice_end) {
			  $num_overlaps ++;
			} else {
			  last;
			}
			
		  }
		}
		
		# update the start coordinate
		$start = $new_end + 1;
	  }
	  
	  #this maybe returning from an undef cache
	  #Need to sub this out even more?
	  return shift @feat_cache;
	};

  return Bio::EnsEMBL::Utils::Iterator->new($coderef);
}


=head2 fetch_Iterator_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Optional string: logic name of analysis
  Arg [3]    : Optional int: Chunk size to iterate over. Default is 500000
  Example    : my $slice_iter = $feature_adaptor->fetch_Iterator_by_Slice($slice);

               while(my $feature = $slice_iter->next && defined $feature){
                 #Do something here
               }

  Description: Creates an Iterator which chunks the query Slice to facilitate
               large Slice queries which would have previously run out of memory
  Returntype : Bio::EnsEMBL::Utils::Iterator
  Exceptions : None
  Caller     : general
  Status     : at risk

=cut

sub fetch_Iterator_by_Slice{
  my ($self, $slice, $logic_name, $chunk_size) = @_;

  my $method_ref = $self->can('fetch_all_by_Slice');

  return $self->fetch_Iterator_by_Slice_method($method_ref, [$slice, $logic_name], 0, $chunk_size);
}


=head2 fetch_all_by_Slice_and_score

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Arg [2]    : (optional) float $score
               lower bound of the the score of the features retrieved
  Arg [3]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Example    : $fts = $a->fetch_all_by_Slice_and_score($slice,90,'Swall');
  Description: Returns a list of features created from the database which are 
               are on the Slice defined by $slice and which have a score 
               greater than $score. If $logic_name is defined, 
               only features with an analysis of type $logic_name will be 
               returned. 
  Returntype : listref of Bio::EnsEMBL::SeqFeatures in Slice coordinates
  Exceptions : none
  Caller     : Bio::EnsEMBL::Slice
  Status     : Stable

=cut

sub fetch_all_by_Slice_and_score {
  my ( $self, $slice, $score, $logic_name ) = @_;

  my $constraint;
  if ( defined($score) ) {
    # Get the synonym of the primary_table
    my @tabs = $self->_tables();
    my $syn  = $tabs[0]->[1];

    $constraint = sprintf( "%s.score > %s",
                $syn,
                $self->dbc()->db_handle()->quote( $score, SQL_FLOAT ) );
  }

  return
    $self->fetch_all_by_Slice_constraint( $slice, $constraint,
                                          $logic_name );
}


=head2 fetch_all_by_Slice_constraint

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Arg [2]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Arg [3]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Example    : $fs = $a->fetch_all_by_Slice_constraint($slc, 'perc_ident > 5');
  Description: Returns a listref of features created from the database which 
               are on the Slice defined by $slice and fulfill the SQL 
               constraint defined by $constraint. If logic name is defined, 
               only features with an analysis of type $logic_name will be 
               returned. 
  Returntype : listref of Bio::EnsEMBL::SeqFeatures in Slice coordinates
  Exceptions : thrown if $slice is not defined
  Caller     : Bio::EnsEMBL::Slice
  Status     : Stable

=cut

sub fetch_all_by_Slice_constraint {
  my ( $self, $slice, $constraint, $logic_name ) = @_;

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
  $constraint =
    $self->_logic_name_to_constraint( $constraint, $logic_name );

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
      #IMPORTANT: NEVER EVER RETURN A COPY OF THE DATA STRUCTURE.
      #           This will hold arrays of values. Since we've been doing
      #           this for so long people are expecting multiple calls
      #           to fetch_by_SliceXXXXX() methods to return the same
      #           array reference.
      return $cache->{$key};
    }
  } ## end if ( !( defined( $self...)))


  my $proj_ref = $self->_get_and_filter_Slice_projections($slice);
  my $bounds = $self->_generate_feature_bounds($slice); 

  # fetch features for the primary slice AND all symlinked slices
  foreach my $seg (@{$proj_ref}) {
    # re-bind the params
    $self->_bind_param_generic_fetch($bind_params); 
    my $offset    = $seg->from_start();
    my $seg_slice = $seg->to_Slice();
    my $features =
      $self->_slice_fetch( $seg_slice, $constraint );

    # If this was a symlinked slice offset the feature coordinates as
    # needed.
    if ( $seg_slice->name() ne $slice->name() ) {

    FEATURE:
      foreach my $f ( @{$features} ) {
        if ( $offset != 1 ) {
          $f->{'start'} += $offset - 1;
          $f->{'end'}   += $offset - 1;
        }

        # discard boundary crossing features from symlinked regions
        foreach my $bound (@{$bounds}) {
          if ( $f->{'start'} < $bound && $f->{'end'} >= $bound ) {
            next FEATURE;
          }
        }

        $f->{'slice'} = $slice;
        push( @result, $f );
      }
    } else {
      push( @result, @{$features} );
    }
  } ## end foreach my $seg (@proj)

  # Will only use feature_cache when set attribute no_cache in DBAdaptor
  if ( defined($key) ) {
    $cache->{$key} = \@result;
  }

  return \@result;
} ## end sub fetch_all_by_Slice_constraint


=head2 fetch_all_by_logic_name

  Arg [1]    : string $logic_name
               the logic name of the type of features to obtain
  Example    : $fs = $a->fetch_all_by_logic_name('foobar');
  Description: Returns a listref of features created from the database.
               only features with an analysis of type $logic_name will
               be returned.  If the logic name is invalid (not in the
               analysis table), a reference to an empty list will be
               returned.
  Returntype : listref of Bio::EnsEMBL::SeqFeatures
  Exceptions : thrown if no $logic_name
  Caller     : General
  Status     : Stable

=cut

sub fetch_all_by_logic_name {
  my ( $self, $logic_name ) = @_;

  if ( !defined($logic_name) ) {
    throw("Need a logic_name");
  }

  my $constraint = $self->_logic_name_to_constraint( '', $logic_name );

  if ( !defined($constraint) ) {
    warning("Invalid logic name: $logic_name");
    return [];
  }

  return $self->generic_fetch($constraint);
}

=head2 fetch_all_by_stable_id_list

  Arg [1]    : string $logic_name
               the logic name of the type of features to obtain
  Arg [2]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Example    : $fs = $a->fetch_all_by_stable_id_list(["ENSG00001","ENSG00002", ...]);
  Description: Returns a listref of features identified by their stable IDs.
               This method only fetches features of the same type as the calling
               adaptor. 
               Results are constrained to a slice if the slice is provided.
  Returntype : listref of Bio::EnsEMBL::Feature
  Exceptions : thrown if no stable ID list is provided.
  Caller     : General
  Status     : Stable

=cut

# Adapted from BaseAdaptor->uncached_fetch_all_by_dbID_list
sub fetch_all_by_stable_id_list {
  my ( $self, $id_list_ref, $slice ) = @_;

  return $self->_uncached_fetch_all_by_id_list($id_list_ref,$slice,"stable_id");
}

# Method that creates an object.  Called by the _objs_from_sth() method
# in the sub-classes (the various feature adaptors).  Overridden by the
# feature collection classes.

sub _create_feature {
  my ( $self, $feature_type, $args ) = @_;
  return $feature_type->new( %{$args} );
}

# This is the same as the above, but calls the new_fast() constructor of
# the feature type.

sub _create_feature_fast {
  my ( $self, $feature_type, $args ) = @_;
  return $feature_type->new_fast($args);
}

=head2 count_by_Slice_constraint

    Arg [1]     : Bio::EnsEMBL::Slice
    Arg [2]     : String Custom SQL constraint
    Arg [3]     : String Logic name to search by
    Description : Finds all features with at least partial overlap to the given
                  slice and sums them up.
                  Explanation of workings with projections:
                  
                  |-------------------------Constraint Slice---------------------------------|
                  |             |                                           |                |
                  |--Segment 1--|                                           |                |
                  |             |--Segment 2, on desired Coordinate System--|                |
                  |             |                                           |---Segment 3----|
              #Feature 1#    #Feature 2#                               #Feature 3#           |
                  |         #####################Feature 4####################               |
                  | #Feature 5# |                                           |                |
                  
                  Feature 1 is overlapping the original constraint. Counted in Segment 1
                  Feature 2,3 and 4  are counted when inspecting Segment 2
                  Feature 5 is counted in Segment 1
                  
    Returntype  : Integer
=cut

sub count_by_Slice_constraint {
    my ($self, $slice, $constraint, $logic_name) = @_;
    
    assert_ref($slice, 'Bio::EnsEMBL::Slice', 'slice');
    
    my $count = 0;
    
    #Remember the bind params
    my $bind_params = $self->bind_param_generic_fetch();
    
    #Table synonym
    my @tables = $self->_tables;
    my ($table_name, $table_synonym) = @{ $tables[0] };
    
    # Basic constraints limit the feature hits to the starting Slice constraint
    $constraint ||= '';
    $constraint = $self->_logic_name_to_constraint( $constraint, $logic_name );
  
    return $count if ! defined $constraint;
  
    #Query logic
    my $sa = $slice->adaptor();
    my $projections = $self->_get_and_filter_Slice_projections($slice);

    my $segment_count = @{$projections};
    #Manual loop to support look-ahead/behind
    for (my $i = 0; $i < $segment_count; $i++) {
        my $seg = $projections->[$i];
        my $seg_slice = $seg->to_Slice();
        $self->_bind_param_generic_fetch($bind_params);
    
        # We cannot filter boundary crossing features in code, so we constrain the
        # SQL. We detect when we are not on the original query Slice and manually
        # filter by the segment Slice's start and end. If we are on "earlier" (5')
        # projection segments, we only limit by the *end* and when we are on the later
        # projection segments, we filter by the *start*.
        
        # This is the same as fetch_all_by_Slice_constraint()'s in-memory filtering
        # except we need to alter projected features in that code with an offset
        my $local_constraint = $constraint;
        if ( $seg_slice->name() ne $slice->name() ) {
            my ($limit_start, $limit_end) = (1,1); #fully constrained by default, flags are reset every iteration
            if($i == 0) {
                $limit_start = 0;
            } elsif ($i == ($segment_count - 1)) {
                $limit_end = 0; #don't check end as we are on the final projection
            }
            $local_constraint .= ' AND ' if $local_constraint;
            my @conditionals;
      
            if($limit_end) {
            # Restrict features to the end of this projection segment when after the named segment 
                push(@conditionals, sprintf('%1$s.seq_region_end <= %2$d', $table_synonym, $seg_slice->end));
            }
            if($limit_start) {
            #Do not cross the start boundary on this projection segment
                 push(@conditionals, sprintf('%1$s.seq_region_start >= %2$d', $table_synonym, $seg_slice->start));
            }
      
            $local_constraint .= join(q{ AND }, @conditionals);
        }
    
  	    my $count_array = $self->_get_by_Slice($seg_slice, $local_constraint, 'count');
  	    $count += $_ for @$count_array;
    }
	
	return $count;
}

=head2 _get_and_filter_Slice_projections

    Arg [1]     : Bio::EnsEMBL::Slice
    Description : Delegates onto SliceAdaptor::fetch_normalized_slice_projection() 
                  with filtering on
    Returntype  : ArrayRef Bio::EnsEMBL::ProjectionSegment; Returns an array
                  of projected segments
=cut

sub _get_and_filter_Slice_projections {
  my ($self, $slice) = @_;
  my $sa = $slice->adaptor();
  my $filter_projections = 1;
  return $sa->fetch_normalized_slice_projection($slice, $filter_projections);
}

=head2 _generate_feature_bounds

    Arg [1]     : Bio::EnsEMBL::Slice
    Description : Performs a projection of Slice and records the bounds
                  of that projection. This can be used later on to restrict
                  Features which overlap into unwanted areas such as
                  regions which exist on another HAP/PAR region.
                  
                  Bounds are defined as projection_start - slice_start + 1.
    Example     : my $bounds = $self->_generate_feature_bounds($slice);
    Returntype  : ArrayRef Integer; Returns the location of the bounds.
=cut

sub _generate_feature_bounds {
  my ($self, $slice) = @_;
  my $sa = $slice->adaptor();
  # construct list of Hap/PAR boundaries for entire seq region
  my @bounds;

  my $sr_id = $slice->get_seq_region_id();
  my $ent_slice = $sa->fetch_by_seq_region_id($sr_id);
  if ( $slice->strand() == -1 ) {
    $ent_slice = $ent_slice->invert();
  }

  my @ent_proj = @{ $sa->fetch_normalized_slice_projection($ent_slice) };
  shift(@ent_proj);    # skip first; 1st does not have bounds normally; may change if we ever have a patch a pos 1 

  @bounds = map { $_->from_start() - $slice->start() + 1 } @ent_proj;
  
  return \@bounds;
}

=head2 _get_by_Slice
    Arg [0]    : Bio::EnsEMBL::Slice to find all the features within
    Arg [1]    : SQL constraint string
    Arg [2]    : Type of query to run. Default behaviour is to select, but 
                 'count' is also valid
    Description: Abstracted logic from _slice_fetch
    Returntype : Listref of Bio::EnsEMBL::Feature, or integers for counting mode
=cut

sub _get_by_Slice {
    my ($self, $slice, $orig_constraint, $query_type) = @_;
    
    # features can be scattered across multiple coordinate systems
    my @tables = $self->_tables;
    my ($table_name, $table_synonym) = @{ $tables[0] };
    my $mapper;
    my @feature_coord_systems;
    
    my $meta_values = $self->db->get_MetaContainer->list_value_by_key( $table_name."build.level");
    if ( @$meta_values and $slice->is_toplevel() ) {
        push @feature_coord_systems, $slice->coord_system();
    } else {
        @feature_coord_systems = @{ $self->db->get_MetaCoordContainer->fetch_all_CoordSystems_by_feature_type($table_name)};
    }
	
    my $assembly_mapper_adaptor = $self->db->get_AssemblyMapperAdaptor();
    my @pan_coord_features;
        
COORD_SYSTEM: foreach my $coord_system (@feature_coord_systems) {
        my @query_accumulator;
        # Build up a combination of query constraints that will quickly establish the result set
        my $constraint = $orig_constraint;
        if ( $coord_system->equals( $slice->coord_system ) ) {
            my $max_len = $self->_max_feature_length
                || $self->db->get_MetaCoordContainer
                    ->fetch_max_length_by_CoordSystem_feature_type( $coord_system,$table_name );
                       
            my $seq_region_id;
            if ( $slice->adaptor ) {
                $seq_region_id = $slice->adaptor->get_seq_region_id($slice);
            } else {
                $seq_region_id = $self->db->get_SliceAdaptor->get_seq_region_id($slice);
            }
            
            my @seq_region_ids = ($seq_region_id);
            while (1) {
                my $ext_seq_region_id = $self->get_seq_region_id_external($seq_region_id);
        
                if ( $ext_seq_region_id == $seq_region_id ) { last }
        
                push( @seq_region_ids, $ext_seq_region_id );
                $seq_region_id = $ext_seq_region_id;
            }
            
            $constraint .= " AND " if ($constraint);

            $constraint .= ${table_synonym}.".seq_region_id IN (". join( ',', @seq_region_ids ) . ") AND ";
            
            #faster query for 1bp slices where SNP data is not compressed
            if ( $self->start_equals_end && $slice->start == $slice->end ) {
                $constraint .= " AND ".$table_synonym.".seq_region_start = ".$slice->end .
                  " AND ".$table_synonym.".seq_region_end = ".$slice->start;
            
            } else {
                if ( !$slice->is_circular() ) {
                    # Deal with the default case of a non-circular chromosome.
                    $constraint .= $table_synonym.".seq_region_start <= ".$slice->end." AND "
                                   .$table_synonym.".seq_region_end >= ".$slice->start;
            
                    if ( $max_len ) {
                        my $min_start = $slice->start - $max_len;
                        $constraint .= " AND ".$table_synonym.".seq_region_start >= ".$min_start;
                    }
            
                } else {
                    # Deal with the case of a circular chromosome.
                    if ( $slice->start > $slice->end ) {
                        $constraint .= " ( ".$table_synonym.".seq_region_start >= ".$slice->start
                            . " OR ".$table_synonym.".seq_region_start <= ".$slice->end
                            . " OR ".$table_synonym.".seq_region_end >= ".$slice->start
                            . " OR ".$table_synonym.".seq_region_end <= ".$slice->end
                            . " OR ".$table_synonym.".seq_region_start > ".$table_synonym.".seq_region_end)";
                    } else {
                        $constraint .= " ((".$table_synonym.".seq_region_start <= ".$slice->end
                            . " AND ".$table_synonym.".seq_region_end >= ".$slice->start.") "
                            . "OR (".$table_synonym.".seq_region_start > ".$table_synonym.".seq_region_end"
                            . " AND (".$table_synonym.".seq_region_start <= ".$slice->end
                            . " OR ".$table_synonym.".seq_region_end >= ".$slice->start.")))";
                  }
              }
           }
           push @query_accumulator, [$constraint,undef,$slice]; # $mapper intentionally absent here.
           
        } else { 
        	#coordinate systems do not match
            $mapper = $assembly_mapper_adaptor->fetch_by_CoordSystems( $slice->coord_system(), $coord_system );
            next unless defined $mapper;

            # Get list of coordinates and corresponding internal ids for
            # regions the slice spans
            my @coords = $mapper->map( $slice->seq_region_name, $slice->start, $slice->end,
                                    $slice->strand, $slice->coord_system );

            @coords = grep { !$_->isa('Bio::EnsEMBL::Mapper::Gap') } @coords;

            next COORD_SYSTEM if ( !@coords );

            my @ids = map { $_->id() } @coords;
            #coords are now id rather than name
            
            if ( @coords > $MAX_SPLIT_QUERY_SEQ_REGIONS && ! $slice->isa('Bio::EnsEMBL::LRGSlice') 
                    && $slice->coord_system->name() ne 'lrg') {
                $constraint = $orig_constraint;
                my $id_str = join( ',', @ids );
                $constraint .= " AND " if ($constraint);
                $constraint .= $table_synonym.".seq_region_id IN ($id_str)";
                
                push @query_accumulator, [$constraint,$mapper,$slice];
            } else {
                my $max_len = (
                    $self->_max_feature_length()
                    || $self->db->get_MetaCoordContainer
                       ->fetch_max_length_by_CoordSystem_feature_type($coord_system, $table_name) 
                );

                my $length = @coords;
                for ( my $i = 0; $i < $length; $i++ ) {
                    $constraint = $orig_constraint;
                    $constraint .= " AND " if ($constraint);
                    $constraint .= $table_synonym.".seq_region_id = "
                        . $ids[$i] . " AND "
                        . $table_synonym.".seq_region_start <= "
                        . $coords[$i]->end() . " AND "
                        . $table_synonym.".seq_region_end >= "
                        . $coords[$i]->start();

                    if ($max_len) {
                        my $min_start = $coords[$i]->start() - $max_len;
                        $constraint .= " AND ".$table_synonym.".seq_region_start >= ".$min_start;
                    }
                    
                    push @query_accumulator, [$constraint,$mapper,$slice];
                } # end multi-query cycle
        } # end else
            
     } # end else (coord sytems not matching)
     
     #Record the bind params if we have to do multiple queries
     my $bind_params = $self->bind_param_generic_fetch();
     
     foreach my $query (@query_accumulator) {
         my ($local_constraint,$local_mapper,$local_slice) = @$query;
         $self->_bind_param_generic_fetch($bind_params);
         if ($query_type and $query_type eq 'count') {
           push @pan_coord_features, $self->generic_count($local_constraint);
         } 
         else {
             my $features = $self->generic_fetch( $local_constraint, $local_mapper, $local_slice );
             $features = $self->_remap( $features, $local_mapper, $local_slice );
             push @pan_coord_features, @$features;
         }
     }
     $mapper = undef;
    } # End foreach
    $self->{_bind_param_generic_fetch} = undef;
    return \@pan_coord_features;
}
#
# helper function used by fetch_all_by_Slice_constraint method
#
sub _slice_fetch {
  my ( $self, $slice, $orig_constraint ) = @_;

  my $features = $self->_get_by_Slice($slice,$orig_constraint);

  return $features;
} ## end sub _slice_fetch


#for a given seq_region_id, gets the one used in an external database, if present, otherwise, returns the internal one
sub get_seq_region_id_external {
  my ( $self, $sr_id ) = @_;
  my $cs_a = $self->db()->get_CoordSystemAdaptor();
  return ( exists( $cs_a->{'_internal_seq_region_mapping'}->{$sr_id} )
           ? $cs_a->{'_internal_seq_region_mapping'}->{$sr_id}
           : $sr_id );
}

#for a given seq_region_id and coord_system, gets the one used in the internal (core) database
sub get_seq_region_id_internal{
  my ( $self, $sr_id ) = @_;
  my $cs_a = $self->db()->get_CoordSystemAdaptor();
  return (  exists $cs_a->{'_external_seq_region_mapping'}->{$sr_id} 
            ? $cs_a->{'_external_seq_region_mapping'}->{$sr_id} 
            : $sr_id);
}

#
# Helper function containing some common feature storing functionality
#
# Given a Feature this will return a copy (or the same feature if no changes 
# to the feature are needed) of the feature which is relative to the start
# of the seq_region it is on. The seq_region_id of the seq_region it is on
# is also returned.
#
# This method will also ensure that the database knows which coordinate
# systems that this feature is stored in.
#

sub _pre_store {
  my $self    = shift;
  my $feature = shift;

  if(!ref($feature) || !$feature->isa('Bio::EnsEMBL::Feature')) {
    throw('Expected Feature argument.');
  }
  my $slice = $feature->slice();

  $self->_check_start_end_strand($feature->start(),$feature->end(),
                                 $feature->strand(), $slice);


  my $db = $self->db();

  my $slice_adaptor = $db->get_SliceAdaptor();

  if(!ref($slice) || !($slice->isa('Bio::EnsEMBL::Slice') or $slice->isa('Bio::EnsEMBL::LRGSlice'))  ) {
    throw('Feature must be attached to Slice to be stored.');
  }

  # make sure feature coords are relative to start of entire seq_region

  if($slice->start != 1 || $slice->strand != 1) {
    #move feature onto a slice of the entire seq_region
    $slice = $slice_adaptor->fetch_by_region($slice->coord_system->name(),
                                             $slice->seq_region_name(),
                                             undef, #start
                                             undef, #end
                                             undef, #strand
                                             $slice->coord_system->version());

    $feature = $feature->transfer($slice);

    if(!$feature) {
      throw('Could not transfer Feature to slice of ' .
            'entire seq_region prior to storing');
    }
  }

  # Ensure this type of feature is known to be stored in this coord system.
  my $cs = $slice->coord_system;

  my ($tab) = $self->_tables();
  my $tabname = $tab->[0];

  my $mcc = $db->get_MetaCoordContainer();

  $mcc->add_feature_type($cs, $tabname, $feature->length);

  my $seq_region_id = $slice_adaptor->get_seq_region_id($slice);

  if(!$seq_region_id) {
    throw('Feature is associated with seq_region which is not in this DB.');
  }

  return ($feature, $seq_region_id);
}


# The same function as _pre_store
# This one is used to store user uploaded features in XXX_userdata db

sub _pre_store_userdata {
  my $self    = shift;
  my $feature = shift;

  if(!ref($feature) || !$feature->isa('Bio::EnsEMBL::Feature')) {
    throw('Expected Feature argument.');
  }

  my $slice = $feature->slice();
  my $slice_adaptor = $slice->adaptor;
  
  $self->_check_start_end_strand($feature->start(),$feature->end(),
                                 $feature->strand(), $slice);


  if(!ref($slice) || !($slice->isa('Bio::EnsEMBL::Slice') or $slice->isa('Bio::EnsEMBL::LRGSlice')) ) {
    throw('Feature must be attached to Slice to be stored.');
  }

  # make sure feature coords are relative to start of entire seq_region

  if($slice->start != 1 || $slice->strand != 1) {
    #move feature onto a slice of the entire seq_region
    $slice = $slice_adaptor->fetch_by_region($slice->coord_system->name(),
                                             $slice->seq_region_name(),
                                             undef, #start
                                             undef, #end
                                             undef, #strand
                                             $slice->coord_system->version());

    $feature = $feature->transfer($slice);

    if(!$feature) {
      throw('Could not transfer Feature to slice of ' .
            'entire seq_region prior to storing');
    }
  }

  # Ensure this type of feature is known to be stored in this coord system.
  my $cs = $slice->coord_system;

  my ($tab) = $self->_tables();
  my $tabname = $tab->[0];

  my $db = $self->db;
  my $mcc = $db->get_MetaCoordContainer();

  $mcc->add_feature_type($cs, $tabname, $feature->length);

  my $seq_region_id = $slice_adaptor->get_seq_region_id($slice);

  if(!$seq_region_id) {
    throw('Feature is associated with seq_region which is not in this DB.');
  }

  return ($feature, $seq_region_id);
}


#
# helper function used to validate start/end/strand and 
# hstart/hend/hstrand etc.
#
sub _check_start_end_strand {
  my ($self, $start, $end, $strand, $slice) = @_;

  #
  # Make sure that the start, end, strand are valid
  #
  if(int($start) != $start) {
    throw("Invalid Feature start [$start].  Must be integer.");
  }
  if(int($end) != $end) {
    throw("Invalid Feature end [$end]. Must be integer.");
  }
  if(int($strand) != $strand || $strand < -1 || $strand > 1) {
    throw("Invalid Feature strand [$strand]. Must be -1, 0 or 1.");
  }
  if($end < $start && ! (defined $slice && $slice->is_circular())) {
    throw("Invalid Feature start/end [$start/$end]. Start must be less " .
          "than or equal to end.");
  }

  return 1;
}


#
# Given a list of features checks if they are in the correct coord system
# by looking at the first features slice.  If they are not then they are
# converted and placed on the slice.
#
sub _remap {
  my ( $self, $features, $mapper, $slice ) = @_;

  #check if any remapping is actually needed
  if(@$features && (!$features->[0]->isa('Bio::EnsEMBL::Feature') ||
                    $features->[0]->slice == $slice)) {
    return $features;
  }

  #remapping has not been done, we have to do our own conversion from
  #to slice coords

  my @out;

  my $slice_start = $slice->start();
  my $slice_end   = $slice->end();
  my $slice_strand = $slice->strand();
  my $slice_cs    = $slice->coord_system();

  my ($seq_region_id, $start, $end, $strand);

  my $slice_seq_region_id = $slice->get_seq_region_id();

  foreach my $f (@$features) {
    #since feats were obtained in contig coords, attached seq is a contig
    my $fslice = $f->slice();
    if(!$fslice) {
      throw("Feature does not have attached slice.\n");
    }
    my $fseq_region_id = $fslice->get_seq_region_id();
    my $fcs = $fslice->coord_system();

    if(!$slice_cs->equals($fcs)) {
      #slice of feature in different coord system, mapping required
      # $seq_region comes out as an integer not a string
      ($seq_region_id, $start, $end, $strand) =
        $mapper->fastmap($fslice->seq_region_name(),$f->start(),$f->end(),$f->strand(),$fcs);

      # undefined start means gap
      next if(!defined $start);
    } else {
      $start          = $f->start();
      $end            = $f->end();
      $strand         = $f->strand();
      $seq_region_id  = $fseq_region_id;
    }

    # maps to region outside desired area
    next if ($start > $slice_end) || ($end < $slice_start) || 
      ($slice_seq_region_id != $seq_region_id);

    #shift the feature start, end and strand in one call
    if($slice_strand == -1) {
      $f->move( $slice_end - $end + 1, $slice_end - $start + 1, $strand * -1 );
    } else {
      $f->move( $start - $slice_start + 1, $end - $slice_start + 1, $strand );
    }

    $f->slice($slice);

    push @out,$f;
  }

  return \@out;
}

#
# get seq region boundary (start|end) for a feature
# the method attempts to retrieve the boundary directly from the db
# return undef if cannot determine the associated feature table
#
sub _seq_region_boundary_from_db {
  my ($self, $feature, $boundary) = @_;

  if(!ref($feature) || !$feature->isa('Bio::EnsEMBL::Feature')) {
    throw('Expected Feature argument.');
  }

  throw "Undefined boundary"
    unless defined $boundary;

  $boundary eq 'start' or $boundary eq 'end'
    or throw "Wrong boundary: select start|end";

  $boundary = 'seq_region_' . $boundary;

  my $sql_helper = $self->dbc->sql_helper;
  throw "Unable to get SqlHelper instance" unless defined $sql_helper;

  my $feature_table = ($self->_tables)[0]->[0];
  warning (sprintf "Unable to get %s for %s instance\nCould not find associated feature table, returning undef", $boundary, ref $feature)
    and return undef unless defined $feature_table;

  my $db_id = $feature->dbID;
  my $attrib_id = $feature_table . '_id';
  my $query = "SELECT ${boundary} from ${feature_table} WHERE ${attrib_id} = ${db_id}"; 

  return $sql_helper->execute_single_result(-SQL => $query);
}

=head2 store

  Arg [1]    : list of Bio::EnsEMBL::SeqFeature
  Example    : $adaptor->store(@feats);
  Description: ABSTRACT  Subclasses are responsible for implementing this 
               method.  It should take a list of features and store them in 
               the database.
  Returntype : none
  Exceptions : thrown method is not implemented by subclass
  Caller     : general
  Status     : At Risk
             : throws if called.

=cut

sub store{
  my $self = @_;

  throw("Abstract method store not defined by implementing subclass\n");
}


=head2 remove

  Arg [1]    : A feature $feature 
  Example    : $feature_adaptor->remove($feature);
  Description: This removes a feature from the database.  The table the
               feature is removed from is defined by the abstract method
               _tablename, and the primary key of the table is assumed
               to be _tablename() . '_id'.  The feature argument must 
               be an object implementing the dbID method, and for the
               feature to be removed from the database a dbID value must
               be returned.
  Returntype : none
  Exceptions : thrown if $feature arg does not implement dbID(), or if
               $feature->dbID is not a true value
  Caller     : general
  Status     : Stable

=cut


sub remove {
  my ($self, $feature) = @_;

  if(!$feature || !ref($feature) || !$feature->isa('Bio::EnsEMBL::Feature')) {
    throw('Feature argument is required');
  }

  if(!$feature->is_stored($self->db)) {
    throw("This feature is not stored in this database");
  }

  my @tabs = $self->_tables;
  my ($table) = @{$tabs[0]};

  my $sth = $self->prepare("DELETE FROM $table WHERE ${table}_id = ?");
  $sth->bind_param(1,$feature->dbID,SQL_INTEGER);
  $sth->execute();

  #unset the feature dbID ad adaptor
  $feature->dbID(undef);
  $feature->adaptor(undef);

  return;
}


=head2 remove_by_Slice

  Arg [1]    : Bio::Ensembl::Slice $slice
  Example    : $feature_adaptor->remove_by_Slice($slice);
  Description: This removes features from the database which lie on a region
               represented by the passed in slice.  Only features which are
               fully contained by the slice are deleted; features which overlap
               the edge of the slice are not removed.
               The table the features are removed from is defined by
               the abstract method_tablename.
  Returntype : none
  Exceptions : thrown if no slice is supplied
  Caller     : general
  Status     : Stable

=cut

sub remove_by_Slice {
  my ($self, $slice) = @_;

  if(!$slice || !ref($slice) || !($slice->isa('Bio::EnsEMBL::Slice') or $slice->isa('Bio::EnsEMBL::LRGSlice')) ) {
    throw("Slice argument is required");
  }

  my @tabs = $self->_tables;
  my ($table_name) = @{$tabs[0]};

  my $seq_region_id = $self->db->get_SliceAdaptor->get_seq_region_id($slice);
  my $start = $slice->start();
  my $end   = $slice->end();

  #
  # Delete only features fully on the slice, not overlapping ones
  #
  my $sth = $self->prepare("DELETE FROM $table_name " .
                           "WHERE seq_region_id = ? " .
                           "AND   seq_region_start >= ? " .
                           "AND   seq_region_end <= ?");

  $sth->bind_param(1,$seq_region_id,SQL_INTEGER);
  $sth->bind_param(2,$start,SQL_INTEGER);
  $sth->bind_param(3,$end,SQL_INTEGER);
  $sth->execute();
  $sth->finish();
}


#
# Internal function. Allows the max feature length which is normally
# retrieved from the meta_coord table to be overridden.  This allows
# for some significant optimizations to be put in when it is known
# that requested features will not be over a certain size.
#
sub _max_feature_length {
  my $self = shift;
  return $self->{'_max_feature_length'} = shift if(@_);
  return $self->{'_max_feature_length'};
}


#
# Lists all seq_region_ids that a particular feature type is found on.
# Useful e.g. for finding out which seq_regions have genes.
# Returns a listref of seq_region_ids.
#
sub _list_seq_region_ids {
  my ($self, $table) = @_;
  
  my @out;
  
  my $sql = qq(
  SELECT DISTINCT
            sr.seq_region_id
  FROM      seq_region sr,
            $table a,
            coord_system cs
  WHERE     sr.seq_region_id = a.seq_region_id
    AND     sr.coord_system_id = cs.coord_system_id
    AND     cs.species_id = ?);

  my $sth = $self->prepare($sql);

  $sth->bind_param( 1, $self->species_id(), SQL_INTEGER );

  $sth->execute();

  while (my ($id) = $sth->fetchrow) {
    push(@out, $id);
  }

  $sth->finish;

  return \@out;
}


sub remove_by_analysis_id {
  my ($self, $analysis_id) = @_;

  $analysis_id or throw("Must call with analysis id");

  my @tabs = $self->_tables;
  my ($tablename) = @{$tabs[0]};

  my $sql = "DELETE FROM $tablename WHERE analysis_id = $analysis_id";
#  warn "SQL : $sql";
      
  my $sth = $self->prepare($sql);
  $sth->execute();
  $sth->finish();
}

sub remove_by_feature_id {
  my ($self, $features_list) = @_;

  my @feats = @$features_list or throw("Must call store with features");

  my @tabs = $self->_tables;
  my ($tablename) = @{$tabs[0]};

  my $sql = sprintf "DELETE FROM $tablename WHERE ${tablename}_id IN (%s)", join ', ', @feats;
#  warn "SQL : $sql";
      
  my $sth = $self->prepare($sql);
  $sth->execute();
  $sth->finish();
}

=head2 fetch_nearest_by_Feature

  Arg [1]    : Reference Feature to start the search from
  Description: Searches iteratively outward from the starting feature until a nearby Feature is found
               If you require more than one result or more control of which features are returned, see
               fetch_all_nearest_by_Feature and fetch_all_by_outward_search. fetch_nearest_by_Feature
               is a convenience method.
  ReturnType : Bio::EnsEMBL::Feature
=cut

sub fetch_nearest_by_Feature {
  my $self = shift;
  my $feature = shift;
  my @hits = @{ $self->fetch_all_by_outward_search(-FEATURE => $feature, -MAX_RANGE => 1000000, -RANGE => 1000, -LIMIT => 1) };
  return $hits[0] if (scalar @hits > 0);
  return;
}


=head2 fetch_all_by_outward_search

  Arguments the same as fetch_all_nearest_by_Feature
  Arg [0]    : -MAX_RANGE : Set an upper limit on the search range, defaults to 10000 bp 
  Arg [1]    : -FEATURE ,Bio::EnsEMBL::Feature : 'Source' Feature to anchor the search for nearest Features
  Arg [2]    : -SAME_STRAND, Boolean (optional)  : Respect the strand of the source Feature with ref, only 
               returning Features on the same strand.
  Arg [3]    : -OPPOSITE_STRAND, Boolean (optional) : Find features on the opposite strand of the same
  Arg [4]    : -DOWNSTREAM/-UPSTREAM, (optional) : Search ONLY downstream or upstream from the source Feature.
               Can be omitted for searches in both directions.
  Arg [5]    : -RANGE, Int     : The size of the space to search for Features. Defaults to 1000 as a sensible starting point
  Arg [6]    : -NOT_OVERLAPPING, Boolean (optional) : Do not return Features that overlap the source Feature
  Arg [7]    : -FIVE_PRIME, Boolean (optional) : Determine range to a Feature by the 5' end, respecting strand
  Arg [8]    : -THREE_PRIME, Boolean (optional): Determine range to a Feature by the 3' end, respecting strand
  Arg [9]    : -LIMIT, Int     : The maximum number of Features to return, defaulting to one. Equally near features are all returned

  Description: Searches for features within the suggested -RANGE, and if it finds none, expands the search area
               until it satisfies -LIMIT or hits -MAX_RANGE. Useful if you don't know how far away the features
               might be, or if dealing with areas of high feature density. In the case of Variation Features, it is
               conceivable that a 2000 bp window might contain very many features, resulting in a slow and bloated
               response, thus the ability to explore outward in smaller sections can be useful.
  Returntype : Listref of [$feature,$distance]
=cut

sub fetch_all_by_outward_search {
  my $self = shift;
   my ($ref_feature, $respect_strand, $opposite_strand, $downstream, $upstream, $search_range,
       $limit,$not_overlapping,$five_prime,$three_prime, $max_range) =
        rearrange([qw(FEATURE SAME_STRAND OPPOSITE_STRAND DOWNSTREAM UPSTREAM RANGE LIMIT NOT_OVERLAPPING FIVE_PRIME THREE_PRIME MAX_RANGE)], @_);
  my $factor = 1;
  $limit ||= 1;
  $search_range ||= 1000;
  $max_range ||= 10000;
  my @results;
  while (scalar @results < $limit && $search_range <= $max_range) {
    $search_range = $search_range * $factor;
    @results = @{ 
      $self->fetch_all_nearest_by_Feature(-RANGE => $search_range, 
                                          -FEATURE => $ref_feature, 
                                          -SAME_STRAND => $respect_strand, 
                                          -OPPOSITE_STRAND => $opposite_strand,
                                          -DOWNSTREAM => $downstream,
                                          -UPSTREAM => $upstream,
                                          -LIMIT => $limit,
                                          -NOT_OVERLAPPING => $not_overlapping,
                                          -FIVE_PRIME => $five_prime,
                                          -THREE_PRIME => $three_prime,
                                          )
    };
    $factor++;
  }
  return \@results;
}

=head2 fetch_all_nearest_by_Feature

  Arg [1]    : -FEATURE ,Bio::EnsEMBL::Feature : 'Source' Feature to anchor the search for nearest Features
  Arg [2]    : -SAME_STRAND, Boolean (optional): Respect the strand of the source Feature with ref, only 
                                                 returning Features on the same strand
  Arg [3]    : -OPPOSITE_STRAND, Boolean (optional) : Find features on the opposite strand of the same
  Arg [4]    : -DOWNSTREAM/-UPSTREAM, (optional) : Search ONLY downstream or upstream from the source Feature.
               Can be omitted for searches in both directions.
  Arg [5]    : -RANGE, Int     : The size of the space to search for Features. Defaults to 1000 as a sensible starting point
  Arg [6]    : -NOT_OVERLAPPING, Boolean (optional) : Do not return Features that overlap the source Feature
  Arg [7]    : -FIVE_PRIME, Boolean (optional) : Determine range to a Feature by its 5' end, respecting strand
  Arg [8]    : -THREE_PRIME, Boolean (optional): Determine range to a Feature by its 3' end, respecting strand
  Arg [9]    : -LIMIT, Int     : The maximum number of Features to return, defaulting to one. Equally near features are all returned
  Example    : #To fetch the gene(s) with the nearest 5' end:
               $genes = $gene_adaptor->fetch_all_nearest_by_Feature(-FEATURE => $feat, -FIVE_PRIME => 1);

  Description: Gets the nearest Features to a given 'source' Feature. The Feature returned and the format of the result
               are non-obvious, please read on.

               When looking beyond the boundaries of the source Feature, the distance is measured to the nearest end 
               of that Feature to the nearby Feature's nearest end.
               If Features overlap the source Feature, then they are given a distance of zero but ordered by
               their proximity to the centre of the Feature.
               
               Features are found and prioritised within 1000 base pairs unless a -RANGE is given to the method. Any overlap with
               the search region is included, and the results can be restricted to upstream, downstream, forward strand or reverse

               The -FIVE_PRIME and -THREE_PRIME options allow searching for specific ends of nearby features, but still needs
               a -DOWN/UPSTREAM value and/or -NOT_OVERLAPPING to fulfil its most common application.


  Returntype : Listref containing an Arrayref of Bio::EnsEMBL::Feature objects and the distance
               [ [$feature, $distance] ... ]
  Caller     : general

=cut
sub fetch_all_nearest_by_Feature{
    my $self = shift;
    my ($ref_feature, $respect_strand, $opposite_strand, $downstream, $upstream, $search_range,$limit,$not_overlapping,$five_prime,$three_prime) =
        rearrange([qw(FEATURE SAME_STRAND OPPOSITE_STRAND DOWNSTREAM UPSTREAM RANGE LIMIT NOT_OVERLAPPING FIVE_PRIME THREE_PRIME)], @_);
    if ( !defined($search_range)) {
      $search_range ||= 1000;
    }
    $limit ||= 1;

    unless (defined($ref_feature) && $ref_feature->isa('Bio::EnsEMBL::Feature')) {
      throw ('fetch_all_nearest_by_Feature method requires a valid Ensembl Feature object to operate');
    }
    throw ('Do not specify both -upstream and -downstream. This is the same as no arguments') if ($upstream && $downstream);
    throw ('Do not use both -SAME_STRAND and -OPPOSITE_STRAND in one query') if ($opposite_strand && $respect_strand);
    my ($region_start,$region_end);
    
    # Define search box. Features overlapping the boundaries are found and considered
    my $slice = $ref_feature->feature_Slice;

    my $five_prime_expansion = 0;
    my $three_prime_expansion = 0;
    if ($upstream) {
      $five_prime_expansion = $search_range;
      $three_prime_expansion = $slice->length * -1;
      # shrinks the slice to exclude the features that do not exist beyond the end of the
      # considered feature
    } 
    elsif ($downstream) {
      $three_prime_expansion = $search_range;
      $five_prime_expansion = $slice->length * -1;
    } 
    else {
      $five_prime_expansion = $search_range;
      $three_prime_expansion = $search_range;
    }

    $slice = $slice->expand( $five_prime_expansion, $three_prime_expansion);
    my $sa = $self->db->get_SliceAdaptor;
    my @candidates; # Features in the search region
    @candidates = @{$self->fetch_all_by_Slice($slice)};
    if ($respect_strand) { @candidates = grep {$_->strand == $ref_feature->strand} @candidates }
    if ($opposite_strand) { @candidates = grep {$_->strand != $ref_feature->strand} @candidates }
    # When user chooses UPSTREAM/DOWNSTREAM && FIVE_PRIME/THREE_PRIME, discount any features with
    # their other end in the overlap region. fetch_all_by_Slice cannot filter these out.
    # fetch_all_by_Slice_constraint could in principle with the correct constraint.
    if (($upstream && $ref_feature->strand == 1) || ($downstream && $ref_feature->strand == -1)) {
      if ($five_prime) {
        @candidates = grep { ($_->strand == 1) ? $_->seq_region_start < $ref_feature->start : $_->seq_region_end < $ref_feature->start} @candidates;
      } elsif ($three_prime) {
        @candidates = grep { ($_->strand == 1) ? $_->seq_region_end < $ref_feature->start : $_->seq_region_start < $ref_feature->start} @candidates;
      }
    }
    if (($downstream && $ref_feature->strand == 1)|| ($upstream && $ref_feature->strand == -1)) {
      if ($five_prime) {
        @candidates = grep { ($_->strand == 1) ? $_->seq_region_start > $ref_feature->end : $_->seq_region_end > $ref_feature->end } @candidates;
      } elsif ($three_prime) {
        @candidates = grep { ($_->strand == 1) ? $_->seq_region_end > $ref_feature->end : $_->seq_region_start > $ref_feature->end } @candidates;
      }
    }
    # Then sort and prioritise the candidates
    my $finalists; # = [[feature, distance, centre-weighted distance, length, dbID],..]
    $finalists = $self->select_nearest($ref_feature,\@candidates,$limit,$not_overlapping,$five_prime,$three_prime);
    $finalists = [ map { [ splice @$_,0,2 ]} @$finalists ]; # Remove the ugly bits from the sight of users.
    return $finalists;
}


=head2 select_nearest

  Arg [1]    : Bio::Ensembl::Feature, a Feature to find the nearest neighbouring feature to.
  Arg [2]    : Listref of Features to be considered for nearness.
  Arg [3]    : Integer, limited number of Features to return. Equally near features are all returned in spite of this limit
  Arg [4]    : Boolean, Overlapping prohibition. Overlapped Features are forgotten
  Arg [5]    : Boolean, use the 5' ends of the nearby features for distance calculation
  Arg [6]    : Boolean, use the 3' ends of the nearby features for distance calculation
  Example    : $feature_list = $feature_adaptor->select_nearest($ref_feature,\@candidates,$limit,$not_overlapping)
  Description: Take a list of possible features, and determine which is nearest. Nearness is a
               tricky concept. Beware of using the distance between Features, as it may not be the number you think
               it should be.
  Returntype : listref of Features ordered by proximity
  Caller     : BaseFeatureAdaptor->fetch_all_nearest_by_Feature

=cut

sub select_nearest {
    my $self = shift;
    my $ref_feature = shift;
    my $candidates = shift;
    my $limit = shift;
    my $not_overlapping = shift;
    my $five_prime = shift;
    my $three_prime = shift;
  
    # Convert circular coordinates to linear ones for distance calculation
    my $ref_start = ($ref_feature->start < $ref_feature->end) ? $ref_feature->start : $ref_feature->start - $ref_feature->length; # Not ->end, in case circular
    my $ref_end = $ref_feature->end;
    my $ref_midpoint = $self->_compute_midpoint($ref_feature);

    my $position_matrix = [];
    my $shortest_distance;
    my $adjusted_distance; # for when features overlap and we need another metric to order by

    foreach my $neighbour (@$candidates) {
        # rearrange coordinates into forward-stranded orientation, plus flatten circular features.
        my $neigh_start = ($neighbour->seq_region_start < $neighbour->seq_region_end) 
          ? $neighbour->seq_region_start : $neighbour->seq_region_start - $neighbour->length; # Not ->end, in case it is circular
        my $neigh_midpoint = $self->_compute_midpoint($neighbour);
        my $neigh_end = $neighbour->seq_region_end;

        # discard overlaps early if not required.
        next if ( $not_overlapping
          && (
               ( $neigh_start >= $ref_start && $neigh_end <= $ref_end )
            || ( $neigh_end <= $ref_end && $neigh_end >= $ref_start )
          )
        );
        my @args = ($ref_start,$ref_midpoint,$ref_end,$neigh_start,$neigh_midpoint,$neigh_end,$neighbour->strand);
        if ($five_prime || $three_prime) {
          my $five = $five_prime ? 1 : 0; # Choose 5' or 3' end
          ($shortest_distance,$adjusted_distance) = $self->_compute_prime_distance($five,@args);
          next unless ($shortest_distance || $adjusted_distance);
        }
        else {
          ($shortest_distance,$adjusted_distance) = $self->_compute_nearest_end(@args);
        }
        push @$position_matrix,[ $neighbour, $shortest_distance, $adjusted_distance, $neighbour->length, $neighbour->display_id] unless ($not_overlapping && $shortest_distance == 0);
    }
    
    # Order by distance, then centre-to-centre distance, then smallest feature first, then an arbitrary ID.
    # $position_matrix looks like this:
    # [ [ $feature, closest measure of distance, size, dbID ] ] 
    # Not proud of this line. It re-emits the original array in sorted order to favour better measures of order.
    my @ordered_matrix = map { [$_->[0],$_->[1],$_->[2],$_->[3],$_->[4]] } 
      sort { abs($a->[1]) <=> abs($b->[1]) || abs($a->[2]) <=> abs($b->[2]) || $a->[3] <=> $b->[3] || $a->[4] <=> $b->[4]} @$position_matrix;
    # Knock off unwanted hits. Equal distances are included, so more can be returned than asked for.
    @ordered_matrix = $self->_discard_excess_features_from_matrix(\@ordered_matrix,$limit);
    return \@ordered_matrix;
}

=head2 _compute_nearest_end

  Arg [1]    : Reference feature start
  Arg [2]    : Reference feature mid-point
  Arg [3]    : Reference feature end
  Arg [4]    : Considered feature start
  Arg [5]    : Considered feature mid-point
  Arg [6]    : Considered feature end
  Example    : $distance = $feature_adaptor->_compute_nearest_end($ref_start,$ref_midpoint,$ref_end,$f_start,$f_midpoint,$f_end)
  Description: For a given feature, calculate the smallest legitimate distance to a reference feature
               Calculate by mid-points to accommodate overlaps
  Returntype : Integer distance in base pairs
  Caller     : BaseFeatureAdaptor->select_nearest()

=cut
sub _compute_nearest_end {
    my ($self,$ref_start,$ref_midpoint,$ref_end,$f_start,$f_midpoint,$f_end) = @_;
    my ($effective_distance,$weighted_distance);    
    if ($ref_start > $f_end ) {
      # upstream, no overlap
      $effective_distance = $f_end - $ref_start; # negative result for upstream
    } elsif ( $ref_end < $f_start ) {
      # downstream
      $effective_distance = $f_start - $ref_end; # positive result for downstream
    } else {
      # overlap,
      $effective_distance = 0;
      $weighted_distance = $f_midpoint - $ref_midpoint;
    }
    return $effective_distance, $weighted_distance;
}

=head2 _compute_prime_distance

  Arg [1]    : Reference feature start
  Arg [2]    : Reference feature mid-point
  Arg [3]    : Reference feature end
  Arg [4]    : Considered feature start
  Arg [5]    : Considered feature mid-point
  Arg [6]    : Considered feature end
  Arg [7]    : Considered feature strand
  Example    : $distance,$weighted_centre_distance = $feature_adaptor->_compute_prime_distance($ref_start,$ref_midpoint,$ref_end,$f_start,$f_midpoint,$f_end,$f_strand)
  Description: Calculate the smallest distance to the 5' end of the considered feature
  Returntype : Integer distance in base pairs or a string warning that the result doesn't mean anything.
               Nearest 5' and 3' features shouldn't reside inside the reference Feature
  Caller     : BaseFeatureAdaptor->select_nearest()

=cut
sub _compute_prime_distance {
    my ($self,$five_prime,$ref_start,$ref_midpoint,$ref_end,$f_start,$f_midpoint,$f_end,$f_strand) = @_;
    
    my $f_coord;
    my $effective_distance;
    my $weighted_distance;

    if ($five_prime) { 
      $f_coord = ($f_strand == 1) ? $f_start : $f_end;
    } else { # three prime end required
      $f_coord = ($f_strand ==1 ) ? $f_end : $f_start;
    }

    if ($ref_start > $f_coord) {
      # n' end of feature upstream of reference
      $effective_distance = $f_coord - $ref_start;
      if ($f_end < $ref_start) {
        $weighted_distance = $f_coord - $ref_midpoint;
      }
    } elsif ($ref_end < $f_coord && $f_start > $ref_end) {
      # n' end of feature downstream of reference
      $effective_distance = $ref_end - $f_coord;
      if ($f_start > $ref_end) {
        $weighted_distance = $f_coord - $ref_midpoint;
      }
    } elsif ($ref_start < $f_coord && $ref_end > $f_coord) {
      # total overlap condition
      $effective_distance = 0;
      $weighted_distance = $f_coord - $ref_midpoint;
    }
    # If feature's end falls outside of the reference feature, compute distance to reference midpoint.
    return $effective_distance, $weighted_distance;
}

=head2 _compute_midpoint

  Arg [1]    : Bio::EnsEMBL::Feature
  Example    : $middle = $feature_adaptor->_compute_midpoint($feature);
  Description: Calculate the mid-point of a Feature. Used for comparing Features that overlap each other
               and determining a canonical distance between two Features for the majority of use cases.
  Returntype : Integer coordinate rounded down.
  Caller     : BaseFeatureAdaptor->select_nearest()

=cut
sub _compute_midpoint {
    my $self = shift;
    my $feature = shift;
    my $start = $feature->seq_region_start;
    my $end = $feature->seq_region_end;
    my $midpoint;
    # consider circular slice
    if ($start > $end) {
        $midpoint = int($end + ($start - $end) / 2);
    } else {
        $midpoint = int($start + ($end - $start) / 2);
    }
    return $midpoint;
}

sub _discard_excess_features_from_matrix {
  my $self = shift;
  my $list = shift;
  my @ordered_matrix = @$list;
  my $limit = shift;
  return @ordered_matrix if $#ordered_matrix == 0 || !defined($limit) || $limit > scalar @ordered_matrix;
  # cut off excess elements
  my @spares = splice @ordered_matrix, $limit, scalar @ordered_matrix - $limit;
  # Check nearest distance of the last member of ordered_matrix against spares and include those that match
  # This prevents first past the post.
  if (scalar @ordered_matrix > 0) {
    my $threshold_distance = $ordered_matrix[-1]->[1];
    my $i = 0;
    while ($i <= $#spares && $spares[$i]->[1] == $threshold_distance) {
      # printf "Considering option %s, %s\n",$spares[$i]->[0]->stable_id,$spares[$i]->[1];
      push @ordered_matrix, $spares[$i];
      $i++;
    }
    
  }
  return @ordered_matrix;
}

1;
