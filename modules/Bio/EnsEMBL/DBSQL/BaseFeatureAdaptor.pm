=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

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
			
			if ($feat_end > $slice_end) {
			  $feat_end  = $feat_cache[$i]->end;
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

  Arg [3]    : string $logic_name
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
                  slice and sums them up
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
	
	#Constraints
	$constraint ||= '';
  $constraint = $self->_logic_name_to_constraint( $constraint, $logic_name );
  # If the logic name was invalid, undef was returned
  return $count if ! defined $constraint;
  
  #Query logic
  my $sa = $slice->adaptor();
  my $proj_ref = $self->_get_and_filter_Slice_projections($slice);

  my $proj_ref_length = @{$proj_ref};
  #Manual loop to support look-ahead/behind
  for (my $i = 0; $i < $proj_ref_length; $i++) {
    my $seg = $proj_ref->[$i];
    my $seg_slice = $seg->to_Slice();
    $self->_bind_param_generic_fetch($bind_params);
    
    # Because we cannot filter boundary crossing features in code we need to
    # do it in SQL. So we detect when we are not on the original query Slice
    # we do manual filtering by the seg_Slice's start and end. If we are on
    # the *1st* section we only limit by the *end* and when we are on the *last*
    # we filter by the *start*
    #
    # This is the same as fetch_all_by_Slice_constraint()'s in-memory filtering
    # except we need to alter projected features in that code with an offset
    my $local_constraint = $constraint;
    if ( $seg_slice->name() ne $slice->name() ) {
      my ($limit_start, $limit_end) = (1,1); #limit both by default
      if($i == 0) {
        $limit_start = 1; #don't check start as we are on the final projection
      }
      elsif($i == ($proj_ref_length - 1)) {
        $limit_end = 1; #don't check end as we are on the final projection
      }
      $local_constraint .= ' AND ' if $local_constraint;
      my @conditionals;
      
      #Do not cross the start boundary so our feature must be less than slice end on all counts
      if($limit_start) {
        push(@conditionals, sprintf('%1$s.seq_region_start <= %2$d AND %1$s.seq_region_end <= %2$d', $table_synonym, $seg_slice->end));
      }
      #Do not cross the start boundary so our feature must be larger than slice start on all counts
      if($limit_end) {
        push(@conditionals, sprintf('%1$s.seq_region_start >= %2$d AND %1$s.seq_region_end >= %2$d', $table_synonym, $seg_slice->start));
      }
      
      $local_constraint .= join(q{ AND }, @conditionals);
    }
    
  	my $count_array = $self->_get_by_Slice($seg_slice, $local_constraint, 'count');
  	#Data comes out as an array
  	$count += $_ for @$count_array;
  }
	
	return $count;
}

=head2 _get_and_filter_Slice_projections

    Arg [1]     : Bio::EnsEMBL::Slice
    Description : Finds all features with at least partial overlap to the given
                  slice and sums them up
    Example     : my $proj= $self->_get_and_filter_Slice_projections($slice);
    Returntype  : ArrayRef Bio::EnsEMBL::ProjectionSegment; Returns an array
                  of projected segments
=cut

sub _get_and_filter_Slice_projections {
  my ($self, $slice) = @_;
  my $sa = $slice->adaptor();
  my @proj = @{ $sa->fetch_normalized_slice_projection($slice) };
  if ( !@proj ) {
    throw( 'Could not retrieve normalized Slices. '
         . 'Database contains incorrect assembly_exception information.'
    );
  }
  
  # Want to get features on the FULL original slice as well as any
  # symlinked slices.
  
  # Filter out partial slices from projection that are on same
  # seq_region as original slice.

  my $sr_id = $slice->get_seq_region_id();

  @proj = grep { $_->to_Slice->get_seq_region_id() != $sr_id } @proj;

  my $segment = bless( [ 1, $slice->length(), $slice ],
                       'Bio::EnsEMBL::ProjectionSegment' );
  push( @proj, $segment );
  return \@proj;
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
                            . " OR ".$table_synonym.".seq_region_start > ".$table_synonym.".seq_region_end";
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
  my $self = shift;
  my $start = shift;
  my $end   = shift;
  my $strand = shift;
  my $slice = shift;

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
  if($end < $start && !$slice->is_circular()) {
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

  my ($seq_region, $start, $end, $strand);

  my $slice_seq_region_id = $slice->get_seq_region_id();
  my $slice_seq_region = $slice->seq_region_name();

  foreach my $f (@$features) {
    #since feats were obtained in contig coords, attached seq is a contig
    my $fslice = $f->slice();
    if(!$fslice) {
      throw("Feature does not have attached slice.\n");
    }
    my $fseq_region = $fslice->seq_region_name();
    my $fseq_region_id = $fslice->get_seq_region_id();
    my $fcs = $fslice->coord_system();

    if(!$slice_cs->equals($fcs)) {
      #slice of feature in different coord system, mapping required

      ($seq_region, $start, $end, $strand) =
        $mapper->fastmap($fseq_region_id,$f->start(),$f->end(),$f->strand(),$fcs);

      # undefined start means gap
      next if(!defined $start);
    } else {
      $start      = $f->start();
      $end        = $f->end();
      $strand     = $f->strand();
      $seq_region = $f->slice->seq_region_name();
    }
    
    # maps to region outside desired area
    next if ($start > $slice_end) || ($end < $slice_start) || 
      ($slice_seq_region ne $seq_region);

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
# Given a logic name and an existing constraint this will
# add an analysis table constraint to the feature.  Note that if no
# analysis_id exists in the columns of the primary table then no
# constraint is added at all
#
sub _logic_name_to_constraint {
  my $self = shift;
  my $constraint = shift;
  my $logic_name = shift;

  return $constraint if(!$logic_name);

  #make sure that an analysis_id exists in the primary table
  my ($prim_tab) = $self->_tables();
  my $prim_synonym = $prim_tab->[1];

  my $found_analysis=0;
  foreach my $col ($self->_columns) {
    my ($syn,$col_name) = split(/\./,$col);
    next if($syn ne $prim_synonym);
    if($col_name eq 'analysis_id') {
      $found_analysis = 1;
      last;
    }
  }

  if(!$found_analysis) {
    warning("This feature is not associated with an analysis.\n" .
            "Ignoring logic_name argument = [$logic_name].\n");
    return $constraint;
  }

  my $aa = $self->db->get_AnalysisAdaptor();
  my $an = $aa->fetch_by_logic_name($logic_name);

  if ( !defined($an) ) {
    return undef;
  }

  my $an_id = $an->dbID();

  $constraint .= ' AND' if($constraint);
  $constraint .= " ${prim_synonym}.analysis_id = $an_id";
  return $constraint;
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


=head1 DEPRECATED METHODS

=cut


=head2 fetch_all_by_RawContig_constraint

  Description: DEPRECATED use fetch_all_by_RawContig_constraint instead

=cut

sub fetch_all_by_RawContig_constraint {
  my $self = shift;
  deprecate('Use fetch_all_by_Slice_constraint() instead.');
  return $self->fetch_all_by_slice_constraint(@_);
}

=head2 fetch_all_by_RawContig

  Description: DEPRECATED use fetch_all_by_Slice instead

=cut

sub fetch_all_by_RawContig {
  my $self = shift;
  deprecate('Use fetch_all_by_Slice() instead.');
  return $self->fetch_all_by_Slice(@_);
}

=head2 fetch_all_by_RawContig_and_score

  Description: DEPRECATED use fetch_all_by_Slice_and_score instead

=cut

sub fetch_all_by_RawContig_and_score{
  my $self = shift;
  deprecate('Use fetch_all_by_Slice_and_score() instead.');
  return $self->fetch_all_by_Slice_and_score(@_);
}

=head2 remove_by_RawContig

  Description: DEPRECATED use remove_by_Slice instead

=cut

sub remove_by_RawContig {
  my $self = shift;
  deprecate("Use remove_by_Slice instead");
  return $self->remove_by_Slice(@_);
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


1;
