package Bio::EnsEMBL::Utils::Iterator;

=head1 LICENSE

 Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

  Bio::EnsEMBL::Utils::Iterator

=head1 SYNOPSIS

  my $variation_iterator =
    $variation_adaptor->fetch_Iterator_by_VariationSet($1kg_set);

  while ( my $variation = $variation_iterator->next ) {
    # operate on variation object
    print $variation->name, "\n";
  }

=head1 DESCRIPTION
  
  Some adaptor methods may return more objects than can fit in memory at once, in these cases 
  you can fetch an iterator object instead of the usual array reference. The iterator object 
  allows you to iterate over the set of objects (using the next() method) without loading the
  entire set into memory at once. You can tell if an iterator is exhausted with the has_next()
  method. The peek() method allows you to fetch the next object from the iterator without
  advancing the iterator - this is useful if you want to check some property of en element in 
  the set while leaving the iterator unchanged.

  You can filter and transform an iterator in an analogous way to using map and grep on arrays 
  using the provided map() and grep() methods. These methods return another iterator, and only 
  perform the filtering and transformation on each element as it is requested, so again these 
  can be used without loading the entire set into memory.

  Iterators can be combined together with the append() method which merges together the 
  iterator it is called on with the list of iterators passed in as arguments. This is 
  somewhat analogous to concatenating arrays with the push function. append() returns a new 
  iterator which iterates over each component iterator until it is exhausted before moving
  on to the next iterator, in the order in which they are supplied to the method.

  An iterator can be converted to an array (reference) containing all the elements in the 
  set with the to_arrayref() method, but note that this array may consume a lot of memory if 
  the set the iterator is iterating over is large and it is recommended that you do not call 
  this method unless there is no way of working with each element at a time.

=head1 METHODS

=cut

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw);

=head2 new

  Args [1]:  Either:
                a) A CODE reference representing the Iterator, or a Slice
                   fetch method to be iterated over using sub Slices.
                b) An ARRAY reference to be iterated over
             If the argument is not defined then we return an 'empty' 
             iterator that immediately returns undef

  Args [2]: Optional: Feature Adaptor which was the source of the Slice 
            fetch CODE ref

  Args [3]: Optional: Slice fetch method params ARRAY ref
  
  Args [4]: Optional: Slice parameter index in param arrays. Default=0

  Args [5]: Optional: sub Slice chunk size. Default=500000

  Example :

    dbID Iterator:

    my @dbIDs = fetch_relevant_dbIDs();

    my $iterator = Bio::EnsEMBL::Utils::Iterator->new(
        sub { return $self->fetch_by_dbID(shift @dbIDs) }
    );

    NB: this is a very simple example showing how to call the constructor
    that would be rather inefficient in practice, real examples should 
    probably be smarter about batching up queries to minimise trips to
    the database. See examples in the Variation API.

    Slice Iterator:

    my $iter = Bio::EnsEMBL::Utils::Iterator->new
	                             ($feat_adaptor->can('fetch_all_by_Slice'),
	                              $feat_adaptor,
	                              $fetch_all_by_Slice_params_ref,
	                              0,#Slice idx
 	                              #500 #chunk length
	                             );

  Description: Constructor, creates a new iterator object
  Returntype : Bio::EnsEMBL::Utils::Iterator instance
  Exceptions : thrown if the supplied argument is not the expected 
  Caller     : general
  Status     : at risk

=cut

sub new{
    my ($class, $ref,  $adaptor, $params_ref, $slice_param_idx, $chunk_size) = @_;

    my $coderef;

    if (! defined $ref) {
	  # if the user doesn't supply an argument, we create a
	  # simple 'empty' iterator that immediately returns undef
   	  $coderef = sub { return undef };
    } 
	elsif (ref($ref) eq 'ARRAY') {
	  # if the user supplies an arrayref as an argument, we 
	  # create an iterator over this array
 	  $coderef = sub { return shift @$ref };
    } 
	elsif (ref($ref) eq 'CODE') {
	  
	  if (! defined $adaptor) {	   # Standard Iterator
		$coderef = $ref;
	  } 
	  else {					   # Slice chunk Iterator
	  
		if (! (ref($adaptor) &&
			   $adaptor->isa('Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor'))
		   ) {
		  throw('You must pass a valid Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor');
		}
		
		if (! ($params_ref && 
			   ref($params_ref) eq 'ARRAY')) {
		  #Don't need to check size here so long as we have valid Slice
		  throw('You must pass a method params ARRAYREF');
		}
	  
		$slice_param_idx = 0 if(! defined $slice_param_idx);
		my $slice = $params_ref->[$slice_param_idx];
		
		if (! (defined $slice &&
			   ref($slice) eq 'Bio::EnsEMBL::Slice')
		   ) {
		  throw('You must pass a valid Bio::EnsEMBL::Slice '.
				'in your method params and a valid slice param idx arg');
		}
		
		$chunk_size ||= 1000000;
		
		my @feat_cache;
		my $finished     = 0;
		my $start        = 1;	#local coord for sub slice
		my $end          = $slice->length;
		my $num_overlaps = 0;
		
		$coderef = 
		  sub {
			
			while (scalar(@feat_cache) == 0 &&
				   ! $finished) {
			  
			  my $new_end = $start + $chunk_size;
			
			  if ($new_end >= $end) {
				# this is our last chunk
				$new_end = $end;
				$finished = 1;  
			  }
			 
			  #Chunk by sub slicing
			  my $sub_slice                   = $slice->sub_Slice($start, $new_end);
			  $params_ref->[$slice_param_idx] = $sub_slice;
			  @feat_cache = @{ $ref->($adaptor, @$params_ref)};
		
			  
			  #Remove & count overlapping features
			  splice(@feat_cache, 0, $num_overlaps) if($num_overlaps);
			  my $i;
		
			  if (scalar(@feat_cache) > 0) {

				my $feat_end  = $feat_cache[$#feat_cache]->end;
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
	  }
	} else {
	  throw("The supplied argument does not look like an arrayref or a coderef ".ref($ref));
	}
	
	return bless {sub => $coderef}, $class;
}

=head2 next

  Example    : $obj = $iterator->next
  Description: returns the next object from this iterator, or undef if the iterator is exhausted
  Returntype : the return type will depend on what this iterator is iterating over
  Exceptions : none
  Caller     : general
  Status     : Experimental

=cut

sub next {
    my $self = shift;

    $self->{next} = $self->{sub}->() unless defined $self->{next};
    
    return delete $self->{next};
}

=head2 has_next

  Example    : if ($iterator->has_next) { my $obj = $iterator->next }
  Description: returns true if this iterator has more elements to fetch, false when
               it is exhausted
  Returntype : boolean
  Exceptions : none
  Caller     : general
  Status     : Experimental

=cut

sub has_next {
    my $self = shift;

    $self->{next} = $self->{sub}->() unless defined $self->{next};

    return defined $self->{next}; 
}

=head2 peek

  Example    : $obj = $iterator->peek
  Description: returns the next object from this iterator, or undef if the iterator is exhausted,
               much like next but does not advance the iterator (so the same object will be 
               returned on the following call to next or peek)
  Returntype : the return type will depend on what this iterator is iterating over
  Exceptions : none
  Caller     : general
  Status     : Experimental

=cut

sub peek {
    my $self = shift;

    $self->{next} = $self->{sub}->() unless defined $self->{next};

    return $self->{next};
}

=head2 grep

  Example    : my $filtered_iterator = $original_iterator->grep(sub {$_->name =~ /^rs/});
  Description: filter this iterator, returning another iterator
  Argument   : a coderef which returns true if the element should be included in the
               filtered set, or false if the element should be filtered out. $_ will be 
               set locally to each element in turn so you should be able to write a block 
               in a similar way as for the perl grep function (although it will need to be 
               preceded with the sub keyword). Otherwise you can pass in a reference to a 
               subroutine which expects a single argument with the same behaviour.
  Returntype : Bio::EnsEMBL::Utils::Iterator
  Exceptions : thrown if the argument is not a coderef
  Caller     : general
  Status     : Experimental

=cut

sub grep {
    my ($self, $coderef) = @_;

    throw('Argument should be a coderef') unless ref $coderef eq 'CODE';

    return Bio::EnsEMBL::Utils::Iterator->new(sub {
        while ($self->has_next) {
            local $_ = $self->next;
            return $_ if $coderef->($_);
        }
        return undef;
    });
}

=head2 map

  Example    : my $transformed_iterator = $original_iterator->map(sub {$_->name});
  Description: transform the elements of this iterator, returning another iterator
  Argument   : a coderef which returns the desired transformation of each element.
               $_ will be set locally set to each original element in turn so you 
               should be able to write a block in a similar way as for the perl map 
               function (although it will need to be preceded with the sub keyword). 
               Otherwise you can pass in a reference to a subroutine which expects a
               single argument with the same behaviour.
  Returntype : Bio::EnsEMBL::Utils::Iterator
  Exceptions : thrown if the argument is not a coderef
  Caller     : general
  Status     : Experimental

=cut


sub map {
    my ($self, $coderef) = @_;
    
    throw('Argument should be a coderef') unless ref $coderef eq 'CODE';

    return Bio::EnsEMBL::Utils::Iterator->new(sub {
        local $_ = $self->next;
        return defined $_ ? $coderef->($_) : undef;
    });
}


=head2 each

  Example    : $iterator->each(sub { print $_->name, "\n"; });
  Description: Performs a full iteration of the current iterator instance.
  Argument   : a coderef which returns the desired transformation of each element.
               $_ will be set locally set to each element.
  Returntype : None
  Exceptions : thrown if the argument is not a coderef
  Caller     : general
  Status     : Experimental

=cut


sub each {
    my ($self, $coderef) = @_;
    throw('Argument should be a coderef') unless ref $coderef eq 'CODE';
    while($self->has_next()) {
        local $_ = $self->next();
        $coderef->($_);
    }
    return;
}


=head2 to_arrayref

  Example    : my $arrayref = $iterator->to_arrayref;
  Description: return a reference to an array containing all elements from the 
               iterator. This is created by simply iterating over the iterator 
               until it is exhausted and adding each element in turn to an array. 
               Note that this may consume a lot of memory for iterators over 
               large collections
  Returntype : arrayref
  Exceptions : none
  Caller     : general
  Status     : Experimental

=cut

sub to_arrayref {
    my ($self) = @_;
    
    my @array;

    while ($self->has_next) {
        push @array, $self->next;
    }

    return \@array;
}

=head2 append

  Example    : my $combined_iterator = $iterator1->append($iterator2, $iterator3);
  Description: return a new iterator that combines this iterator with the others
               passed as arguments, this new iterator will iterate over each
               component iterator (in the order supplied here) until it is 
               exhausted and then move on to the next iterator until all are
               exhausted
  Argument   : an array of Bio::EnsEMBL::Utils::Iterator objects
  Returntype : Bio::EnsEMBL::Utils::Iterator
  Exceptions : thrown if any of the arguments are not iterators
  Caller     : general
  Status     : Experimental

=cut

sub append {
    my ($self, @queue) = @_;

    for my $iterator (@queue) {
        throw("Argument to append doesn't look like an iterator")
            unless UNIVERSAL::can($iterator, 'has_next') && UNIVERSAL::can($iterator, 'next');
    }

    # push ourselves onto the front of the queue
    unshift @queue, $self;

    return Bio::EnsEMBL::Utils::Iterator->new(sub {
        # shift off any exhausted iterators
        while (@queue && not $queue[0]->has_next) {
            shift @queue;
        }

        # and return the next object from the iterator at the 
        # head of the queue, or undef if the queue is empty
        return @queue ? $queue[0]->next : undef;
    });
}

=head2 take

  Example    : my $limited_iterator = $iterator->take(5);
  Description: return a new iterator that only iterates over the
               first n elements of this iterator
  Argument   : a positive integer
  Returntype : Bio::EnsEMBL::Utils::Iterator
  Exceptions : thrown if the argument is negative
  Caller     : general
  Status     : Experimental

=cut

sub take {
    my ($self, $n) = @_;
    
    throw("Argument cannot be negative") if $n < 0;

    my $cnt = 0;

    return Bio::EnsEMBL::Utils::Iterator->new(sub {
        return $cnt++ >= $n ? undef : $self->next;
    });
}

=head2 skip

  Example    : my $limited_iterator = $iterator->skip(5);
  Description: skip over the first n elements of this iterator (and then return
               the same iterator for your method chaining convenience)
  Argument   : a positive integer
  Returntype : Bio::EnsEMBL::Utils::Iterator
  Exceptions : thrown if the argument is negative
  Caller     : general
  Status     : Experimental

=cut

sub skip {
    my ($self, $n) = @_;
    
    throw("Argument cannot be negative") if $n < 0;

    $self->next for (0 .. $n-1);

    return $self;
}

=head2 reduce

  Example    : my $tot_length = $iterator->reduce(sub { $_[0] + $_[1]->length }, 0);
  Description: reduce this iterator with the provided coderef, using the (optional)
               second argument as the initial value of the accumulator
  Argument[1]: a coderef that expects 2 arguments, the current accumulator
               value and the next element in the set, and returns the next
               accumulator value. Unless the optional second argument is
               provided the first accumulator value passed in will be the 
               first element in the set
  Argument[2]: (optional) an initial value to use for the accumulator instead 
               of the first value of the set
  Returntype : the same as the return type of the coderef 
  Exceptions : thrown if the argument is not a coderef
  Caller     : general
  Status     : Experimental

=cut

sub reduce {
    my ($self, $coderef, $init_val) = @_;

    throw('Argument should be a coderef') unless ref $coderef eq 'CODE';

    my $result = defined $init_val ? $init_val : $self->next;

    while ($self->has_next) { 
        $result = $coderef->($result, $self->next);
    }

    return $result;
}

1;

