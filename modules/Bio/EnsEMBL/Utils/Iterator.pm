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

  Argument : a coderef representing the iterator, this anonymous subroutine
             is assumed to return the next object in the set when called,
             and to return undef when the set is exhausted. If the argument
             is not defined then we return an 'empty' iterator that immediately
             returns undef

  Example :

    my @dbIDs = fetch_relevant_dbIDs();

    my $iterator = Bio::EnsEMBL::Utils::Iterator->new(
        sub { return $self->fetch_by_dbID(shift @dbIDs) }
    );

    NB: this is a very simple example showing how to call the constructor
    that would be rather inefficient in practice, real examples should 
    probably be smarter about batching up queries to minimise trips to
    the database. See examples in the Variation API.

  Description: Constructor, creates a new iterator object
  Returntype : Bio::EnsEMBL::Utils::Iterator instance
  Exceptions : thrown if the supplied argument is not a coderef
  Caller     : general
  Status     : Experimental

=cut

sub new {
    my $class = shift;

    my $coderef = shift;
    
    # if the user doesn't supply a coderef, we create a
    # simple 'empty' iterator that immediately returns undef
    if (not defined $coderef) {
        $coderef = sub { return undef };
    }
    elsif (ref $coderef ne 'CODE'){
        throw("The supplied argument does not look like an coderef")
    }

    my $self = {sub => $coderef};

    return bless $self, $class;
}

=head2 next

  Example    : $obj = $iterator->next()
  Description: returns the next object from this iterator, or undef if the iterator is exhausted
  Returntype : object reference (the type will depend on what this iterator is iterating over)
  Exceptions : none
  Caller     : general
  Status     : Experimental

=cut

sub next {
    my $self = shift;

    # if someone has called has_next or peek, there might be a cached value we can return
    $self->{next} ||= $self->{sub}->();
    
    return delete $self->{next};
}

=head2 has_next

  Example    : if ($iterator->has_next) { my $obj = $iterator->next }
  Description: returns true if this iterator has more objects to fetch, false when it is exhausted
  Returntype : boolean
  Exceptions : none
  Caller     : general
  Status     : Experimental

=cut

sub has_next {
    my $self = shift;

    $self->{next} ||= $self->{sub}->();

    return defined $self->{next}; 
}

=head2 peek

  Example    : $obj = $iterator->peek()
  Description: returns the next object from this iterator, or undef if the iterator is exhausted,
               much like next() but does not advance the iterator (so the same object will be 
               returned on the following call to next() or peek())
  Returntype : object reference (the type will depend on what this iterator is iterating over)
  Exceptions : none
  Caller     : general
  Status     : Experimental

=cut

sub peek {
    my $self = shift;

    $self->{next} ||= $self->{sub}->();

    return $self->{next};
}

=head2 grep

  Example    : my $filtered_iterator = $original_iterator->grep(sub {$_->name =~ /^rs/});
  Description: filter this iterator, returning another iterator
  Argument   : a coderef which returns true if the element should be included in the
               filtered set, or false if the element should be filtered out. $_ will be 
               set locally to each element in turn so you should be able to write a block 
               in a similar way as for the perl grep function (although it will need to be 
               preceded with the sub keyword). Otherwise you can pass in a reference to an 
               existing subroutine with the same behaviour.
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
            return $_ if $coderef->();
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
               Otherwise you can pass in a reference to an existing subroutine with 
               the same behaviour.
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
        return defined $_ ? $coderef->() : undef;
    });
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

    while (my $obj = $self->next) {
        push @array, $obj;
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

1;
