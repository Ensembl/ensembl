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
    $variation_adaptor->fetch_iterator_by_VariationSet($1kg_set);

  while ( my $variation = $variation_iterator->next ) {
    # operate on variation object
    print $variation->name, "\n";
  }


=head1 DESCRIPTION
  
  Some adaptor methods may return more objects than can fit in memory at once, in these cases 
  you can fetch an iterator object instead of the usual list reference. The iterator object 
  allows you to iterate over the set of objects (using the next() method) without loading the
  entire set into memory at once. You can tell if an iterator is exhausted with the has_next()
  method.

  You can map and grep an iterator in an analogous way to using map and grep on arrays using
  the provided map and grep methods. These methods return another iterator, and only perform
  the filtering and transformation on each element as it is requested, so again these can be 
  used without loading the entire set into memory.

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Iterator;

use strict;
use warnings;

=head2 new

  Argument : a coderef representing the iterator, this anonymous subroutine
          is assumed to return the next object in the set when called,
          and to return undef when the set is exhausted. If the argument
          is not defined then we return an 'empty' iterator that immediately
          returns undef

  Example :

    my $iterator = Bio::EnsEMBL::Utils::Iterator->new(
        sub { return $self->fetch_by_dbID(shift @dbIDs) }
    );

    NB: this is a very simple example showing how to call the constructor
    that would be rather inefficient in practice, real examples should 
    probably be smarter about batching up queries to minimise trips to
    the database. See examples in the Variation API.

  Description: Constructor, creates a new iterator object
  Returntype : Bio::EnsEMBL::Utils::Iterator instance
  Exceptions : dies if the supplied argument is not a coderef
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
    else {
        die "The supplied argument does not look like an coderef"
            unless ref $coderef eq 'CODE';
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

    # if someone has called has_next, there might be a cached value we can return

    if ($self->{next}) {
        return delete $self->{next};
    }
    else {
        return $self->{sub}->();
    }
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

    $self->{next} = $self->{sub}->();

    return defined $self->{next}; 
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
  Exceptions : dies if the argument is not a coderef
  Caller     : general
  Status     : Experimental

=cut

sub grep {
    my ($self, $coderef) = @_;

    die "Argument should be a coderef" unless ref $coderef eq 'CODE';

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
  Exceptions : dies if the argument is not a coderef
  Caller     : general
  Status     : Experimental

=cut

sub map {
    my ($self, $coderef) = @_;
    
    die "Argument should be a coderef" unless ref $coderef eq 'CODE';

    return Bio::EnsEMBL::Utils::Iterator->new(sub {
        local $_ = $self->next;
        return defined $_ ? $coderef->() : undef;
    });
}

1;
