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

Bio::EnsEMBL::DBSQL::Support::BaseCache - Base cache code 

=head1 SYNOPSIS

  package Cache;
  
  use base qw/Bio::EnsEMBL::DBSQL::Support::BaseCache/;
  
  sub build_cache {
    return {}; #sends back a very basic cache which is a hash
  }
  
  1;
  
  #In use
  $cache->put(1, 'a');
  is($cache->get(1), 'a');
  is($cache->put(1, 'b'), 'a'); #put returns the object it replaced
  is($cache->delete(1), 'b');   #delete returns the object it removed
  
  is_deeply([$cache->cache_keys()], [1]);
  
  $cache->clear_cache();
  
  is($cache->size(), 0);

=head1 DESCRIPTION

A base class used for holding methods common to all cache implemnetations. 
Never use this class to do direct caching instead use one of the following

=over 8

=item C<Bio::EnsEMBL::DBSQL::Support::LruIdCache>

=item C<Bio::EnsEMBL::DBSQL::Support::FullIdCache>

=back

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::Support::BaseCache;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw/throw/;
use Bio::EnsEMBL::Utils::Scalar qw/assert_ref/;
use Scalar::Util qw/weaken/;

=head2 new

  Arg [1]    : Bio::EnsEMBL::DBSQL::BaseAdaptor $db_adaptor
  Example    : $cache = new CacheInheritedFromBaseCache($db_adaptor);
  Description: Creates a new cache class which handles all the basics of
               working with a cache apart from what that cache implementation
               is (apart from a hash)
  Returntype : Bio::EnsEMBL::DBSQL::Support::BaseCache
  Exceptions : none
  Caller     : BaseAdaptors
  Status     : Beta

=cut

sub new {
  my ($class, $adaptor) = @_;
  
  $class = ref($class) || $class;
  my $self = bless({}, $class);
  
  throw "Need an adaptor instance to delegate calls to" unless $adaptor;
  $self->adaptor($adaptor);
  
  return $self;
}

=head2 adaptor

  Arg [1]    : Bio::EnsEMBL::DBSQL::BaseAdaptor $db_adaptor
  Description: Getter setter for the adaptor this serves as an ID cacher for
  Returntype : Bio::EnsEMBL::DBSQL::BaseAdaptor
  Exceptions : none
  Caller     : BaseAdaptors
  Status     : Beta

=cut


sub adaptor {
  my ($self, $adaptor) = @_;
  if(defined $adaptor) {
    assert_ref($adaptor, 'Bio::EnsEMBL::DBSQL::BaseAdaptor', 'adaptor');
  	$self->{'adaptor'} = $adaptor;
  	weaken($self->{'adaptor'});
  }
  return $self->{'adaptor'};
}

=head2 cache

  Description: Returns back a Hash implementing object and also calls 
               C<build_cache()> for an initialise on demand system
  Returntype : Hash
  Exceptions : none
  Caller     : BaseAdaptors and internal
  Status     : Beta

=cut

sub cache {
  my ($self) = @_;
  if(! defined $self->{cache}) {
    my $cache = $self->build_cache();
    $self->{cache} = $cache;
  }
  return $self->{cache};
}

=head2 delete_cache

  Example    : $cache->delete_cache();
  Description: Deletes the cache. Normally used to trigger a C<build_cache()> 
               call
  Returntype : none
  Exceptions : none
  Caller     : BaseAdaptors
  Status     : Beta

=cut

sub delete_cache {
  my ($self) = @_;
  delete $self->{cache};
  return;
}

=head2 get

  Arg [1]    : String key to retrieve
  Example    : $cache->put(1,'a'); is($cache->get(1), 'a');
  Description: Retrieves the value held in the cache. Can return undef if the
               value could not be found
  Returntype : Scalar value held in the cache or nothing
  Exceptions : If key was undefined
  Caller     : BaseAdaptors
  Status     : Beta

=cut

sub get {
  my ($self, $key) = @_;
  throw "No key given" unless defined $key;
  my $cache = $self->cache();
  if(exists $cache->{$key}) {
    return $cache->{$key};
  }
  return;
}

=head2 get_by_list

  Arg [1]    : ArrayRef keys to retrieve
  Example    : is($cache->get_by_list([1,2]), ['a','b']);
  Description: Retrieves the values held in the cache. If a key cannot be
               found you get no entry in the resulting array returned.
  Returntype : ArrayRef of found values
  Exceptions : None
  Caller     : BaseAdaptors
  Status     : Beta

=cut

sub get_by_list {
  my ($self, $list) = @_;
  assert_ref($list, 'ARRAY', 'list');
  my @output;
  my $cache = $self->cache();
  foreach my $key (@{$list}) {
    next if ! defined $key;
    if(exists $cache->{$key}) {
      push(@output, $cache->{$key});
    }
  }
  return \@output;
}

=head2 put

  Arg [1]    : String key to store
  Arg [2]    : Scalar value to store
  Example    : $cache->put(2, 'b');
  Description: Stores a value in the cache. Returns the previous value held 
               under that key if one existed
  Returntype : Scalar of the previously held value if one existed
  Exceptions : None
  Caller     : BaseAdaptors
  Status     : Beta

=cut

sub put {
  my ($self, $key, $new) = @_;
  throw "No key given" unless defined $key;
  my $cache = $self->cache();
  my $old;
  if(exists $cache->{$key}) {
    $old = $cache->{$key};
  }
  $cache->{$key} = $new;
  return $old if $old;
  return;
}

=head2 remove

  Arg [1]    : String key to remove
  Example    : is($cache->remove(1), 'a');
  Description: Removes the supplied key from the cache
  Returntype : Scalar value held in the cache or nothing
  Exceptions : If key was undefined
  Caller     : BaseAdaptors
  Status     : Beta

=cut

sub remove {
  my ($self, $key) = @_;
  throw "No key given" unless defined $key;
  my $cache = $self->cache();
  my $old;
  if(exists $cache->{$key}) {
    $old = $cache->{$key};
    delete $cache->{$key};
  }
  return $old if $old;
  return;
}

=head2 clear_cache

  Example    : $cache->clear_cache();
  Description: Removes all values from the cache but does not delete the cache
               instance
  Returntype : None
  Exceptions : None
  Caller     : BaseAdaptors
  Status     : Beta

=cut

sub clear_cache {
  my ($self) = @_;
  my $cache = $self->cache();
  %{$cache} = ();
  return;
}

=head2 cache_keys

  Example    : my @keys = $cache->cache_keys();
  Description: Returns all keys held by the cache
  Returntype : List of all available keys
  Exceptions : None
  Caller     : BaseAdaptors
  Status     : Beta

=cut

sub cache_keys {
  my ($self) = @_;
  return keys %{$self->cache()};
}

=head2 size

  Example    : $cache->size();
  Description: Returns the number of keys currrently held by the cache
  Returntype : Int $size
  Exceptions : None
  Caller     : BaseAdaptors
  Status     : Beta

=cut

sub size {
  my ($self) = @_;
  return scalar($self->cache_keys());
}

=head2 build_cache

  Description: Implement to return the required type of cache
  Returntype : Hash
  Exceptions : None
  Caller     : BaseAdaptors
  Status     : Beta

=cut

sub build_cache {
  my ($self) = @_;
  die 'Unimplemented';
}

1;