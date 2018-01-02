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

  #Try using SQL - cache will be consulted accordingly
  my $ids = $cache->get_by_sql('select dbid from table where val like ?', ['someval%']);

=head1 DESCRIPTION

A base class used for holding methods common to all cache implementations. 
Never use this class to do direct caching instead use one of the following

=over 8

=item C<Bio::EnsEMBL::DBSQL::Support::LruIdCache>

=item C<Bio::EnsEMBL::DBSQL::Support::FullIdCache>

=back

To provide exta functionality to the caches you should override one of the above
classes and extend. Caches work when you use inheritence by composition in their
target adaptor.

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
  Example    : $cache = CacheInheritedFromBaseCache->new($db_adaptor);
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
  Exceptions : If the given ArrayRef was undefined
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

=head2 get_by_sql

  Arg [1]    : String SQL to execute. Should return the key of this cache in column 1
  Arg [2]    : ArrayRef Parameters to bind to the specified query
  Example    : $cache->get_by_sql('select id from tables where p like ?', ['val%']);
  Description: Executes the given SQL statement against the construnction adaptor's
               backing database. The found IDs are then passed into C<get_by_list()>
               where the elements are returned should the cache hold them.

               Remember if the cache is un-aware of the key or the specific 
               implementation used cannot perform database lookups based on cache misses
               you will not be able to retrieve the object in question.
  Returntype : ArrayRef of found values
  Exceptions : Thrown if SQL and parameters are empty and not the expected types. All
               other exceptions come from DBI/SQL operations.
  Caller     : BaseAdaptors
  Status     : Beta

=cut

sub get_by_sql {
  my ($self, $sql, $params) = @_;
  throw "No SQL given" unless $sql;
  assert_ref($params, 'ARRAY', 'params');
  my $helper = $self->adaptor()->dbc()->sql_helper();
  my $ids = $helper->execute_simple(-SQL => $sql, -PARAMS => $params);
  return $self->get_by_list($ids);
}

=head2 put

  Arg [1]    : String key to store
  Arg [2]    : Scalar value to store
  Example    : $cache->put(2, 'b');
  Description: Stores a value in the cache. Returns the previous value held 
               under that key if one existed
  Returntype : Scalar of the previously held value if one existed
  Exceptions : If key was undefined
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
  #Clear the cache
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


=head2 cached_values

  Example    : my @values = $cache->cached_values();
  Description: Returns all values held by the cache
  Returntype : List of all available values
  Exceptions : None
  Caller     : BaseAdaptors
  Status     : Beta

=cut

sub cached_values {
  my ($self) = @_;
  return values %{$self->cache()};
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
