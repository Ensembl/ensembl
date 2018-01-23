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

Bio::EnsEMBL::DBSQL::Support::LruIdCache - ID based caching using an LRU backed cache  

=head1 SYNOPSIS

  my $cache = Bio::EnsEMBL::DBSQL::Support::LruIdCache->new($adaptor, 2);
  my $obj = $cache->put(1, 'a');
  my $obj = $cache->put(2, 'b');
  my $obj = $cache->put(3, 'c');
  
  is_deeply([$cache->cache_keys()], ['b','c']); #ID 1 was ejected under LRU rules

=head1 DESCRIPTION

An implementation of caching which uses the oldest accessed key as the 
value to be ejected from the cache when the maximum size has been hit. See
the following page for more information about this algorithm

http://en.wikipedia.org/wiki/Cache_algorithms#Least_Recently_Used

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::Support::LruIdCache;

use strict;
use warnings;
use base qw/Bio::EnsEMBL::DBSQL::Support::BaseCache/;
use Bio::EnsEMBL::Utils::Cache;
use Bio::EnsEMBL::Utils::Scalar qw/assert_ref/;

my $DEFAULT = 1000;

=head2 new

  Arg [1]    : Bio::EnsEMBL::DBSQL::BaseAdaptor $db_adaptor
  Arg [2]    : Int size of the cache. Defaults to 1000
  Example    : my $cache = Bio::EnsEMBL::DBSQL::Support::LruIdCache->new($db_adaptor, 10);
  Description: Creates a new cache class instance
  Returntype : Bio::EnsEMBL::DBSQL::Support::BaseCache
  Exceptions : none
  Caller     : BaseAdaptors
  Status     : Beta

=cut

sub new {
  my ($class, $adaptor, $lru_size) = @_;
  my $self = $class->SUPER::new($adaptor);
  $lru_size ||= $DEFAULT;
  $self->lru_size($lru_size);
  return $self;
}

=head2 lru_size

  Arg [1]    : Int size of the cache
  Example    : $cache->lru_size(10);
  Description: Resets the size of the cache and forces a rebuild of the cache 
               object to apply this new sizing. Also functions as a getter
  Returntype : Int
  Exceptions : none
  Caller     : BaseAdaptors
  Status     : Beta

=cut

sub lru_size {
  my ($self, $lru_size) = @_;
  if(defined $lru_size) {
    $self->{'lru_size'} = $lru_size;
    $self->delete_cache();
  }
  return $self->{'lru_size'};
}

=head2 build_cache

  Description: Returns an instance of C<Bio::EnsEMBL::Utils::Cache> and sized
               according to the return value in C<lru_cache()>
  Returntype : Bio::EnsEMBL::Utils::Cache
  Exceptions : none
  Caller     : BaseAdaptors
  Status     : Beta

=cut

sub build_cache {
  my ($self) = @_;
  tie my %cache, 'Bio::EnsEMBL::Utils::Cache', $self->lru_size();
  return \%cache;
}

=head2 get

  Arg [1]    : String key to retrieve
  Example    : $is($cache->get(1), 'a');
  Description: Retrieves the value held in the cache. If the value is not in
               the cache we will retrieve the value from 
               C<_uncached_fetch_by_dbID> and then store that value
  Returntype : Scalar value held in the cache or nothing if the ID was not present
  Exceptions : If key was undefined
  Caller     : BaseAdaptors
  Status     : Beta

=cut

sub get {
  my ($self, $key) = @_;
  my $value = $self->SUPER::get($key);
  if(! defined $value) {
    my $new_value = $self->adaptor()->_uncached_fetch_by_dbID($key);
    if(defined $new_value) {
      $self->put($key, $new_value);
      $value = $new_value;
    }
  }
  return $value if defined $value;
  return;
}

=head2 get_by_list

  Arg [1]    : ArrayRef keys to retrieve
  Arg [2]    : Bio::EnsEMBL::Slice optional attribute for 
               C<_uncached_fetch_all_by_dbID_list()> delegation
  Example    : is($cache->get_by_list([1,2]), ['a','b']);
  Description: Attempts to retrieve all values currently available in the cache,
               fetches any remaining values and stores these in the cache.
  Returntype : ArrayRef of found values
  Exceptions : None
  Caller     : BaseAdaptors
  Status     : Beta

=cut

sub get_by_list {
  my ($self, $list, $slice) = @_;
  assert_ref($list, 'ARRAY', 'list');
  my $results = $self->SUPER::get_by_list($list);
  my %available = map { $_->dbID(), $_ } @{$results};
  my @to_fetch;
  foreach my $input (@{$list}) {
    next if ! defined $input;
    push(@to_fetch, $input) if ! exists $available{$input};
  }
  
  my $fetched = $self->adaptor()->_uncached_fetch_all_by_dbID_list(\@to_fetch, $slice);
  foreach my $obj (@{$fetched}) {
    push(@{$results}, $obj);
    $self->put($obj->dbID(), $obj);
  }
  
  return $results;
}

1;
