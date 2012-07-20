package main;

use strict;
use warnings;
use Test::More;
use Test::Exception;
use Bio::EnsEMBL::DBSQL::GeneAdaptor;
use Bio::EnsEMBL::DBSQL::Support::FullIdCache;
use Bio::EnsEMBL::Test::MultiTestDB;

my $STORE_TESTS = 1;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $multi->get_DBAdaptor('core');
my $gene_adaptor = $dba->get_GeneAdaptor();

my $gene_ids = [18256, 18257, 18258];
my $genes = $gene_adaptor->fetch_all_by_dbID_list($gene_ids);

sub BEGIN {
  no strict 'refs';
  *Bio::EnsEMBL::DBSQL::GeneAdaptor::_build_id_cache = sub {
    my ($self) = @_;
    return Bio::EnsEMBL::DBSQL::Support::FullIdCache->new($self);
  };
  no warnings 'redefine';
  my $original_store = \&Bio::EnsEMBL::DBSQL::GeneAdaptor::store; 
  *Bio::EnsEMBL::DBSQL::GeneAdaptor::store = sub {
    my ($self, @args) = @_;
    my $id = $original_store->($self, @args);
    $self->clear_cache();
    return $id if $id;
    return;
  };
  my $original_update = \&Bio::EnsEMBL::DBSQL::GeneAdaptor::update; 
  *Bio::EnsEMBL::DBSQL::GeneAdaptor::update = sub {
    my ($self, @args) = @_;
    $original_update->($self, @args);
    $self->clear_cache();
    return;
  };
}

{
  $multi->save('core', 'gene', 'transcript');
  
  my $adaptor = $gene_adaptor;
  my $cache = $adaptor->_id_cache();
  
  dies_ok { $cache->get(undef) } 'Getting an undef from the cache should throw an error';
  dies_ok { $cache->put(undef, undef) } 'Putting an undef into the cache should throw an error';
  
  my $genes = $gene_adaptor->fetch_all();
  
  is_deeply([sort $cache->cache_keys()], [sort map { $_->dbID() } @{$genes} ], 'Cache holds all genes');
  
  cmp_ok($cache->size(), '>', 3, 'Checking we have more than 3 values in our cache');
  
  my $cached_ids = $adaptor->fetch_all_by_dbID_list($gene_ids);
  my $refetched_cached_ids = $adaptor->fetch_all_by_dbID_list($gene_ids);
  
  is_deeply($refetched_cached_ids, $cached_ids, 'Ensuring two calls return the same values');
  my $first_gene_refetch = $adaptor->fetch_by_dbID($gene_ids->[0]);
  is_deeply($first_gene_refetch, $refetched_cached_ids->[0], 'fetch_by_dbID() should return the same object as fetch_all_by_dbID_list()');
  ok(! defined $adaptor->fetch_by_dbID(1), 'Fetching with a bad ID returns nothing');
  is(scalar(@{$adaptor->fetch_all_by_dbID_list([1])}), 0, 'Fetching with a bad ID returns an empty array');
  
  #Turn off caching; we should get a fresh object out of the cache
  my $cached_obj = $adaptor->fetch_by_dbID($gene_ids->[0]);
  is($adaptor->fetch_by_dbID($gene_ids->[0]), $cached_obj, 'Checking two objects are the same i.e. cached');
  is($adaptor->fetch_all_by_dbID_list([$gene_ids->[0]])->[0], $cached_obj, 'Checking two objects are the same i.e. cached');
  $adaptor->db()->no_cache(1);
  isnt($adaptor->fetch_by_dbID($gene_ids->[0]), $cached_obj, 'Checking two objects are no longer the same as the cache is off');
  isnt($adaptor->fetch_all_by_dbID_list([$gene_ids->[0]])->[0], $cached_obj, 'Checking two objects are no longer the same as the cache is off');
  $adaptor->db()->no_cache(0);
  
  #hit seq region cache & check we clear that one out as well
  $adaptor->fetch_all_by_Slice($genes->[0]->slice());
  $adaptor->clear_cache();
  
  is(scalar(keys(%{$gene_adaptor->_slice_feature_cache()})), 0, 'clear_cache() should have hit the slice feature cache as well');
  ok(! defined $adaptor->{_id_cache}->{cache}, 'Cache clear has deleted the hash');
  
  #Quick repopulate of the cache then store
  if($STORE_TESTS) {
    $adaptor->fetch_all_by_dbID_list($gene_ids);
    my $old_gene = $cached_ids->[0];
    my $new_gene = Bio::EnsEMBL::Gene->new(
      -SLICE => $old_gene->slice(), -START => $old_gene->start(), -END => $old_gene->end(), -STRAND => $old_gene->strand(), -ANALYSIS => $old_gene->analysis()
    );
    $new_gene->add_Transcript($old_gene->get_all_Transcripts()->[0]);
    my $new_id = $adaptor->store($new_gene);
    ok(! defined $adaptor->{_id_cache}->{cache}, 'Cache clear has deleted the hash');
    ok(scalar(grep { $_ == $new_id} $cache->cache_keys()), 'cache already has the new ID');
    
    #Repopulate cache & then update a gene
    $adaptor->fetch_all_by_dbID_list($gene_ids);
    $adaptor->update($cached_ids->[0]);
    ok(! defined $adaptor->{_id_cache}->{cache}, 'Cache clear has deleted the hash');
    ok(scalar(grep { $_ == $new_id} $cache->cache_keys()), 'cache already has the new ID');
  }
  
  $multi->restore('core');
}

done_testing();
