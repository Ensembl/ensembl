# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

## no critic (RequireFilenameMatchesPackage)

package main;

use strict;
use warnings;
use Test::More;
use Test::Warnings;
use Test::Exception;
use Bio::EnsEMBL::DBSQL::GeneAdaptor;
use Bio::EnsEMBL::DBSQL::Support::LruIdCache;
use Bio::EnsEMBL::Test::MultiTestDB;

my $STORE_TESTS = 1;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $multi->get_DBAdaptor('core');
my $gene_adaptor = $dba->get_GeneAdaptor();

my $gene_ids = [18256, 18257, 18258];
my $genes = $gene_adaptor->fetch_all_by_dbID_list($gene_ids);

sub BEGIN {
  no strict 'refs'; ##no critic
  *Bio::EnsEMBL::DBSQL::GeneAdaptor::_build_id_cache = sub {
    my ($self) = @_;
    return Bio::EnsEMBL::DBSQL::Support::LruIdCache->new($self, 3);
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
  
  my $cached_ids = $adaptor->fetch_all_by_dbID_list($gene_ids);
  my $refetched_cached_ids = $adaptor->fetch_all_by_dbID_list($gene_ids);
  
  is_deeply($refetched_cached_ids, $cached_ids, 'Ensuring two calls which do not exceed the cache returns the same values');
  my $first_gene_refetch = $adaptor->fetch_by_dbID($gene_ids->[0]);
  is_deeply($first_gene_refetch, $refetched_cached_ids->[0], 'fetch_by_dbID() should return the same object as fetch_all_by_dbID_list()');
  is_deeply([sort $cache->cache_keys()], $gene_ids, 'Cache still holds original values');
  ok(! defined $adaptor->fetch_by_dbID(1), 'Fetching with a bad ID returns nothing');
  is_deeply([sort $cache->cache_keys()], $gene_ids, 'Cache still holds original values');
  ok(defined $adaptor->fetch_by_dbID(18259), 'New fetch still returns values');
  is(scalar($cache->cache_keys()), 3, 'Cache can only hold 3 elements');
  cmp_ok([sort $cache->cache_keys()]->[2], '!=', $gene_ids->[2], 'Last element is not the same as our original set');
    
  my $final_refetch = $adaptor->fetch_all_by_dbID_list($gene_ids);
  cmp_ok([sort $cache->cache_keys()]->[2], '==', $gene_ids->[2], 'Refetching brings the expected gene back into the cache');
  
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
  is(scalar($cache->cache_keys()), 0, 'Cache clear has resulted in no more elements being held');

  #### Quick test of the get_by_sql method
  {
    my $genes_refetch = $adaptor->fetch_all_by_dbID_list($gene_ids);
    is($cache->size(), 3, 'We should have only persisted 3 values in the cache');
    my $original_keys = join(q{!=!}, $cache->cache_keys());
    my $sql = 'select gene_id from gene where biotype =?';
    my $pc_sql_genes = $cache->get_by_sql($sql, ['protein_coding']);
    isnt(join(q{!=!}, $cache->cache_keys()), $original_keys, 'Using SQL genes should have reset the IDs in the cache');
    $adaptor->clear_cache();
  }
  
  #Quick repopulate of the cache then store
  if($STORE_TESTS) {
    $adaptor->fetch_all_by_dbID_list($gene_ids);
    my $old_gene = $cached_ids->[0];
    my $new_gene = Bio::EnsEMBL::Gene->new(
      -SLICE => $old_gene->slice(), -START => $old_gene->start(), -END => $old_gene->end(), -STRAND => $old_gene->strand(), -ANALYSIS => $old_gene->analysis()
    );
    $new_gene->add_Transcript($old_gene->get_all_Transcripts()->[0]);
    my $new_id = $adaptor->store($new_gene);
    is(scalar($cache->cache_keys()), 0, 'store() has caused a clear_cache() call');
    my $id = $adaptor->fetch_by_dbID($new_id);
    is_deeply([$cache->cache_keys()], [$new_id], 'fetch on the new ID has worked');
    
    #Repopulate cache & then update a gene
    $adaptor->fetch_all_by_dbID_list($gene_ids);
    $adaptor->update($cached_ids->[0]);
    is(scalar($cache->cache_keys()), 0, 'update() has caused a clear_cache() call');
    $adaptor->fetch_by_dbID($new_id);
    is_deeply([$cache->cache_keys()], [$new_id], 'fetch on the new ID has worked');
  }
  $multi->restore('core');
}

done_testing();
