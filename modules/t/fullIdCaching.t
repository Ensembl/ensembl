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

package fullIdCaching;

use base qw/Bio::EnsEMBL::DBSQL::Support::FullIdCache/;
use strict;
use warnings;

sub support_additional_lookups {
  return 1;
}

sub compute_keys {
  my ($self, $object) = @_;
  return { 
    biotype => $object->biotype(), 
    logic_name => $object->analysis()->logic_name(), 
    dbID => $object->dbID()
  };
}

1;

#######################

package main;

use strict;
use warnings;
use Test::More;
use Test::Warnings;
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
  no strict 'refs'; ## no critic
  *Bio::EnsEMBL::DBSQL::GeneAdaptor::_build_id_cache = sub {
    my ($self) = @_;
    return fullIdCaching->new($self);
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

  ############ Test additional lookup code
  {
    my $protein_coding_genes = $cache->get_all_by_additional_lookup('biotype', 'protein_coding');
    is_deeply([sort map { $_->dbID() } @{$protein_coding_genes} ], [sort map { $_->dbID() } @{$genes} ], 'Protein coding lookup returns all genes');

    my $wibble_genes = $cache->get_all_by_additional_lookup('biotype', 'wibble');
    is_deeply([], $wibble_genes, 'biotype is a valid lookup but wibble is a bad key. Empty array return');

    my $bad_lookup = $cache->get_all_by_additional_lookup('wibble', 'wibble');
    is_deeply([], $wibble_genes, 'wibble is a bad lookup. Empty array return');

    my $individual_gene = $cache->get_by_additional_lookup('dbID', $gene_ids->[0]);
    is($individual_gene->dbID, $gene_ids->[0], 'Lookup of dbID returns a single value');
    is($cache->get($gene_ids->[0]), $individual_gene, 'Same object is returned from the main get() method and the lookup');

    dies_ok { $cache->get_by_additional_lookup('biotype', 'protein_coding') } 'Expect to die as the query will return more than one value';

    # Add the same gene to the cache and check it's not there twice
    $cache->put($individual_gene->dbID(), $individual_gene);
    my $same_genes = $cache->get_all_by_additional_lookup('dbID', $gene_ids->[0]);
    is(scalar(@$same_genes), 1, 'The gene is still in 1 copy in the lookup');

    #Clear the cache and make sure that we can still retrieve by additional values
    {
      $cache->clear_cache();
      my $db_id = $protein_coding_genes->[0]->dbID();
      my $fetched_gene = $cache->get_by_additional_lookup('dbID', $db_id);
      ok(defined $fetched_gene, 'Checking we got a gene back if we fetched using additional lookups alone');
      cmp_ok($fetched_gene->dbID(), '==', $db_id, 'Checking we got a gene back if we fetched using additional lookups alone');
    }


    $cache->remove($individual_gene->dbID());
    my $new_protein_coding_genes = $cache->get_all_by_additional_lookup('biotype', 'protein_coding');

    ok(! defined $cache->get_by_additional_lookup('dbID', $gene_ids->[0]), 'Removed gene so lookups can no longer return an object');
    ok(! exists $cache->_additional_lookup()->{dbID}->{$gene_ids->[0]}, 'Removed the resulting array from the lookup hash');
    ok(exists $cache->_additional_lookup()->{biotype}->{protein_coding}, 'Biotype lookup still exists for protein_coding');
    is(scalar(@{$new_protein_coding_genes}), (scalar(@{$protein_coding_genes}) -1), 'Reduced the returned number of protein coding genes by one');

    $cache->put($individual_gene->dbID(), $individual_gene);
    ok(defined $cache->get_by_additional_lookup('dbID', $gene_ids->[0]), 'Added the gene back in. Everything is OK again');

    #Checking DBSQL based lookup works
    my $sql_protein_coding_genes = $cache->get_by_sql('select gene_id from gene where biotype =?', ['protein_coding']);
    is_deeply([sort map { $_->dbID() } @{$sql_protein_coding_genes} ], [sort map { $_->dbID() } @{$protein_coding_genes} ], 'SQL based protein_coding lookup returns all genes');
  }
  
  #Turn off caching; we should get a fresh object out of the cache
  my $cached_obj = $adaptor->fetch_by_dbID($gene_ids->[0]);
  is($adaptor->fetch_by_dbID($gene_ids->[0]), $cached_obj, 'Checking two objects are the same i.e. cached');
  is($adaptor->fetch_all_by_dbID_list([$gene_ids->[0]])->[0], $cached_obj, 'Checking two objects are the same i.e. cached');
  $adaptor->db()->no_cache(1);
  isnt($adaptor->fetch_by_dbID($gene_ids->[0]), $cached_obj, 'Checking two objects are no longer the same as the cache is off');
  isnt($adaptor->fetch_all_by_dbID_list([$gene_ids->[0]])->[0], $cached_obj, 'Checking two objects are no longer the same as the cache is off');
  $adaptor->db()->no_cache(0);
  
  # Negating no_cache
  $cached_obj = $adaptor->fetch_by_dbID($gene_ids->[0]);
  $adaptor->db()->no_cache(1);
  $adaptor->ignore_cache_override(1);
  is($adaptor->fetch_by_dbID($gene_ids->[0]), $cached_obj, 'Checking no_cache can be ignored.');
  $adaptor->db()->no_cache(0);
  $adaptor->ignore_cache_override(0);
  
  
  #hit seq region cache & check we clear that one out as well
  $adaptor->fetch_all_by_Slice($genes->[0]->slice());
  $adaptor->clear_cache();
  is(scalar(keys(%{$gene_adaptor->_slice_feature_cache()})), 0, 'clear_cache() should have hit the slice feature cache as well');
  ok(! defined $adaptor->{_id_cache}->{cache}, 'Cache clear has deleted the hash');
  
  #Quick repopulate of the cache then store
  if($STORE_TESTS) {
    $multi->save('core', 'gene', 'transcript');
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
    
    $multi->restore('core');
  }
}

done_testing();
