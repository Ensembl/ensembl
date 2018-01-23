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

use strict;
use warnings;

use Test::More;
use Test::Warnings;
use Data::Dumper;
$|=1;

BEGIN {	
    use_ok('Bio::EnsEMBL::Utils::Iterator');
}

my ($it, $num, $i);

# check that a coderef works as an argument to the constructor

$it = Bio::EnsEMBL::Utils::Iterator->new(sub {return 3});

is($it->next, 3, "got expected value from iterator created from coderef");

# from now on we create an iterator from an arrayref because it's simpler

$it = Bio::EnsEMBL::Utils::Iterator->new([1,2,3]);

ok($it->has_next, 'iterator still has elements');

$i = 0;

while ($num = $it->next) {
    is($num, ++$i, "$num is next element");
}

is($i, 3, 'got right number of elements');

# test empty iterator

my $empty_it = Bio::EnsEMBL::Utils::Iterator->new;

ok((not $empty_it->has_next), 'empty iterator has no elements');

# test grepping and mapping

$it = Bio::EnsEMBL::Utils::Iterator->new([1,2,3,4]);

my $even_it = $it->grep(sub {$_ % 2 == 0});

$i = 0;

while ($num = $even_it->next) {
    $i++;
    ok($num % 2 == 0, "$num is even");
}

is($i, 2, 'got right number of grepped elements back');

$it = Bio::EnsEMBL::Utils::Iterator->new([1,2,3]);

my $doubled_it = $it->map(sub {$_ * 2});

$i = 0;

while ($num = $doubled_it->next) {
    is($num, ++$i * 2, "$num = $i * 2");
}

is($i, 3, 'got right number of mapped elements back');

$it = Bio::EnsEMBL::Utils::Iterator->new([1,2,3]);

my $filtered_mapped_it = $it->map(sub {$_ * 3})->grep(sub {$_ % 2 == 0});

is($filtered_mapped_it->next, 6, "filtering and mapping together work");

ok(!$filtered_mapped_it->has_next, "only got one result from filtering and mapping");

# test combining iterators

$it = Bio::EnsEMBL::Utils::Iterator->new([1,2,3]);

$empty_it = Bio::EnsEMBL::Utils::Iterator->new;

my $it2 = Bio::EnsEMBL::Utils::Iterator->new([4,5,6]);

my $combined_it = $it->append($empty_it, $it2);

$i = 0;

while ($num = $combined_it->next) {
    $i++;
}

is($i, 6, 'got right number of elements from combined iterators');

$it = Bio::EnsEMBL::Utils::Iterator->new([1,2,3]);

my $el = $it->peek;

is($el, $it->next, 'peek did not remove element');

# test converting iterator to an arrayref

$it = Bio::EnsEMBL::Utils::Iterator->new([1,2,3]);

is_deeply($it->to_arrayref, [1,2,3], 'to_arrayref returned expected array');

# test take and skip

$it = Bio::EnsEMBL::Utils::Iterator->new([1,2,3,4])->take(2);

is_deeply($it->to_arrayref, [1,2], 'took correct elements');

$it = Bio::EnsEMBL::Utils::Iterator->new([1,2,3,4])->skip(2);

is_deeply($it->to_arrayref, [3,4], 'skipped expected elements');

# test reduce

$num = Bio::EnsEMBL::Utils::Iterator->new([1,2,3,4])->reduce(sub {$_[0] + $_[1]});

is($num, 10, "reduce calculated sum correctly");

# test reduce with initial accumulator value

$num = Bio::EnsEMBL::Utils::Iterator->new([1,2,3,4])->reduce(sub {$_[0] + $_[1]}, 10);

is($num, 20, "reduce with initial value calculated sum correctly");

$num = 0;
Bio::EnsEMBL::Utils::Iterator->new([1..12])->each(sub { $num++; });
is($num, 12, 'each iterates over all elements');



# test BaseFeatureAdaptor Iterator methods

use Bio::EnsEMBL::Test::MultiTestDB;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $db = $multi->get_DBAdaptor( "core" );
#debug( "Test database instatiated" );
ok( $db , 'Test database instantiated');

my $sa        = $db->get_SliceAdaptor;
my $ga        = $db->get_GeneAdaptor;



#Do this via direct SQL as unlikely to change
my $sql = 'SELECT sr.name, g.seq_region_start from gene g, seq_region sr where g.seq_region_id = sr.seq_region_id and sr.name = "MT_NC_001807"';
my ($sr_name, $sr_start) = @{$db->dbc->db_handle->selectrow_arrayref($sql)};
my $slice     = $sa->fetch_by_region('chromosome', $sr_name, 1, ($sr_start + 1000000));
my @genes     = @{$ga->fetch_all_by_Slice($slice)};
my $num_genes = scalar(@genes);
ok($num_genes, "Found $num_genes Gene(s) on test Slice");


SKIP:{
  #test fetch_Iterator_by_Slice and implicitly fetch_Iterator_by_Slice_method
  my $gi;

  if(! $num_genes){
	skip 'Skipping fetch_Iterator_by_Slice test due to lack of genes on slice, '.
	  'please redefined slice in test', 1;
  }
  else{
	#define chunk size such that the first chunk has no genes 
	#and the first gene crosses the 2nd and 3rd chunk;
	my $chunk_size = ($genes[0]->seq_region_start + 10)/2;
	$gi = $ga->fetch_Iterator_by_Slice($slice, undef, $chunk_size);
  }  

  my $got_gi = isa_ok($gi, 'Bio::EnsEMBL::Utils::Iterator', 'Gene Iterator');
	
 SKIP:{
	my $gene_cnt = 0;

	if(! $got_gi){
	  skip 'Skipping fetch_Iterator_by_Slice test as failed to fetch Gene Iterator', 1;
	}
	else{
	  my $gene;
	  
	  while($gene = $gi->next){
	  print $gene->stable_id."\n";
	  		$gene_cnt++;		
	  }			
	}

	ok(($gene_cnt == $num_genes), "fetch_Iterator_by_Slice returned correct number of features($gene_cnt == $num_genes)");	
  }
}

done_testing();

