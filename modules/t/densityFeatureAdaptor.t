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

use Bio::EnsEMBL::Utils::Exception qw(verbose);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::DensityTypeAdaptor;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::DensityFeature;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

use Test::More;
use Test::Warnings;
our $verbose = 0;
verbose('WARNING');

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db = $multi->get_DBAdaptor('core');


my $dfa = $db->get_DensityFeatureAdaptor();

ok(ref($dfa) && $dfa->isa('Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor'));

#
# Teststore
#

$multi->save('core', 'analysis');
$multi->save('core', 'density_type');
$multi->save('core', 'density_feature');



my $aa  = $db->get_AnalysisAdaptor();
my $analysis = new Bio::EnsEMBL::Analysis (-program     => "densityFeatureAdaptor.t",
					   -database    => "ensembl",
					   -gff_source  => "densityFeatureAdaptor.t",
					   -gff_feature => "density",
					   -logic_name  => "GeneDensityTest");
ok(!$analysis->is_stored($db));
$aa->store($analysis);
ok($analysis->is_stored($db));


my $block_size = 6000000;

my $dta = $db->get_DensityTypeAdaptor();
my $dt = Bio::EnsEMBL::DensityType->new(-analysis   => $analysis,
					-block_size => $block_size,
					-value_type => 'sum');

ok(!$dt->is_stored($db));
$dta->store($dt);	
ok($dt->is_stored($db));


my $slice_adaptor = $db->get_SliceAdaptor();
my $slice = $slice_adaptor->fetch_by_region('chromosome', '20', 1, ($block_size*10));

my $start = $slice->start();
my $end = ($start + $block_size)-1;
my $term = $slice->start+$slice->length;
  
my @density_features=();
while($start < $term){
  my $sub_slice = $slice_adaptor->fetch_by_region('chromosome','20',$start,$end);
  my $count =0;
  foreach my $gene (@{$sub_slice->get_all_Genes()}){
    if($gene->analysis()->logic_name() ne "pseudogene" and $gene->start >=1 ){
      $count++
    }
  }
	
#  print $count."\n";

  push @density_features, Bio::EnsEMBL::DensityFeature->new(-seq_region    => $slice,
							     -start         => $start,
							     -end           => $end,
							     -density_type  => $dt,
							     -density_value => $count);

  $start = $end+1;
  $end   = ($start + $block_size)-1;
}

ok(scalar( @density_features) == 10); 

ok(!$density_features[0]->is_stored($db));
$dfa->store(@density_features);		

#
# get back from database and check
#

my @filtered_features = grep { $_->density_value() != 0 } @density_features;

my @stored_features = @{$dfa->fetch_all_by_Slice($slice,'GeneDensityTest', 10)};

for (my $i=0; $i< scalar(@filtered_features); $i++){
  ok($filtered_features[0]->density_value() == $stored_features[0]->density_value());
}



#
# Now some density features with region_feature count set on density type
# Lets say we want 300 features on our seq_regions
# 


debug( "Region Features enabled densities" );

$dt = Bio::EnsEMBL::DensityType->new(-analysis   => $analysis,
				     -region_features => 300,
				     -value_type => 'sum');

ok(!$dt->is_stored($db));
$dta->store($dt);	
ok($dt->is_stored($db));

@density_features = ();
my $chr_20_slice = $slice_adaptor->fetch_by_region( "chromosome", 20 );
my $step = POSIX::ceil( $chr_20_slice->length() / 300);
$start = 1;
while( $start < $chr_20_slice->length() ) {
  my $end = $start+$step -1;
  if( $end > $chr_20_slice->length ) { $end = $chr_20_slice->length();}
  push @density_features, Bio::EnsEMBL::DensityFeature->new(-seq_region    => $chr_20_slice,
							     -start         => $start,
							     -end           => $end,
							     -density_type  => $dt,
							     -density_value => 5 );
  $start += $step;
}

$dfa->store( @density_features );
ok($density_features[-1]->is_stored($db));
debug( "Created ".scalar( @density_features )." density features on chr 20" );

@stored_features = @{$dfa->fetch_all_by_Slice($chr_20_slice,'GeneDensityTest', 100, "interpolate" )};
ok( scalar( @stored_features ) == 100 );
debug( "Interpolated retrieved ".scalar( @stored_features ));

ok( abs( $stored_features[0]->density_value() - 15) < 0.0001 );
debug( "Density value = ".$stored_features[0]->density_value() );


# test the retreival of the right sized features
# first yet another density size
# comes to about 1000bp on chr 20

$dt = Bio::EnsEMBL::DensityType->new(-analysis   => $analysis,
				     -region_features => 30_000,
			             -value_type => 'sum');
$dta->store($dt);	

# need some features 
@density_features = ();
$step = POSIX::ceil( $chr_20_slice->length() / 30_000);
for my $arr ( [1,1000,1],[1001,2000,2],[2001,3000,3],[31_000_000,31_000_999,31.0],[31_500_000, 31_500_999,31.5] ) {
  push @density_features, Bio::EnsEMBL::DensityFeature->new(-seq_region    => $chr_20_slice,
							     -start         => $arr->[0],
							     -end           => $arr->[1],
							     -density_type  => $dt,
							     -density_value => $arr->[2] );
}
$dfa->store( @density_features );

# now check for retrieval
my $sub_Slice = $chr_20_slice->sub_Slice( 1, 1000_000 );
@stored_features = @{$dfa->fetch_all_by_Slice( $sub_Slice, 'GeneDensityTest', 2 )};
ok( $stored_features[0]->length() > 10000 );
@stored_features = @{$dfa->fetch_all_by_Slice( $sub_Slice, 'GeneDensityTest', 10 )};
ok( $stored_features[0]->length() > 10000 );
@stored_features = @{$dfa->fetch_all_by_Slice( $sub_Slice, 'GeneDensityTest', 100 )};
print $stored_features[0]->length() . "\n";
ok( $stored_features[0]->length() == 1000 );

# Check for fetching all features
@stored_features = @{$dfa->fetch_all()};
print "\n";
is( @stored_features, 318, 'Number of stored features');
@stored_features = @{$dfa->fetch_all('GeneDensityTest')};
ok( $stored_features[0]->length() > 1000 );
is( @stored_features, 306, "Number of stored gene densities"); 
@stored_features = @{$dfa->fetch_all('rubbish')};
ok( @stored_features == 0);





$multi->restore('core', 'analysis');
$multi->restore('core', 'density_type');
$multi->restore('core', 'density_feature');

done_testing();
