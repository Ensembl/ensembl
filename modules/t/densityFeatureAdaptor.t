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

BEGIN { $| = 1;
	use Test;
	plan tests => 18;
}



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
ok($density_features[0]->is_stored($db));

#
# get back from database and check
#

my @stored_features = @{$dfa->fetch_all_by_Slice($slice,'GeneDensityTest', 10)};

for (my $i=0; $i< scalar(@density_features); $i++){
  ok($density_features[0]->density_value() == $stored_features[0]->density_value());
}



$multi->restore('core', 'analysis');
$multi->restore('core', 'density_type');
$multi->restore('core', 'density_feature');

