

use strict;
use warnings;

use Bio::EnsEMBL::DensityFeature;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Analysis;

use Bio::EnsEMBL::Test::TestUtils;


our $verbose = 0; #set to 1 to turn on debug printouts



BEGIN { $| = 1;
	use Test;
	plan tests => 4;
}

use Bio::EnsEMBL::Test::TestUtils;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db = $multi->get_DBAdaptor('core');


my $analysis =  new Bio::EnsEMBL::Analysis (-program     => "densityFeature.t",
					   -database    => "ensembl",
					   -gff_source  => "densityFeature.t",
					   -gff_feature => "density",
					   -logic_name  => "GeneDensityTest");

my $dt = Bio::EnsEMBL::DensityType->new(-analysis   => $analysis,
					-block_size => 600,
					-value_type => 'sum');

my $slice_adaptor = $db->get_SliceAdaptor();
my $slice = $slice_adaptor->fetch_by_region('chromosome', '20', 1, 600);


#
#test the constructor
#
my $feat = Bio::EnsEMBL::DensityFeature->new(-seq_region    => $slice,
				             -start         => 1,
					     -end           => 300,
					     -density_type  => $dt,
					     -density_value => 123);

ok($feat && ref $feat && $feat->isa('Bio::EnsEMBL::DensityFeature'));


#
# Test the getter setter functions;
#

ok(&test_getter_setter($feat, 'start', 100));
ok(&test_getter_setter($feat, 'end', 500));
ok(&test_getter_setter($feat, 'density_value', 456));


