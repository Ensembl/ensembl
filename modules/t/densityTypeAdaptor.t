
use strict;

BEGIN { $| = 1;
	use Test;
	plan tests => 20;
}

use Bio::EnsEMBL::Test::MultiTestDB;

use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
ok( $multi );

my $db = $multi->get_DBAdaptor( "core" );

my $dta = $db->get_DensityTypeAdaptor();

ok(ref($dta));

#
# test fetch_all_by_logic_name
#
my @dts = @{$dta->fetch_all_by_logic_name('RepeatCoverage')};

ok(@dts == 1);

ok($dts[0]->dbID == 2);
ok($dts[0]->analysis->logic_name eq 'RepeatCoverage');
ok($dts[0]->block_size == 100);
ok($dts[0]->value_type eq 'ratio');


#
# test fetch_by_dbID
#

my $dt = $dta->fetch_by_dbID(1);
ok($dt->dbID == 1);
ok($dt->analysis->logic_name eq 'SNPDensity');
ok($dt->block_size == 100);
ok($dt->value_type eq 'sum');

#
# test fetch_all
#
@dts = @{$dta->fetch_all()};
ok(@dts == 2);

#
# test store
#

$multi->save('core', 'density_type', 'analysis');

my $analysis = Bio::EnsEMBL::Analysis->new
  (-program     => "test",
   -database    => "ensembl",
   -gff_source  => "densityFeature.t",
   -gff_feature => "density",
   -logic_name  => "GeneDensityTest");


#
# test constructor
#
$dt = Bio::EnsEMBL::DensityType->new
  (-analysis   => $analysis,
   -block_size => 600,
   -value_type => 'sum');


$dta->store($dt);

ok($dt->adaptor == $dta);
ok($dt->dbID);

$dt = $dta->fetch_by_dbID($dt->dbID);
ok($dt->analysis->logic_name eq 'GeneDensityTest');
ok($dt->block_size == 600);
ok($dt->value_type eq 'sum');

# make sure analysis was stored too
ok($dt->analysis->dbID && $dt->analysis->adaptor);

my $dbID = $dt->dbID();
my $rows = count_rows($db, 'density_type');

# try to store the same density type a second time
# should not be entered in the db twice
$dt = Bio::EnsEMBL::DensityType->new
  (-analysis => $analysis,
   -block_size => 600,
   -value_type => 'sum');


$dta->store($dt);

ok($dt->dbID == $dbID);
ok(count_rows($db, 'density_type') == $rows);


$multi->restore('core', 'density_type', 'analysis');


