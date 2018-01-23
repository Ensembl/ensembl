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

use Test::More;
use Test::Warnings;
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

done_testing();
