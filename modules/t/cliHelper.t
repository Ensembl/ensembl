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
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Operon;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Utils::CliHelper;

debug("Startup test");
ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $dba = $multi->get_DBAdaptor("core");

SKIP: {
skip 'CliHelper not supported for SQLite yet', 1 if $dba->dbc()->driver() eq 'SQLite';

my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();

debug("Checking default options");
my $opts = { host   => $dba->dbc()->host(),
             user   => $dba->dbc()->username(),
             pass   => $dba->dbc()->password(),
             port   => $dba->dbc()->port(),
             dbname => $dba->dbc()->dbname(), };

my $dba_args = $cli_helper->get_dba_args_for_opts($opts);

is( scalar(@$dba_args),        1 );
is( $dba_args->[0]->{-HOST},   $opts->{host} );
is( $dba_args->[0]->{-USER},   $opts->{user} );
is( $dba_args->[0]->{-PASS},   $opts->{pass} );
is( $dba_args->[0]->{-PORT},   $opts->{port} );
is( $dba_args->[0]->{-DBNAME}, $opts->{dbname} );
ok( !defined $dba_args->[0]->{-SPECIES} );
is( $dba_args->[0]->{-SPECIES_ID},1 );
is( $dba_args->[0]->{-MULTISPECIES_DB}, 0);

$opts->{species_id} = 1;
$opts->{species} = "homo_sapiens";
$dba_args = $cli_helper->get_dba_args_for_opts($opts);

is( scalar(@$dba_args),        1 );
is( $dba_args->[0]->{-HOST},   $opts->{host} );
is( $dba_args->[0]->{-USER},   $opts->{user} );
is( $dba_args->[0]->{-PASS},   $opts->{pass} );
is( $dba_args->[0]->{-PORT},   $opts->{port} );
is( $dba_args->[0]->{-DBNAME}, $opts->{dbname} );
is( $dba_args->[0]->{-SPECIES}, $opts->{species} );
is( $dba_args->[0]->{-SPECIES_ID}, $opts->{species_id} );
is( $dba_args->[0]->{-MULTISPECIES_DB}, 0 );


my $srcopts = { srchost   => $dba->dbc()->host(),
                srcuser   => $dba->dbc()->username(),
                srcpass   => $dba->dbc()->password(),
                srcport   => $dba->dbc()->port(),
                srcdbname => $dba->dbc()->dbname(), };

debug("Checking prefix options without single species specified");
my $src_dba_args =
  $cli_helper->get_dba_args_for_opts( $srcopts, 0, "src" );

is( scalar(@$src_dba_args),        1 );
is( $src_dba_args->[0]->{-HOST},   $srcopts->{srchost} );
is( $src_dba_args->[0]->{-USER},   $srcopts->{srcuser} );
is( $src_dba_args->[0]->{-PASS},   $srcopts->{srcpass} );
is( $src_dba_args->[0]->{-PORT},   $srcopts->{srcport} );
is( $src_dba_args->[0]->{-DBNAME}, $srcopts->{srcdbname} );
ok( defined $src_dba_args->[0]->{-SPECIES} );
is( $src_dba_args->[0]->{-SPECIES_ID},1 );
is( $src_dba_args->[0]->{-MULTISPECIES_DB}, 0 );

debug("Checking prefix options with single species specified");
$src_dba_args =
  $cli_helper->get_dba_args_for_opts( $srcopts, 1, "src" );

is( scalar(@$src_dba_args),        1 );
is( $src_dba_args->[0]->{-HOST},   $srcopts->{srchost} );
is( $src_dba_args->[0]->{-USER},   $srcopts->{srcuser} );
is( $src_dba_args->[0]->{-PASS},   $srcopts->{srcpass} );
is( $src_dba_args->[0]->{-PORT},   $srcopts->{srcport} );
is( $src_dba_args->[0]->{-DBNAME}, $srcopts->{srcdbname} );
is( $src_dba_args->[0]->{-SPECIES_ID},1 );
ok( !defined $src_dba_args->[0]->{-SPECIES} );
is( $src_dba_args->[0]->{-MULTISPECIES_DB}, 0 );

debug("Checking prefix options with single species undefined");
$src_dba_args =
  $cli_helper->get_dba_args_for_opts( $srcopts, undef, "src" );

is( scalar(@$src_dba_args),        1 );
is( $src_dba_args->[0]->{-HOST},   $srcopts->{srchost} );
is( $src_dba_args->[0]->{-USER},   $srcopts->{srcuser} );
is( $src_dba_args->[0]->{-PASS},   $srcopts->{srcpass} );
is( $src_dba_args->[0]->{-PORT},   $srcopts->{srcport} );
is( $src_dba_args->[0]->{-DBNAME}, $srcopts->{srcdbname} );
is( $src_dba_args->[0]->{-SPECIES_ID},1 );
ok( !defined $src_dba_args->[0]->{-SPECIES} );
is( $src_dba_args->[0]->{-MULTISPECIES_DB}, 0 );

$multi->restore('core');

debug("Test database restored");
ok($dba);

my $nameless = Bio::EnsEMBL::Test::MultiTestDB->new("nameless");
$dba = $nameless->get_DBAdaptor("core");

debug("Nameless test database instantiated");
ok($dba);
debug("Checking default options");
$opts = { host   => $dba->dbc()->host(),
             user   => $dba->dbc()->username(),
             pass   => $dba->dbc()->password(),
             port   => $dba->dbc()->port(),
             dbname => $dba->dbc()->dbname(), };

$dba_args = $cli_helper->get_dba_args_for_opts($opts);

is( scalar(@$dba_args),        1 );
is( $dba_args->[0]->{-HOST},   $opts->{host} );
is( $dba_args->[0]->{-USER},   $opts->{user} );
is( $dba_args->[0]->{-PASS},   $opts->{pass} );
is( $dba_args->[0]->{-PORT},   $opts->{port} );
is( $dba_args->[0]->{-DBNAME}, $opts->{dbname} );
ok( !defined $dba_args->[0]->{-SPECIES} );
is( $dba_args->[0]->{-SPECIES_ID},1 );
is( $dba_args->[0]->{-MULTISPECIES_DB}, 0 );

my $collection = Bio::EnsEMBL::Test::MultiTestDB->new('test_collection');

$dba = $collection->get_DBAdaptor("core");

debug("Collection test database instantiated");
ok($dba);
debug("Checking default options");
$opts = { host   => $dba->dbc()->host(),
             user   => $dba->dbc()->username(),
             pass   => $dba->dbc()->password(),
             port   => $dba->dbc()->port(),
             dbname => $dba->dbc()->dbname(), };

$dba_args = $cli_helper->get_dba_args_for_opts($opts);

is(scalar(@$dba_args),2);
is( $dba_args->[0]->{-HOST},   $opts->{host} );
is( $dba_args->[0]->{-USER},   $opts->{user} );
is( $dba_args->[0]->{-PASS},   $opts->{pass} );
is( $dba_args->[0]->{-PORT},   $opts->{port} );
is( $dba_args->[0]->{-DBNAME}, $opts->{dbname} );
ok( defined $dba_args->[0]->{-SPECIES} );
is( $dba_args->[0]->{-SPECIES_ID},1 );
is( $dba_args->[0]->{-MULTISPECIES_DB}, 1 );
is( $dba_args->[1]->{-HOST},   $opts->{host} );
is( $dba_args->[1]->{-USER},   $opts->{user} );
is( $dba_args->[1]->{-PASS},   $opts->{pass} );
is( $dba_args->[1]->{-PORT},   $opts->{port} );
is( $dba_args->[1]->{-DBNAME}, $opts->{dbname} );
ok( defined $dba_args->[1]->{-SPECIES} );
is( $dba_args->[1]->{-SPECIES_ID},2 );
is( $dba_args->[1]->{-MULTISPECIES_DB}, 1 );

} # SKIP for SQLite

done_testing();
