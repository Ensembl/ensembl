# Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

use Config;
use Test::More;
use File::Temp qw/tempfile/;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils qw/warns_like/;

my $threads;
if($Config{useithreads} && ! $ENV{ENS_FORCE_NOTHREADS}) {
  note 'Using threaded tests';
  require threads;
  $threads = 1;
}
else {
  note 'Using non-threaded tests';
}

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $multi_db->get_DBAdaptor('core');
my $dbc = $db->dbc();

my $reg = 'Bio::EnsEMBL::Registry';

my $registry_template = <<'TMPL';
{
  package Reg;
  use Bio::EnsEMBL::DBSQL::DBAdaptor;
  Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -HOST => '%s',
    -PORT => %d,
    -USER => '%s',
    -PASSWORD => '%s',
    -DBNAME => '%s',
    -DRIVER => '%s',
    -SPECIES => 'new'
  );
}
1;
TMPL

#Testing threaded re-loads of the registry
{
  my ($fh, $filename) = tempfile();
  my $final = sprintf($registry_template,
                      $dbc->host() || '',
                      $dbc->port() || 0,
                      $dbc->username() || '',
                      $dbc->password() || '',
                      $dbc->dbname(),
                      $dbc->driver(),
      );
  print $fh $final;
  close $fh;
  
  my $call = sub {
    my @results;
    foreach my $inc (0..9) {
      push(@results, $reg->load_all($filename));
    }
    return \@results;
  };
  
  my @results;
  my @expected;
  if($threads) {
    $ENV{RUNTESTS_HARNESS_NORESTORE} = 1;
    my @thrds;
    foreach my $thr (0..9) {
      push(@thrds, threads->create($call));
    }
    foreach my $thr (@thrds) {
      my $results = $thr->join();
      my $msg = sprintf('THREAD %s: Checking first call loaded 1 DBAdaptor and the remaining 9 did nothing', $thr->tid());
      is_deeply($results, [1, ((0)x9)], $msg);
    }
  }
  else {
    foreach my $itr (1..10) {
      my $res = $call->();
      is_deeply($res, [1, (0)x9], "Testing iteration 1 where we can load an adaptor") if $itr == 1;
      is_deeply($res, [(0)x10], "Testing iteration $itr where we can no longer load adaptors") if $itr > 1;
    }
    ok("Calling of single-threaded load went off without any problems");
  }
}

#Testing auto-correction of arguments for common 1st line methods
{
  my $tester = sub {
    my ($misspelling) = @_;
    my %params = (-HOST => $dbc->host(), -PORT => $dbc->port(), -USER => $dbc->username());
    $params{-PASS} = $dbc->password() if $dbc->password();
    my $db_version = -2;
    $params{"-${misspelling}"} = $db_version;
    warns_like( sub { $reg->load_registry_from_db(%params) }, qr/${misspelling}.+mis-spelling/, "Testing that param -${misspelling} succeeded");
    return;
  };
  $tester->('dbversion');
  $tester->('version');
  $tester->('verion');
  $tester->('verison');
}

# Test get_all_species

my @species = $reg->get_all_species();
ok(scalar(@species) == 1, "get_all_species");
ok(scalar(@{ $reg->get_all_species('cahoona') }) == 0, "get_all_species with bogus data.");

done_testing();
