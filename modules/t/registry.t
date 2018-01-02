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

use Config;
use Test::More;
use Test::Deep;
use Test::Warnings;
use Test::Exception;
use File::Temp qw/tempfile/;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::ApiVersion qw( software_version );
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

my $reg = 'Bio::EnsEMBL::Registry';

# Check whilst the Registry is empty
is_deeply($reg->get_all_DBAdaptors(), [], 'get_all_DBAdaptors() returns an array-ref even when the Registry is empty');

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $multi_db->get_DBAdaptor('core');
my $dbc = $db->dbc();

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

my @species = @{ $reg->get_all_species() };
ok(scalar(@species) == 1, "get_all_species");
ok(scalar(@{ $reg->get_all_species('cahoona') }) == 0, "get_all_species with bogus data.");

# Test get_all_DBAdaptors
my $registry_register_dba = $Bio::EnsEMBL::Registry::registry_register{'_DBA'};
is( scalar(@{$reg->get_all_DBAdaptors()}), scalar(@{$registry_register_dba}), "get_all_DBAdaptors() on all species and groups" );
cmp_deeply( $reg->get_all_DBAdaptors(), shallow($registry_register_dba), "get_all_DBAdaptors() on all species and groups: comparing references" );
ok(scalar(@{$reg->get_all_DBAdaptors(-SPECIES => $species[0])}), "get_all_DBAdaptors() on a valid species");
warns_like(
    sub { is(scalar(@{$reg->get_all_DBAdaptors(-SPECIES => 'cahoona')}), 0, "get_all_DBAdaptors() on a non-existing species"); },
    qr/cahoona is not a valid species name/,
    q{Warns that the species doesn't exist},
);

dies_ok { $reg->load_all('i really hope there is no file named this way', undef, undef, undef, 1) } 'Pointing to a non-existing file should throw an error (if the option is switched on)';
is($reg->load_all('i really hope there is no file named this way'), 0, 'Pointing to a non-existing file does not throw an error if the option is switched off');

dies_ok { $reg->add_DBAdaptor() } 'add_DBAdaptor() must get a valid species';

# [ENSCORESW-2509]. Test correct handling of multi-species databases
SKIP: {
  skip 'Tests for collection DBs are done on MySQL engine', 5 unless $dbc->driver eq 'mysql';
  my @collection_dbs = create_collection_dbs();

  my %params = (-HOST => $dbc->host(), -PORT => $dbc->port(), -USER => $dbc->username());
  $params{-PASS} = $dbc->password() if $dbc->password();
  $reg->load_registry_from_db(%params);

  @species = grep { !/new/ } @{$reg->get_all_species()};
  is(scalar @species, 4, "Number of species loaded from collection DBs");
  ok(grep(/escherichia_coli_dh1/, @species), "Loaded ensembl genomes species from collection DB");
  ok(grep(/mus_musculus_aj/, @species), "Loaded ensembl species from collection DB");
  ok($reg->get_DBAdaptor("lactobacillus_iners_lactinv_01v1_a", "core"), "get_DBAdaptor on ensembl genomes species");
  ok($reg->get_DBAdaptor("mus_musculus_balbcj", "core"), "get_DBAdaptor on ensembl species");
  
  destroy_collection_dbs(\@collection_dbs);
}

sub create_collection_dbs {
  my $collection_dbs = 
  [
   {
    name => sprintf("bacteria_9_collection_core_37_%d_1", software_version()),
    species => [[39,2,"species.db_name","escherichia_coli_dh1"],[154,5,"species.db_name","lactobacillus_iners_lactinv_01v1_a"]]
   },
   {
    name => sprintf("mouse_5_collection_core_%d_1", software_version()),
    species => [[27,1,"species.db_name","mus_musculus_aj"],[91,2,"species.db_name","mus_musculus_balbcj"]]
   }
  ];


  # to create meta table 
  my $create_meta_table = <<META;
CREATE TABLE IF NOT EXISTS meta ( 
  meta_id int(11) NOT NULL AUTO_INCREMENT, 
  species_id int(10) unsigned DEFAULT '1', 
  meta_key varchar(40) NOT NULL, 
  meta_value varchar(255) NOT NULL, 
  PRIMARY KEY (meta_id) 
);
META

  my $insert_into_meta_table_prefix = "INSERT INTO meta VALUES";
  
  foreach my $db (@{$collection_dbs}) {
    $dbc->do(sprintf "create database if not exists %s", $db->{name});
    $dbc->do(sprintf "use %s", $db->{name});
    $dbc->do($create_meta_table);

    my $insert_into_meta_table = $insert_into_meta_table_prefix;
    foreach my $meta_value (@{$db->{species}}) {
      $insert_into_meta_table .=
	sprintf " (%d,%d,'%s','%s'),", $meta_value->[0], $meta_value->[1], $meta_value->[2], $meta_value->[3];
    }
    $insert_into_meta_table =~ s/.$//;
    
    $dbc->do($insert_into_meta_table);
  }

  return map { $_->{name} } @{$collection_dbs};
}

sub destroy_collection_dbs {
  my $collection_dbs = shift;
  map { $dbc->do(sprintf "drop database if exists %s", $_) } @{$collection_dbs};
}



done_testing();
