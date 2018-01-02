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

# We respond to the following ENV
#  ENS_REMOTE_SCHEMA - If set to true we will use the remote schema as the reference schema. Otherwise we use the local schema held in sql/table.sql

use strict;
use warnings;

use Test::More;
use Test::Warnings;
use Test::Differences;

use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Utils::Net qw/do_GET/;
use Bio::EnsEMBL::Utils::IO qw/slurp/;
use Bio::EnsEMBL::Test::MultiTestDB;
use File::Find;
use File::Spec::Functions qw/updir catfile catdir/;
use File::Temp qw/tempfile/;
use FindBin qw/$Bin/;

SKIP: {

  # Get last DB version
  my $current_release = software_version();
  my $last_release = $current_release - 1;

  # Get patch location and relevant set of patches
  my $project_dir = catdir($Bin, updir(), updir());
  my $sql_dir = catdir($project_dir, 'sql');
  my @patches;
  find(sub {
    if($_ =~ /^patch_${last_release}_${current_release}[_-]?\w+\.sql$/) {
      push(@patches, $File::Find::name);
    }
  }, $sql_dir);

  # Get the last SQL schema
  my $last_table_sql = eval { get_table_sql($last_release); };

  skip "Skipping DB patch tests as we cannot fetch the SQL for release $last_release", (scalar(@patches)+1) 
    if $@;

  my $db = Bio::EnsEMBL::Test::MultiTestDB->new();
  my $dba = $db->get_DBAdaptor('core');
  my $dbc = $dba->dbc();
  
  # Create last release DB
  my $patched_db_name = $db->create_db_name('schemapatchestemp');
  note 'Creating database ' . $patched_db_name;
  $dba->dbc()->do("create database $patched_db_name");
  
  # Load last release's schema
  my ($fh, $sql_schema_file) = tempfile();
  print $fh $last_table_sql;
  close $fh;
  my $loaded_schema = load_sql($dbc, $patched_db_name, $sql_schema_file);
  
  skip 'Skipping DB patch tests as we cannot load the last release schema into a database', scalar(@patches) 
    unless $loaded_schema;

  # Create last release DB  
  my $current_table_sql;
  if ($ENV{ENS_REMOTE_SCHEMA}) {
    $current_table_sql = eval { get_table_sql($current_release); };
    skip 'Skipping DB patch test as we cannot fetch current release schema file (table.sql)', 1
      if $@;
  } else {
    my $table_sql = catfile($sql_dir, 'table.sql');
    skip 'Skipping DB patch test as we cannot find current release schema file (table.sql)', 1
      unless -e $table_sql;
    skip 'Skipping DB patch test as we current release schema file (table.sql) is not readable', 1
      unless -r $table_sql;
    $current_table_sql = slurp($table_sql);
  }
  
  skip "Skipping DB patch tests as we cannot find the SQL for release $current_release", (scalar(@patches)+1) 
    unless defined $current_table_sql;

  my $current_db_name = $db->create_db_name('schematemp');
  note 'Creating database ' . $current_db_name;
  $dba->dbc()->do("create database $current_db_name");

  # Load current release's schema
  ($fh, $sql_schema_file) = tempfile();
  print $fh $current_table_sql;
  close $fh;
  $loaded_schema = load_sql($dbc, $current_db_name, $sql_schema_file);
  
  skip 'Skipping DB patch tests as we cannot load current release schema into a database', scalar(@patches) 
    unless $loaded_schema;

  # Now apply all current patches in order
  foreach my $patch (sort @patches) {
    # Get the number of patch entries before applying next patch
    my $previous_patches = get_num_patches($dbc, $patched_db_name);

    note "Applying patch $patch";
    load_sql($dbc, $patched_db_name, $patch);
    check_after_patch($dbc, $patched_db_name, $previous_patches);
  }

  # check the two schemas after applying the patch
  compare_after_patches($dbc, $patched_db_name, $current_db_name, $last_release, $current_release);
  
  note 'Dropping database ' . $patched_db_name;
  $dba->dbc()->do("drop database if exists $patched_db_name");
  note 'Dropping database ' . $current_db_name;
  $dba->dbc()->do("drop database if exists $current_db_name");
}

done_testing();

# Compare source schema with target after a series of patches
sub compare_after_patches {
  my ($dbc, $source_schema, $target_schema, $last_release, $current_release) = @_;

  # compare source/target schema type/version
  $dbc->do("use $target_schema");
  my $sql_helper = $dbc->sql_helper;
  my ($target_schema_type, $target_schema_version) =
    ($sql_helper->execute_single_result(-SQL => "select meta_value from meta where meta_key='schema_type'"),
     $sql_helper->execute_single_result(-SQL => "select meta_value from meta where meta_key='schema_version'"));

  $dbc->do("use $source_schema");
  my ($source_schema_type, $source_schema_version) = 
    ($sql_helper->execute_single_result(-SQL => "select meta_value from meta where meta_key='schema_type'"),
     $sql_helper->execute_single_result(-SQL => "select meta_value from meta where meta_key='schema_version'"));

  is($source_schema_type, $target_schema_type, "Schema type after patches");
  is($source_schema_version, $target_schema_version, "Schema version after patches");
  
  # Check if the patch meta value does not contain line breaks
  my $patch_meta_values_with_newlines = 
    $sql_helper->execute_simple(-SQL => "select meta_value from meta where meta_key='patch' and meta_value like '%\n%'");
  is(scalar @{$patch_meta_values_with_newlines}, 0, "No line breaks in patch meta values");

  # check the two schemas contain the same tables
  my $source_tables = get_table_names($dbc, $source_schema);
  my $target_tables = get_table_names($dbc, $target_schema);
  my $diff = (union_intersection_difference($source_tables, $target_tables))[2];
  is(scalar @{$diff}, 0, "Same table set");
  
  # check each table has the same definition in both schemas
  map { eq_or_diff(get_create_table($dbc, $source_schema, $_), 
	   get_create_table($dbc, $target_schema, $_),
	   "Table $_ definition")} 
    @{$source_tables};

  # get target patches
  $dbc->do("use $target_schema");
  my $target_patches = $sql_helper->execute_simple(-SQL => "select meta_value from meta where meta_key='patch'");

  # get target-specific patches in source
  $dbc->do("use $source_schema");
  my $source_patches = $sql_helper->execute_simple(-SQL => "select meta_value from meta where meta_key='patch' and meta_value like 'patch_${last_release}_${current_release}_%'");

  my %source_patches;
  my %target_patches;
  map { $source_patches{$_}++ } @{$source_patches};
  map { $target_patches{$_}++ } @{$target_patches};
  map { ok(exists $source_patches{$_}, "$_ in patched database") } @{$target_patches};
  map { ok(exists $target_patches{$_}, "$_ in original database") } @{$source_patches};



}

# Get the name of all tables of a certain schema
sub get_table_names {
  my ($dbc, $schema_name) = @_;
  $dbc->do("use $schema_name");

  return 
    $dbc->sql_helper->execute_simple(-SQL => 'show tables');
}

# Get the create table SQL statement
sub get_create_table {
  my ($dbc, $schema_name, $table_name) = @_;
  $dbc->do("use $schema_name");
  my $sql_helper = $dbc->sql_helper;

  my $create_table = $sql_helper->execute(
    -SQL => "show create table $table_name",
    -CALLBACK => sub {
      return (shift @_)->[1];
    }
  )->[0];
  
  # stripping AUTOINCREMENT=? definitions in the way since 
  # they are allowed to be different
  $create_table =~ s/AUTO_INCREMENT=\d+//;

  return $create_table;
}

# Check source schema after applying one patch
sub check_after_patch {
  my ($dbc, $source_schema, $num_previous_patches) = @_;

  # see if after patch we gain a meta key corresponding to the applied patch
  my $num_current_patches = get_num_patches($dbc, $source_schema);
  is($num_current_patches, $num_previous_patches + 1, "Source schema gains patch meta key after patch");
}

# Get the number of patches applied
sub get_num_patches {
  my ($dbc, $schema) = @_;
  $dbc->do("use $schema");
  my $sql_helper = $dbc->sql_helper;

  return 
    $sql_helper->execute_single_result(-SQL => "select count(*) from meta where meta_key='patch' and species_id is NULL");
}

# Load SQL subroutine
sub load_sql {
  my ($dbc, $dbname, $path) = @_;
  my %args = ( host => $dbc->host(), port => $dbc->port(), user => $dbc->username(), password => $dbc->password());
  my $cmd_args = join(q{ }, map { "--${_}=$args{$_}" } keys %args);
  my $cmd = "mysql $cmd_args $dbname < $path 2>&1";
  my $output = `$cmd`;
  my $ec = ($? >> 8);
  if($ec != 0) {
    note($output);
    return fail("MySQL command failed with error code '$ec'");
  }
  return pass("MySQL was able to load the file $path core schema");    
}

# Get table.sql for a given Ensembl release
sub get_table_sql {
  my $release = shift;
  $release = $release == software_version() ? 'master' : "release/${release}";
  
  my $url = "https://raw.githubusercontent.com/Ensembl/ensembl/${release}/sql/table.sql";
  return get_url($url);
}

# Assume if ensembl.org is down then there is no point in continuing with the tests (1 shot)
sub test_ensembl {
  my $content = eval { do_GET('http://www.ensembl.org', 1); };
  my $success = 1;
  if($@) {
    note 'ensembl.org is unreachable. Cannot continue with tests';
    $success = 0;
  }
  return $success;
}

sub get_url {
  my ($url) = @_;
  my $content = eval { do_GET($url, 5, 0.5); };
  return $content if defined $content;
  diag $@;
  fail("We do not have access to HTTP::Tiny or LWP. Cannot continue") if $@;
  return;
}

# Computing Union, Intersection, or Difference of Unique Lists
sub union_intersection_difference {
  my ($a, $b) = @_;

  my (@union, @isect, @diff);
  my %count = ();
  foreach my $e (@$a, @$b) { $count{$e}++ }
  
  foreach my $e (keys %count) {
    push(@union, $e);
    if ($count{$e} == 2) {
      push @isect, $e;
    } else {
      push @diff, $e;
    }
  }
  return (\@union, \@isect, \@diff);
}
