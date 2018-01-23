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

use FindBin qw/$Bin/;
use Test::More;
use Test::Warnings;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Utils::IO qw/slurp/;
use File::Spec::Functions qw/updir catfile catdir/;
use File::Temp qw/tempfile/;

my $db = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $db->get_DBAdaptor('core');
my $dbc = $dba->dbc();

# Get file location
my $sql_dir = catdir($Bin, updir(), updir(), 'sql');
my $sql_file = catfile($sql_dir, 'table.sql');

my $result = ok(-f $sql_file, 'Checking SQL schema file exists');

SKIP: {
  skip 'Skipping DB creation tests as schema file cannot be found at '.$sql_file, 1 unless $result;
  
  #Create DB & load schema
  my $new_db_name = $db->create_db_name('schematemp');
  note 'Creating database '.$new_db_name;
  $dba->dbc()->do("create database $new_db_name");

  my $table_string = slurp $sql_file;  
  my @fks = sort grep /FOREIGN[\s\n]+?KEY[\s\n]+?\(.+?\)[\s\n]+?REFERENCES/i, $table_string;
  if(@fks) {
    fail("Definition of foreign keys detected in SQL schema file\n".join(', ',@fks));
  } else {
    pass("SQL schema file does not define foreign keys");
  }

  my %args = ( host => $dbc->host(), port => $dbc->port(), user => $dbc->username(), password => $dbc->password());
  my $cmd_args = join(q{ }, map { "--${_}=$args{$_}" } keys %args);
  my $cmd = "mysql $cmd_args $new_db_name < $sql_file 2>&1";
  my $output = `$cmd`;
  my $ec = ($? >> 8);
  if($ec != 0) {
    note($output);
    fail("MySQL command failed with error code '$ec'");
  }
  else {
    pass("MySQL was able to load the Ensembl core schema");
  }
  
  note 'Dropping database '.$new_db_name;
  $dba->dbc()->do("drop database $new_db_name");

  # Check viability of foreign_keys.sql 

  $table_string =~ s/ENGINE=MyISAM/ENGINE=InnoDB/g;
  my $fk_file = catfile($sql_dir, 'foreign_keys.sql');
  $table_string .= slurp $fk_file;

  my ($temp_fh,$temp_filename) = tempfile;
  {
    local $| = 1;
    print $temp_fh $table_string;
  }
  $new_db_name = $db->create_db_name('fkschematemp');
  note 'Creating database with foreign keys '.$new_db_name;
  $dba->dbc()->do("create database $new_db_name");

  %args = ( host => $dbc->host(), port => $dbc->port(), user => $dbc->username(), password => $dbc->password());
  $cmd_args = join(q{ }, map { "--${_}=$args{$_}" } keys %args);
  $cmd = "mysql $cmd_args $new_db_name < $temp_filename 2>&1";
  $output = `$cmd`;
  $ec = ($? >> 8);
  if($ec != 0) {
    note($output);
    fail("MySQL (InnoDB) command failed with error code '$ec'");
  }
  else {
    pass("MySQL (InnoDB) was able to load the Ensembl core schema with foreign keys");
  }
  
  note 'Dropping database '.$new_db_name;
  $dba->dbc()->do("drop database $new_db_name");
}


done_testing();
