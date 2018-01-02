#!/usr/bin/env perl
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


# Designed to sort a feature table by seq_region_id, start and end

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use POSIX qw/strftime/;

my ($db_name,$db_host,$db_user,$db_pass,$db_port,$help);
my @tables;
my @columns;
my ($optimise, $nobackup, $nolock);

GetOptions ("db_name|dbname|database=s" => \$db_name,
            "db_host|dbhost|host=s" => \$db_host,
            "db_user|dbuser|user|username=s" => \$db_user,
            "db_pass|dbpass|pass|password=s" => \$db_pass,
            "db_port|dbport|port=s" => \$db_port,
            'table|tables=s@' => \@tables,
            'optimise!' => \$optimise,
            'nolock!' => \$nolock,
            'nobackup!' => \$nobackup,
            'columns=s@' => \@columns,
            "h|help!"        => \$help,
);

if ($help) {&usage; exit 0;}
unless ($db_name and $db_host) {print "Insufficient arguments\n"; &usage; exit 1;}

if(!@columns) {
  @columns = qw/seq_region_id seq_region_start seq_region_end/;
}

sub get_adaptor {
  my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -species => 'tmp',
    -group => 'db',
    -dbname => $db_name,
    -host => $db_host,
    -user => $db_user,
    -port => $db_port
  );
  $dba->dbc->password($db_pass) if $db_pass;
  return $dba;
}

# see bottom of file for this method call
sub run {
  my $dba = get_adaptor();
  foreach my $table (@tables) {
    info('Processing %s', $table);
    sort_table($table, $dba);
    optimise_table($table, $dba) if $optimise;
  }
  return;
}

sub sort_table {
  my ($table, $dba) = @_;
  info("Starting sort");
  my $s_table = $table.'_sorted';
  
  if(!$nolock) {
    info("Locking table %s", $table);
    $dba->dbc()->do("lock tables ${table} write");
  }
  
  if(!$nobackup) {
    backup_table($table); 
  }
  
  info("Re-ordering table");
  my $cols = join(',', @columns);
  $dba->dbc()->do("ALTER TABLE ${table} ORDER BY $cols");

  if(!$nolock) {
    info("Unlocking table %s", $table);
    $dba->dbc()->do("unlock tables");
  }
  
  info("Done");
  return;
}

sub optimise_table {
  my ($table, $dba) = @_;
  info("Optimising table");
  $dba->dbc()->do("OPTIMIZE TABLE ${table}");
  info("Done");
  return;
}

sub backup_table {
  my ($table, $dba) = @_;
  my $bak_name = _next_name($table, $dba);
  info("Backing up table to $bak_name");
  $dba->dbc()->do("create table ${bak_name} like ${table}");
  $dba->dbc()->do("alter table ${bak_name} disable keys");
  $dba->dbc()->do("insert into ${bak_name} select * from ${table}");
  $dba->dbc()->do("alter table ${bak_name} enable keys");
  info("Finished backing up");
  return $bak_name;
}

sub _next_name {
  my ($table, $dba) = @_;
  my $count = 0;
  my $new_name;
  while (1) {
    $new_name = $count == 0 ? "${table}_bak" : "${table}_bak_${count}";
    my $r = $dba->dbc()->sql_helper()->execute(
      -SQL => 'show tables like ?', -PARAMS => [$new_name]
    );
    last if scalar(@{$r}) == 0;
    $count++;
  }
  return $new_name;
}

sub info {
  my ($msg, @args) = @_;
  my $m = sprintf $msg, @args;
  my $time = strftime('%c',localtime());
  printf STDERR '[%s] %s', $time, $m;
  print STDERR "\n";
  return;
}

sub usage {
    print "Description:

Sort a feature table by seq_region_id, seq_region_start. Also run optimise
if specified. Holds a backup of the table that has just been sorted
unless told otherwise.

Synopsis:

  perl sort_feature_table.pl -db_name NAME

Options:
    
    -db_name            The DB to add these features to
    -database
    -dbname
    
    -db_host            Hostname for the DB
    -host
    -dbhost
    
    -db_user            Username for the DB
    -user
    -username
    -dbuser
    
    -db_pass            Password for the DB
    -pass
    -password
    -dbpass
    
    -db_port            Port for the DB
    -dbport
    -port
    
    -table              The table(s) to sort. Specify multiple parameters
    -tables
    
    -optimise           Optimise the table post sort
    
    -nobackup           Do not backup the original table
    
    -nolock             Stop the code from applying for table locks
    
    -columns            Specify the columns to sort on. Defaults to
                        seq_region_id, seq_region_start and seq_region_end.
                        Multiple parameters allowed
    
    -help
";
}

run();
