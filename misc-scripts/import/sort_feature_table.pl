#!/usr/bin/env perl

# Designed to sort a feature table by seq_region_id and then start

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use POSIX qw/strftime/;

my ($db_name,$db_host,$db_user,$db_pass,$db_port,$help);
my @tables;
my ($optimise, $drop);

GetOptions ("db_name|dbname|database=s" => \$db_name,
            "db_host|dbhost|host=s" => \$db_host,
            "db_user|dbuser|user|username=s" => \$db_user,
            "db_pass|dbpass|pass|password=s" => \$db_pass,
            "db_port|dbport|port=s" => \$db_port,
            'table|tables=s@' => \@tables,
            'optimise!' => \$optimise,
            'drop!' => \$drop,
            "h|help!"        => \$help,
);

if ($help) {&usage; exit 0;}
unless ($db_name and $db_host) {print "Insufficient arguments\n"; &usage; exit 1;}

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
  
  info("Creating alternative table and disabling keys");
  $dba->dbc()->do("create table ${s_table} like ${table}");
  $dba->dbc()->do("alter table ${s_table} disable keys");
  
  info("Sort/insert");
  $dba->dbc()->do("insert into ${s_table} select * from ${table} order by seq_region_id, seq_region_start");
  
  info("Re-enabling keys");
  $dba->dbc()->do("alter table ${s_table} enable keys");
  
  my $bak_name = _next_name($table, $dba);
  info("Moving %s to %s", $table, $bak_name);
  $dba->dbc()->do("alter table ${table} rename as $bak_name");
  
  info("Moving %s to %s", $s_table, $table);
  $dba->dbc()->do("alter table ${s_table} rename as $table");
  
  if($drop) {
    info("Dropping table %s", $bak_name);
    $dba->dbc()->do("drop table $bak_name");
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
  print "\n";
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
    
    -drop               Drop the original table post sort
    
    -help
";
}

run();