#!/usr/local/bin/perl -w

use strict;
use Getopt::Long;

$| = 1;

my $usage = "\nUsage: $0 input_file

Copy mysql databases between different servers and run myisamchk on the indices when copied.

The input file should have the following format

source_server\tsource_database\tdestination_server\tdestination_database

e.g.

#source_server\tsource_db\tdestination_server\tdestination_db
ecs3d.internal.sanger.ac.uk     homo_sapiens_core_13_31 ecs2d.internal.sanger.ac.uk     homo_sapiens_core_14_31

Lines starting with # are ignored and considered as comments.

RESTRICTIONS:
============
1- The destination_server is also the one on which the script has to be run,
   otherwise the copying process for the corresponding database is skipped
2- This script works only for copy processes from and to ecs nodes, namely
   ecs1[abcdefgh]
   ecs2[abcdef]
   ecs3d only

";

my ($input_file) = @ARGV;
my $help = 0;

GetOptions('h' => \$help);

if ($help || scalar @ARGV == 0) {
  print $usage;
  exit 0;
}

my @dbs_to_copy;

my %mysql_directory_per_svr = ('ecs1a' => "/mysql1a",
			       'ecs1b' => "/mysql2a",
			       'ecs1c' => "/mysql3a",
			       'ecs1d' => "/mysql4a",
			       'ecs1e' => "/mysql5a",
			       'ecs1f' => "/mysql6a",
			       'ecs1g' => "/mysql7a",
			       'ecs1h' => "/mysql_archive",
			       'ecs2a' => "/mysqla",
			       'ecs2b' => "/mysqlb",
			       'ecs2c' => "/mysqlc",
			       'ecs2d' => "/mysqld",
			       'ecs2e' => "/mysqle",
			       'ecs2f' => "/mysqlf",
			       'ecs3d' => "/mysqld");

open F, $input_file ||
  die "Can not open $input_file, $!\n";

while (my $line = <F>) {
  next if ($line =~ /^\#.*$/);
  if ($line =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+)$/) {
    my ($src_srv,$src_db,$dest_srv,$dest_db) = ($1,$2,$3,$4);
    my %hash = ('src_srv' => $src_srv,
		'src_db' => $src_db,
		'dest_srv' => $dest_srv,
		'dest_db' => $dest_db);
    push @dbs_to_copy, \%hash;
  } else {
    warn "
The input file has the wrong format.
source_server\tsource_db\tdestination_server\tdestination_db
EXIT 1
";
    exit 1;
  }
}

close F;

my $working_host = $ENV{'HOST'};
my $copy_executable = "/usr/bin/cp";

foreach my $db_to_copy (@dbs_to_copy) {
  print STDERR "//Starting new copy process\n";
  unless ($db_to_copy->{dest_srv} =~ /^$working_host.*$/) {
    warn "WARN: skipped copy of ".$db_to_copy->{src_db}." from ".$db_to_copy->{src_srv}." to ". $db_to_copy->{dest_srv} . ".
WARN: This script should be run on the destination host ". $db_to_copy->{dest_srv} ."\n";
    next;
  }

  my $source_srv = $db_to_copy->{src_srv};
  $source_srv =~ s/(ecs[123].{1}).*/$1/;
  my $destination_srv = $db_to_copy->{dest_srv};
  $destination_srv =~ s/(ecs[123].{1}).*/$1/;
  
  unless (defined $mysql_directory_per_svr{$source_srv} && 
	  defined $mysql_directory_per_svr{$destination_srv}) {
    warn "This script works only to copy dbs between ecs nodes
$source_srv $destination_srv
EXIT 2
";
    exit 2;
  }
  my $source_db_directory = $mysql_directory_per_svr{$source_srv}."/current/var/".$db_to_copy->{src_db};
  my $destination_db_directory = $mysql_directory_per_svr{$destination_srv}."/current/var/".$db_to_copy->{dest_db};
  my $myisamchk_executable = $mysql_directory_per_svr{$destination_srv}."/current/bin/myisamchk";
    
  $source_srv =~ s/(ecs[123]).*/$1/;
  $destination_srv =~ s/(ecs[123]).*/$1/;

  if ($source_srv eq $destination_srv) {
    $copy_executable = "/usr/bin/cp";
  } else {
    $copy_executable = "/usr/bin/rcp";
  }

  print STDERR "Creating $destination_db_directory directory...";
  if (system("mkdir $destination_db_directory") == 0) {
    print STDERR "Done\n";
  } else {
    print STDERR "Failed
skipped copy of ".$db_to_copy->{src_db}." from ".$db_to_copy->{src_srv}." to ". $db_to_copy->{dest_srv} . "\n";
    next;
  }

  if ($copy_executable eq "/usr/bin/cp") {
    print STDERR "cp Copying $db_to_copy->{src_srv}:$source_db_directory/*...";
    my $command_line = "$copy_executable -r $source_db_directory/* $destination_db_directory";
    if (system("$command_line") == 0) {
      
    } else {
      print STDERR "Failed
skipped copy of ".$db_to_copy->{src_db}." from ".$db_to_copy->{src_srv}." to ". $db_to_copy->{dest_srv} . "\n";
      next;
      
    }
  } elsif ($copy_executable eq "/usr/bin/rcp") {
    print STDERR "rcp Copying $db_to_copy->{src_srv}:$source_db_directory/*...";
    my $command_line = "$copy_executable -r $db_to_copy->{src_srv}:$source_db_directory/* $destination_db_directory";
    
    if (system("$command_line") == 0) {
      print STDERR "Done\n";
    } else {
      print STDERR "Failed
skipped copy of ".$db_to_copy->{src_db}." from ".$db_to_copy->{src_srv}." to ". $db_to_copy->{dest_srv} . "\n";
      next;
      
    }
  }
  
  print STDERR "Checking $destination_db_directory/*.MYI in progress\n";
  if (system("$myisamchk_executable -r $destination_db_directory/*.MYI") == 0) {
    print STDERR "Checking $destination_db_directory/*.MYI Done\n";
  } else {
    print STDERR "Checking $destination_db_directory/*.MYI Failed
skipped checking of " . $db_to_copy->{dest_srv} . "\n";
    next;
  }
}

print STDERR "//End of all copy processes\n";
