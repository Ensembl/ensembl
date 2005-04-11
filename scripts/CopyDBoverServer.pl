#!/usr/local/ensembl/bin/perl -w

use strict;
use Getopt::Long;
use Cwd 'chdir';

$| = 1;

my $usage = "\nUsage: $0 -pass XXXXX [-noflush] --nochk] input_file

Copy mysql databases between different servers and run myisamchk on the indices when copied.

The input file should have the following format

source_server\\tsource_port\\tsource_db\\tdestination_server\\tdestination_port\\tdestination_db

e.g.

#source_server\\tsource_port\\tsource_db\\tdestination_server\\tdestination_port\\tdestination_db
ecs3.internal.sanger.ac.uk 3307 homo_sapiens_core_13_31 ecs2.internal.sanger.ac.uk 3364 homo_sapiens_core_14_31

Lines starting with # are ignored and considered as comments.
Blank (or whitespace-only) lines are ignored.

Note fields can be separated by any number of tabs or spaces.

RESTRICTIONS:
============
1- The destination_server has to match the generic server you are running the script on,
   either ecs1, ecs2, ecs3, ecs4 or ia64[ef] otherwise the copying process for the corresponding database
   is skipped
2- This script works only for copy processes from and to ecs/ia64 nodes, namely
   ecs1[abcdefgh] port 3306
   ecs2 port: 336[1-6]
   ecs3 port: 300[47]
   ecs4 port: 335[0-3]
   ia64[ef] port: 3306
3- -pass is compulsory and is expected to be the mysql password to connect as ensadmin

This script MUST be run as the mysqlens Unix user.

Also it is best to run it on the same node as the destination MySQL instance, e.g. for ecs2:3365 run

caa_stat mysql_3365

... which will reveal that it's running on ecs2e (caa_stat alone will show all instances). 
Note that this only works for Tru64 CAA instances, i.e. ecs2 & ecs4.
\n\n";

my $help = 0;
my $pass;
my $noflush = 0;

GetOptions('h' => \$help,
	   'pass=s' => \$pass,
           'noflush' => \$noflush);

if ($help || scalar @ARGV == 0 || ! defined $pass) {
  print $usage;
  exit 0;
}

my ($input_file) = @ARGV;
my @dbs_to_copy;

	my %mysql_directory_per_svr = ('ecs1a:3306' => "/mysql1a/databases",
                               'ecs1b:3306' => "/mysql2a/databases",
                               'ecs1c:3306' => "/mysql3a/databases",
                               'ecs1d:3306' => "/mysql4a/databases",
                               'ecs1e:3306' => "/mysql5a/databases",
                               'ecs1f:3306' => "/mysql6a/databases",
                               'ecs1g:3306' => "/mysql7a/databases",
                               'ecs1h:3306' => "/mysql_archive/current/var",
			       'ecs2:3361' => "/mysql/data_3361/databases",
			       'ecs2:3362' => "/mysql/data_3362/databases",
			       'ecs2:3363' => "/mysql/data_3363/databases",
			       'ecs2:3364' => "/mysql/data_3364/databases",
			       'ecs2:3365' => "/mysql/data_3365/databases",
			       'ecs2:3366' => "/mysql/data_3366/databases",
			       'ecs3:3307' => "/mysql/current/var",
			       'ecs3:3309' => "/mysqlh/current/var",
			       'ecs3:3304' => "/mysql_archive/current/var",
			       'ecs4:3350' => "/mysql-3350/databases",
			       'ecs4:3351' => "/mysql-3351/databases",
			       'ecs4:3352' => "/mysql-3352/databases",
			       'ecs4:3353' => "/mysql-3353/databases",
			       'ia64e:3306' => "/mysql/data_3306/databases",
                               'ia64f:3306' => "/mysql/data_3306/databases",
                               'ia64g:3306' => "/mysql/data_3306/databases",
                               'ia64h:3306' => "/mysql/data_3306/databases");

my $working_host = $ENV{'HOST'};
my $generic_working_host = $working_host;
$generic_working_host =~ s/(ecs[1234]).*/$1/;
my $working_dir = $ENV{'PWD'};
my %already_flushed;

# parsing/checking the input file

open F, $input_file ||
  die "Can not open $input_file, $!\n";

while (my $line = <F>) {
  next if ($line =~ /^\#.*$/);
  next if ($line =~ /^\s*$/);
  if ($line =~ /^(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s*$/) {
    my ($src_srv,$src_port,$src_db,$dest_srv,$dest_port,$dest_db) = ($1,$2,$3,$4,$5,$6);
    unless ($dest_srv =~ /^$generic_working_host.*$/) {
      my $generic_destination_server = $dest_srv;
      $generic_destination_server =~ s/(ecs[1234]).*/$1/;
      warn "// skipped copy of $src_db from $src_srv to $dest_srv
// this script should be run on a generic destination host $generic_destination_server\n";
      next;
    }
    my $src_srv_ok = 0;
    my $dest_srv_ok = 0;
    foreach my $available_srv_port (keys %mysql_directory_per_svr) {
      my ($srv,$port) = split ":", $available_srv_port;
      if ($src_srv =~ /^$srv.*$/ && $src_port == $port) {
	$src_srv_ok = 1;
      }
      if ($dest_srv =~ /^$srv.*$/ && $dest_port == $port) {
	$dest_srv_ok = 1;
      }
    }
    unless ($src_srv_ok && $dest_srv_ok) {
      warn "// skipped copy of $src_db from $src_srv to $dest_srv
// this script works only to copy dbs between certain ecs_nodes:mysql_port" .
join(", ", keys %mysql_directory_per_svr) ."\n";
      next;
    }
    my %hash = ('src_srv' => $src_srv,
		'src_db' => $src_db,
		'src_port' => $src_port,
		'dest_srv' => $dest_srv,
		'dest_db' => $dest_db,
		'dest_port' => $dest_port,
		'status' => "FAILED");
    push @dbs_to_copy, \%hash;
  } else {
    warn "
The input file has the wrong format,
$line
source_server\\tsource_port\\tsource_db\\tdestination_server\\tdestination_port\\tdestination_db
EXIT 1
";
    exit 1;
  }
}

close F;

# starting copy processes
foreach my $db_to_copy (@dbs_to_copy) {
  my $copy_executable;
  if (-e "/usr/bin/cp") {
    $copy_executable = "/usr/bin/cp";
  } elsif (-e "/bin/cp") {
    $copy_executable = "/bin/cp";
  }
  print STDERR "//
// Starting new copy process
//\n";

  my $source_srv = $db_to_copy->{src_srv};
  $source_srv =~ s/(ecs[1234][a-h]?)\.*.*/$1/;
  $source_srv =~ s/(ia64[ef])\.*.*/$1/;
  
  my $source_port = $db_to_copy->{src_port};

  my $source_db = $mysql_directory_per_svr{$source_srv . ":" . $source_port} . "/" . $db_to_copy->{src_db};

  my $destination_srv = $db_to_copy->{dest_srv};
  $destination_srv =~ s/(ecs[1234][a-h]?)\.*.*/$1/;
  my $destination_port = $db_to_copy->{dest_port};

  my $destination_directory = $mysql_directory_per_svr{$destination_srv . ":" . $destination_port};
  
  my $destination_tmp_directory = $destination_directory;
  $destination_tmp_directory =~ s/\/var//;
  $destination_tmp_directory =~ s/\/databases//;
  $destination_tmp_directory .= "/tmp";

  # checking that destination db does not exist
  if (-e "$destination_directory/$db_to_copy->{dest_db}") {
    print STDERR "// $destination_directory/$db_to_copy->{dest_db} already exists, make sure to
// delete it or use another destination name for the database
// Skipped copy of $db_to_copy->{src_db} from $db_to_copy->{src_srv} to $db_to_copy->{dest_srv}
";
    next;
  }
  
  my $myisamchk_executable = "/usr/local/ensembl/mysql/bin/myisamchk";
  
  $source_srv =~ s/(ecs[1234]).*/$1/;
  $destination_srv =~ s/(ecs[1234]).*/$1/;

  if ($source_srv ne $destination_srv) {
    $copy_executable = "/usr/bin/rcp";
  }

  $source_srv = undef;
  $destination_srv = undef;

  # flush tables; in the source server
  unless (defined $already_flushed{$db_to_copy->{src_srv}} || $noflush) {
    print STDERR "// flushing tables in $db_to_copy->{src_srv}...";
    my $flush_cmd = "echo \"flush tables;\" | mysql -h $db_to_copy->{src_srv} -u ensadmin -p$pass -P$source_port";
    if (system($flush_cmd) == 0) {
      print STDERR "DONE\n";
    } else {
      print STDERR "FAILED
skipped copy of ".$db_to_copy->{src_db}." from ".$db_to_copy->{src_srv}." to ". $db_to_copy->{dest_srv} . "\n";
      next;
    }
    $already_flushed{$db_to_copy->{src_srv}} = 1;
  }

  # cp the db to $destination_tmp_directory in the destination server
  my $copy_cmd;
  if ($copy_executable =~ /\/bin\/cp$/) {
    print STDERR "// cp Copying $db_to_copy->{src_srv}:$source_db...";
    $copy_cmd = "$copy_executable -r $source_db $destination_tmp_directory/$db_to_copy->{dest_db}";
    
  # OR rcp the db to $destination_tmp_directory in the destination server
  } elsif ($copy_executable eq "/usr/bin/rcp") {
    print STDERR "// rcp Copying $db_to_copy->{src_srv}:$source_db...";
    $copy_cmd = "$copy_executable -r $db_to_copy->{src_srv}:$source_db $destination_tmp_directory/$db_to_copy->{dest_db}";
  }

  if (system("$copy_cmd") == 0) {
    print STDERR "DONE\n";
  } else {
    print STDERR "FAILED
skipped copy of $db_to_copy->{src_db} from $db_to_copy->{src_srv} to $db_to_copy->{dest_srv}\n";
    next;
  }

  # checks/fixes the indices
  print STDERR "// Checking $destination_tmp_directory/$db_to_copy->{dest_db}/*.MYI in progress...
//\n";
  chdir "$destination_tmp_directory/$db_to_copy->{dest_db}";
  my $myisamchk_cmd = "ls | grep MYI | xargs $myisamchk_executable -F -f -s --key_buffer_size=2000000000 --sort_buffer_size=2000000000 --read_buffer_size=2000000 --write_buffer_size=2000000";
  if (system("$myisamchk_cmd") == 0) {
    print STDERR "//
// Checking $destination_tmp_directory/$db_to_copy->{dest_db}/*.MYI DONE\n";
    chdir "$working_dir";
  } else {
    print STDERR "//
// Checking $destination_tmp_directory/$db_to_copy->{dest_db}/*.MYI FAILED
skipped checking/copying of $db_to_copy->{dest_db}\n";
    system("rm -rf $destination_tmp_directory/$db_to_copy->{dest_db}");
    chdir "$working_dir";
    next;
  }

  # moves db to mysql directory if checking went fine, skip otherwise
  if (system("mv $destination_tmp_directory/$db_to_copy->{dest_db} $destination_directory") == 0) {
    print STDERR "// moving $destination_tmp_directory/$db_to_copy->{dest_db} to $destination_directory DONE\n";
  } else {
    print STDERR "// moving $destination_tmp_directory/$db_to_copy->{dest_db} to $destination_directory FAILED\n";
    system("rm -rf $destination_tmp_directory/$db_to_copy->{dest_db}");
    next;
  }
  $db_to_copy->{status} = "SUCCEEDED";
}

print STDERR "//
// End of all copy processes
//
// Processes summary\n";


foreach  my $db_to_copy (@dbs_to_copy) {
  print STDERR "// $db_to_copy->{status} copy of $db_to_copy->{src_db} on $db_to_copy->{src_srv} to $db_to_copy->{dest_db} on $db_to_copy->{dest_srv} \n";
}

print STDERR "\n";
