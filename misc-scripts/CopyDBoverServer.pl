#!/usr/local/ensembl/bin/perl -w

use strict;
use Getopt::Long;
use Cwd 'chdir';
use Sys::Hostname qw(hostname);

$| = 1;

my $usage = "\nUsage: $0 -pass XXXXX [-noflush -probe_mapping -xref] input_file

Copy mysql databases between different servers and run myisamchk on the indices when copied.

-dbflush       Flushes only those DBs to be copied, default is to flush entire instance.
-noflush       Skips table flushing
-probe_mapping Only copies the tables relevant to running genomics and transcript probe mapping in isolation
-xref          Only copies the xref tables, to enable running the transcript probe mapping with output to an isolation xref DB.

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
   ia64[efgh] port: 3306
3- -pass is compulsory and is expected to be the mysql password to connect as ensadmin
4- If you add a new instance, remember to check there is a /tmp directory in the /mysql/databases,
   otherwise rcp will complain

This script MUST be run as the mysqlens Unix user.

Also it is best to run it on the same node as the destination MySQL instance, e.g. for ecs2:3365 run

caa_stat mysql_3365

... which will reveal that it's running on ecs2e (caa_stat alone will show all instances). 
Note that this only works for Tru64 CAA instances, i.e. ecs2 & ecs4.
\n\n";


my ($pass, $xref, $probe_mapping, $help, $noflush, $dbflush);


GetOptions(
		   'h'             => \$help,
		   'pass=s'        => \$pass,
           'noflush'       => \$noflush,
		   'probe_mapping' => \$probe_mapping,
		   'xref'          => \$xref,
		   'dbflush'      => \$dbflush,
		  );

if ($help || scalar @ARGV == 0 || ! defined $pass) {
  print $usage;
  exit 0;
}

if($dbflush && $noflush){
  die('Cannot specify mutually exclusive options -noflush and -dbflush');
}

if($probe_mapping && $xref){
  die('Cannot specify mutually exclusive options -probe_mapping and -xref');
}


my ($input_file) = @ARGV;
my @dbs_to_copy;

my %mysql_directory_per_svr = ('genebuild1:3306'   => "/mysql/data_3306/databases",
			       'genebuild2:3306'   => "/mysql/data_3306/databases",
			       'genebuild3:3306'   => "/mysql/data_3306/databases",
			       'genebuild4:3306'   => "/mysql/data_3306/databases",
			       'genebuild5:3306'   => "/mysql/data_3306/databases",
			       'genebuild6:3306'   => "/mysql/data_3306/databases",
							    'genebuild7:3306'   => "/mysql/data_3306/databases",
							   'genebuild7:5306'   => "/mysql/data_3306/databases",
			       'mart1:3306'        => "/mysql/data_3306/databases",
			       'mart2:3306'        => "/mysql/data_3306/databases",
			       'compara1:3306'     => "/mysql/data_3306/databases",
			       'compara2:3306'     => "/mysql/data_3306/databases",
			       'compara2:5316'     => "/mysql/data_5316/databases",
			       'compara3:3306'     => "/mysql/data_3306/databases",
			       'ens-genomics1:3306' => "/mysql/data_3306/databases",
			       'ens-genomics2:3306' => "/mysql/data_3306/databases",
			       'ens-staging:3306'  => "/mysql/data_3306/databases",
							   'ens-staging1:3306'  => "/mysql/data_3306/databases",
							   'ens-staging2:3306'  => "/mysql/data_3306/databases",
							   
			       'ens-livemirror:3306'  => "/mysql/data_3306/databases",
			       'ensdb-archive:3304' => "/mysql/data_3304/databases",
			       'ens-research:3306' => "/mysql/data_3306/databases",
							   'ens-research:3309' => "/mysql/data_3309/databases",
							   #'ensdb-1-11:5317' => '/mysql/data_5317/databases',
							   #mysqlens doesn't have releod privelages here
                               'ensdb-2-12:5106' => '/mysqlv5.1-test/data_5106/databases');

my %tables = (
			  #xref tables are required for transcript/gene generation
			  xref => [('xref', 'object_xref', 'identity_xref', 'go_xref', 
						'external_db', 'external_synonym', 'unmapped_reason', 'unmapped_object')],

			  probe_mapping => [('oligo_array', 'oligo_probe', 'oligo_feature', 'coord_system', 
								 'seq_region', 'seq_region_attrib', 'seq_region_mapping', 'mapping_set', 'assembly_exception', 'attrib_type', 'analysis',
								 'exon', 'exon_stable_id', 'exon_transcript', 'assembly', 'dna',
								 'analysis_description', 'transcript', 'transcript_attrib', 'transcript_stable_id', 
								 'translation', 'translation_stable_id', 'meta', 'meta_coord')],
			  #translation & dna required for generating annotated UTRs
			  			 );

#add xref tables to probe mapping tables
push @{$tables{'probe_mapping'}}, @{$tables{'xref'}};

#Set table sub set
my $table_type = '';
$table_type = 'xref' if $xref;
$table_type = 'probe_mapping' if $probe_mapping;

#Currently this fails if xref or probe_mapping is specified for a non-core DB
#We need to default to normal copy if dbname !~ _core_


my $flush_scope = (defined $dbflush) ? 'src_db' : 'src_srv';
my ($generic_working_host) = (gethostbyname(hostname));
$generic_working_host =~ s/\..*//;
my $working_dir = $ENV{'PWD'};
my %already_flushed;

# parsing/checking the input file

#This does not catch unreadable files!
open F, $input_file ||
  die "Can not open $input_file, $!\n";

while (my $line = <F>) {

  next if ($line =~ /^\#.*$/);
  next if ($line =~ /^\s*$/);

  if ($line =~ /^(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s*$/) {
    my ($src_srv,$src_port,$src_db,$dest_srv,$dest_port,$dest_db) = ($1,$2,$3,$4,$5,$6);
    my ($dest_srv_host) = (gethostbyname($dest_srv));
    $dest_srv_host =~ s/\..*//;

    unless ($generic_working_host =~ /^$dest_srv_host/) {
      warn "// skipped copy of $src_db from $src_srv to $dest_srv
// this script should be run on a generic destination host $dest_srv\n";
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
// this script works only to copy dbs between certain nodes:mysql_port" .
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

my $copy_executable;

if (-e "/usr/bin/cp") {
  $copy_executable = "/usr/bin/cp";
} elsif (-e "/bin/cp") {
  $copy_executable = "/bin/cp";
}

#change STDERR to STDOUT?

# starting copy processes
foreach my $db_to_copy (@dbs_to_copy) {

  print STDERR "//\n// Starting new copy process\n//\n";

  my $time;
  my $source_srv = $db_to_copy->{src_srv};
  $source_srv =~ s/\..*//;
  
  my $source_port = $db_to_copy->{src_port};

  my $source_db = $mysql_directory_per_svr{$source_srv . ":" . $source_port} . "/" . $db_to_copy->{src_db};

  my $destination_srv = $db_to_copy->{dest_srv};
  $destination_srv =~ s/\..*//;
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
  unless (defined $already_flushed{$db_to_copy->{$flush_scope}} || $noflush) {
  

	#Is this sitll flushing the whole instance?
	#Do we need to flush each individual table
	#and reset the query cache?
	print STDERR "// flushing tables in $db_to_copy->{$flush_scope} (".&get_time.")...\t";
	my $flush_cmd;

	if($dbflush){
	  my $tables_cmd = "echo 'show tables;' | mysql -h $db_to_copy->{src_srv} -u ensadmin -p$pass -P$source_port $db_to_copy->{src_db}";	  

	  my @tables = split/\n/, `$tables_cmd`;
	  shift @tables;#remove field header

	  $flush_cmd = "echo 'flush tables ".join(', ', @tables).";' | mysql -h $db_to_copy->{src_srv} -u ensadmin -p$pass -P$source_port $db_to_copy->{src_db}";
	}
	else{
	  $flush_cmd = "echo 'flush tables;' | mysql -h $db_to_copy->{src_srv} -u ensadmin -p$pass -P$source_port";
	  #$flush_cmd .= ' '.$db_to_copy->{src_db} if $dbflush;
	}

	#Now flushes on db specific tables
	#This was introduced to enable copying a an unused DB when others are being heavily used
	#flush tables still reset the query cache, so will this provide an appreciable speed up?
	#Apparently so.  But still slow due to burden on server from other queries
	

    if (system($flush_cmd) == 0) {
      print STDERR "DONE (".&get_time.")\n";
    } else {
      print STDERR "FAILED
skipped copy of ".$db_to_copy->{src_db}." from ".$db_to_copy->{src_srv}." to ". $db_to_copy->{dest_srv} . "\n";
      next;
    }


	#Can we capture a CTRL-C here to exit the whole script
	#Otherwise we exit the flush and the script acrries on trying to copy which can mess up tables.


    $already_flushed{$db_to_copy->{$flush_scope}} = 1;
  }

  # cp the db to $destination_tmp_directory in the destination server
  my $copy_cmd;
  #Need to make the tmp dir first for file copying
  if($table_type && ! -e "$destination_tmp_directory/$db_to_copy->{dest_db}"){
	$copy_cmd = "mkdir $destination_tmp_directory/$db_to_copy->{dest_db};"
  }

  if ($copy_executable =~ /\/bin\/cp$/) {

	if($table_type && ($db_to_copy->{src_db} =~ /_core_/)){#Copy core table subset
	
	  foreach my $table(@{$tables{$table_type}}){
		$copy_cmd .= "$copy_executable  $source_db/$table.* $destination_tmp_directory/$db_to_copy->{dest_db};";
	  }
	  
	}#Full copy
	else{
	  $copy_cmd = "$copy_executable -r $source_db $destination_tmp_directory/$db_to_copy->{dest_db}";
	}
	
	print STDERR "// cp Copying $table_type $db_to_copy->{src_srv}:$source_db...";
 
  } # OR rcp the db to $destination_tmp_directory in the destination server
  elsif ($copy_executable eq "/usr/bin/rcp") {
    
	if($table_type && ($db_to_copy->{src_db} =~ /_core_/)){#Copy core table subset
	
	  foreach my $table(@{$tables{$table_type}}){
		$copy_cmd .= "$copy_executable -r $db_to_copy->{src_srv}:$source_db/$table.* $destination_tmp_directory/$db_to_copy->{dest_db};";
	  }
	  
	}#Full copy
	else{
	  $copy_cmd = "$copy_executable -r $db_to_copy->{src_srv}:$source_db $destination_tmp_directory/$db_to_copy->{dest_db}";
	}

	print STDERR "// rcp Copying $table_type $db_to_copy->{src_srv}:$source_db...";

  }

  if (system("$copy_cmd") == 0) {
    print STDERR "DONE (".&get_time.")\n";
  } else {
    print STDERR "FAILED
skipped copy of $db_to_copy->{src_db} from $db_to_copy->{src_srv} to $db_to_copy->{dest_srv}\n";
print "$copy_cmd\n";
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

print STDERR "//\n// End of all copy processes (".&get_time.")\n//\n// Processes summary\n";


foreach  my $db_to_copy (@dbs_to_copy) {
  print STDERR "// $db_to_copy->{status} copy of $db_to_copy->{src_db} on $db_to_copy->{src_srv} to $db_to_copy->{dest_db} on $db_to_copy->{dest_srv} \n";
}

print STDERR "\n";

sub get_time{

  my ($sec, $min, $hour) = localtime;

  return $hour.':'.$min.':'.$sec;
}

1;
