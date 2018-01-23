#!/usr/bin/perl
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

use Getopt::Long;

# download and import the database
my $database;
my $root;
my $new_database;
my $user;
my $pass;
my $port;
my $host;
my $cleanup = undef;
my $force = undef; # if set ignore checksum dies just wrtie warnings.
my $mysqltmpdir = undef;
my $quiet = 0;

GetOptions ('root=s'                    => \$root,
            'database=s'                => \$database,
	    'new_database=s'            => \$new_database,
	    'host=s'                    => \$host,
            'force'                     => \$force,
	    'cleanup'                   => \$cleanup,
            'port=s'                    => \$port,
            'user=s'                    => \$user,
            'pass=s'                    => \$pass,
	    'mysqltempdir=s'            => \$mysqltmpdir,
	    'quiet'                     => \$quiet,
	    'help'  => sub { usage(); exit(0);}  
	   );


if(defined($database)){
  if(!defined($root)){
    #query database to try and guess root;
    $database =~ /\S+_\S+_\S+_(\d+)_/;
    my $release = $1;
    if(defined($release)){
      $root = "//ftp.ensembl.org/ensembl/pub/release-".$release."/mysql";
      print "Using $root as the root obtained from the database name\n" unless $quiet;
    }
    else{
      die "No root given i.e. ftp.ensembl.org/pub/release-54/mysql and could not guess from the database name $database";
    }
  }
}

if(!defined($root)){
  die "No root given i.e. ftp.ensembl.org/pub/release-54/mysql and no database name given to try and guess root from";
}

if(!defined($new_database)){
  $new_database = $ENV{"USER"}."_".$database;
  print "will create new database $new_database\n" unless $quiet;
}

if(!defined $user or !defined $pass or !defined $host){
  die "Need user, password and host for mysql instance to create new database on\n";
}


my $mysql_options  = "-h$host -u$user -p$pass";
if(defined($port)){
  $mysql_options .= " -P$port";
}

print "rsync --recursive rsync:$root/$database .\n" unless $quiet;
my $line;
#goto SKIP;
if($quiet){
  $line = `rsync --recursive --verbose rsync:$root/$database .`;
}
else{
  $line = `rsync --recursive --quiet rsync:$root/$database .`;
}

print $line unless $quiet;
#SKIP:

#if it does snot exist then so be it just ignore error code
#my $com = "mysql $mysql_options -e'drop database ".$new_database."'";
#$line = `$com`;
# no need to check here as if the databae does not exist it should get an error
# just done to delete if it exists already


##
## generate error to test
##
#$mysql_options =~ s/-uensadmin/-uensro/g;

my $com = "mysql $mysql_options -e'create database $new_database'";
$line = `$com`;
if($? or $line =~ /Error/ or $line =~ /ERROR/){
  print $line;
  die "Error during mysql\n";
}
else{
  print "Created new database $new_database on host $host\n" unless $quiet;
}


$mysql_options .= " $new_database";


#get the database schema and load it.
print "now creating the schema\n" unless $quiet;
system("gunzip  -f $database/$database.sql.gz");
system("mysql $mysql_options < $database/$database.sql");
system("gzip  $database/$database.sql");

system("gunzip -f $database/CHECKSUMS.gz");
print "now parse the checksum\n" unless $quiet;

if(defined($mysqltmpdir)){
  $mysql_options = " --tmpdir $mysqltmpdir ".$mysql_options;
}

open(CHK,"<$database/CHECKSUMS") or die "COuld not open CHECKSUMS for reading???\n";
while (<CHK>){
  chomp;
  my ($e1, $e2, $file) = split;
  my $table;
  my $index = "";
  if($file =~ /(\S+)(.*\d*).txt.gz/){
    $table = $1;
    $index = $2;
  }
  else{
    print "ignore $file\n" unless $quiet;
    next;
  }
  if(!-e "$database/$file"){
    print STDERR "$database/$file does not exist. It is specified in the CHECKSUM file but cannot be found?";
    cleanup(1)
  }	
  $com = "sum $database/$file";
  $line = `$com`;
  if($?){
    print STDERR "$com failed\n";
    print STDERR "with output:".$line."\n";
    print STDERR "and error code $?\n";
    print STDERR "Ending as no checksum could be obtained";
    cleanup(1);
  }
  my ($s1, $s2, @junk) = split (/\s+/,$line);
  if($s1 != $e1 or $s2 != $e2){
    print STDERR "Warning: checksums do not match for file $database/$file\n" unless $quiet;
    print STDERR "         from checksum we have $e1 and $e2\n" unless $quiet;
    print STDERR "         but from sum  we have $s1 and $s2\n" unless $quiet;
    if(defined($force)){
      print "   Force set so carrying on\n" unless $quiet;
    }
    else{
      print STDERR "Checksums do not match which can be a problem.\n";
      print STDERR "But the CHECKSUM file can sometimes be wrong as the database may have been\n";
      print STDERR "updated without the CHECKSUM file being updated\n";
      print STDERR "To continue with just warning use the -force flag in the options\n";
      cleanup(1);
    }
  }
  
  system("gunzip -f $database/$file");

  my $str= "mysqlimport --fields_escaped_by=\\\\ $mysql_options ".$ENV{"PWD"}."/$database/$table$index.txt";
  print "$str\n" unless $quiet;
  $line = `$str`;
  if($line =~ /Error/ or $?){
    print STDERR $line;
    print STDERR "error code $?\n";
    print STDERR  "Error during mysqlimport\n";
    cleanup(1);
  }
  print $line unless $quiet;
  system("gzip $database/$table$index.txt");
  print "\n\n" unless $quiet;

}
close CHK;

cleanup();




sub cleanup{
  my $died = shift;
  if(defined($died) and $died){
    system("gzip $database/CHECKSUMS");
    exit 1;
  }
  if(defined($cleanup)){
    system("rm -Rf $database");
  }
  exit 0;
}


sub usage{
print << "EOH";
This perl script will download (rsync) the necesary ftp files and load them into a new local 
ensembl mysql database. It will check that the checksums match and do all the zipping and 
unzipping of the files.


 load_database_from_ftp.pl -root {root} -database {database} -new_database {database2} 
            -force -cleanup -quiet -help
            -host {host} -port {port} -user {user} -pass {password}
	    -mysqltempdir {dir}

  -root             Root directory for ftp files 

  -database         Database name to get data for

  -new_database     Name of the new database

  -user             User name to access database. Must allow writing.

  -pass             Password for user.

  -host             Database host.

  -port             Database port.

  -force            import data even if the checksums do not match

  -cleanup          remove the downloaded files at the end

  -quiet            No output except for serous error message

  -mysqltmpdir      Mysql may not have enough tmp space so this can be set to another directory

  -help             print this help text



examples:-

1) perl load_database_from_ftp_site.pl -database homo_sapiens_core_54_36p -host mysqlhostname
                                    -user mysqluser -pass mysqlpassword -force

This will download the ftp files for the 54 release of the human core database and create a database 
called <userid>_homo_sapiens_core_59_36p where userid is the login name of the user. To choose you 
own database name use the -new_database option.


2) load_database_from_ftp_site.pl -databases homo_sapiens_core_57_37d -new_database homo_sapiens_core_59_37d
       -host mysqlhostname -user mysqluser -pass mysqlpassword -quiet -cleanup -mysqltmpdir /scratch/

Will load the human core database into the mysql instance on mysqlhostname and use the directory 
/scratch/ to use as the tmp directory for mysql.

EOH

}	
