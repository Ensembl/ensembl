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

use strict;
use warnings;

use Getopt::Long;
use Bio::EnsEMBL::Registry;
my $reg = "Bio::EnsEMBL::Registry";

# download and import the database
my $root;
my $prefix="";
my $release;
my $specieslist;
my $grouplist;
my $user;
my $pass;
my $port=3306;
my $host;
my $cleanup = undef;
my $force = undef; # if set ignore checksum dies just wrtie warnings.
my $mysqltempdir = undef;
my $quiet = 0;
my $run = undef;

GetOptions ('root=s'                    => \$root,
            'prefix=s'                  => \$prefix,
	    'release=s'                 => \$release,
            'species=s'                 => \$specieslist,
	    'groups=s'                  => \$grouplist,
            'host=s'                    => \$host,
            'force'                     => \$force,
            'cleanup'                   => \$cleanup,
            'port=s'                    => \$port,
            'user=s'                    => \$user,
            'pass=s'                    => \$pass,
            'mysqltempdir=s'            => \$mysqltempdir,
            'quiet'                     => \$quiet,
            'run'                       => \$run,
            'help'  => sub { usage(); exit(0);}
           );



my @names;
if(defined($specieslist)){
  @names = split(",",$specieslist);
}
else{
  usage();
  die "No species set?\n";
}

my @types;
if(defined($grouplist)){
  @types = split(",",$grouplist);
}
else{
  usage();
  die "No groups set?\n";
}

my $db_version = undef;

#
#connect to latest databases to get species name
#

$reg->no_version_check(1);
$reg->load_registry_from_db(
			    -host => "ensembldb.ensembl.org",
			    -user => "anonymous",
			    -db_version => 59, # comment out later.
			   );

my @species;
foreach my $sp (@names){
  my $adap = $reg->get_adaptor($sp, "core", "slice");
  if(defined($adap)){
    my $name = $adap->dbc->dbname;
#    print $name."\n";
    if(defined($name)){
      if($name =~ /(\S+_\S+)_core/){
	push @species, $1;
#	print "sp is $1\n";
      }
    }
  }
  else{
   print  "Could not find species $sp so ignoring\n";
  }
}

if(defined($release)){
  if($release =~ /^\d+$/){
    $db_version = $release;
  }
  else{
    die "release must be an integer\n";
  }
}
else{
  $release = $reg->software_version();
}

my $sqlport = 5306;
if($release < 47){
  $sqlport = 3306;
}


my @database_list;

my $sqltemplate = 'mysql -hensembldb.ensembl.org -uanonymous -PPORT --skip-column-names -e\'show databases like "SPECIES%TYPE%RELEASE%"\'';

$sqltemplate =~ s/PORT/$sqlport/;
#print $sqltemplate."\n";

foreach my $sp (@species){
#  print $sp."\n";
  foreach my $ty (@types){
#    print "\t$ty\n";
    my $sql = $sqltemplate;
    $sql =~ s/SPECIES/$sp/;
    $sql =~ s/RELEASE/$release/;
    if($ty eq "all"){
      $sql =~ s/TYPE//;
    }
    else{
      $ty .= "\\_";
      $sql =~ s/TYPE/$ty/;
    }
#    print $sql."\n";
    my $line = `$sql`;
    my @vals = split(/\n/,$line);
    foreach my $db (@vals){
#      print "\t".$db."\n";
      push @database_list, $db;
    }
  }
}
if(!defined($host) or !defined $user){
  usage();
  die " No host or user\n";
}

#
# check mysql instance data to be copoed to.
#
my $com =  "mysql -h$host -u$user -P$port ";
if(defined($pass)){
  $com .= "-p$pass ";
}
$com .=   "-e'show databases like \"justatest\"' ";
#print $com."\n";

my $line = `$com`;
if($?){
  print $com." fails\n";
  die "$line";
}
if($line =~ /ERROR/){
  die "problem with mysql information\n$line\n";
}


use FindBin '$Bin';
my $com_init =  "perl ".$Bin."/load_database_from_ftp_site.pl -host $host -user $user ";
if(defined($force)){
  $com_init .= "-force ";
}
if(defined($cleanup)){
  $com_init .= "-cleanup ";
}
if(defined($pass)){
  $com_init .= "-pass $pass ";
}
if(defined($root)){
  $com_init .= "-root $root ";
}
if(defined($mysqltempdir)){
  $com_init .= "-mysqltempdir $mysqltempdir ";
}
if(defined($quiet)){
  $com_init .= "-quiet ";
}

my $okay="";
my $prob ="";

foreach my $db (@database_list){
  my $com =  "mysql -h$host -u$user -P$port ";
  if(defined($pass)){
    $com .= "-p$pass ";
  }
  $com .=   "-e'show databases like \"$prefix$db\"'";
  
#  print $db."\n";
  $line = `$com`;
#  print $line;
  if($line =~ /$db/ and !defined($force)){
    $prob .= "\t$prefix$db\n";
    next;
  }
  elsif(defined($run)){
    my $cmd = $com_init."-database $db -new_database $prefix$db ";
    
    print STDERR "Copying $db to $host as $prefix$db\n";
    my $output = `$cmd`;
    open(OUT,">$db.OUTPUT");
    print OUT $line;
    close OUT;
  }
  else{
    $okay .= "\t$db to $host $prefix$db\n";
  }

}



if(!defined($run)){
  if(length($prob) > 1){
    print "Problem with the following databases as they already exist on $host\n";
    print $prob;
  }
  if(length($okay) > 1){
    print "The following would be copied:-\n";
    print $okay;
  }
  print "\nYou need to set the flag -run to actually do the data copy\n";
  print "By default it is not done so that this list can be checked first\n";
}
else{
  if(length($prob) > 1){
    print "Problem with the following databases as they already exist on $host so not copied\n";
    print $prob;
  } 
}

sub usage{
print << "EOH";
It uses the Registry from the core API to get the species name to pass on to the script 
load_database_from_ftp.pl. 

 load_multiple_databases.pl -root {root} -prefix {prefix} -release {number}
            -species {s1,s2,s3} -groups {type1,type2} -force -cleanup -quiet -help
            -host {host} -port {port} -user {user} -pass {password}
	    -mysqltempdir {dir} -list

  -root             Root directory for ftp files 

  -prefix           Database name to get data for

  -release          Release version of the dtaabase to get

  -species          Comma separated list of species to get

  -groups           Comma separated list of database types to get 
                    ( from core,variation,funcgen,otherfeatures,vega etc or all)

  -user             User name to access database. Must allow writing.

  -pass             Password for user.

  -host             Database host.

  -port             Database port.

  -force            import data even if the checksums do not match
                    or the new database already exists.

  -cleanup          remove the downloaded files at the end

  -quiet            No output except for serous error message

  -mysqltmpdir      Mysql may not have enough tmp space so this can be set to another directory

  -run              If set will start the download etc else it will just list the databases. 
                    NOTE: Not default as this script does alot so we want to make sure everything 
                          is correct first before starting.

  -help             print this help text



examples:-

1) perl load_multiple_databases.pl -release 54  -groups core -species human -host mysqlhostname
                                    -user mysqluser -pass mysqlpassword -force -run -prefix "copy_"

This will download the ftp files for the 54 release of the human core database and create a database 
called copy_homo_sapiens_core_59_36p.


2) perl load_multiple_databases.pl -release 59 -species mouse -groups all -run
       -host mysqlhostname -user mysqluser -pass mysqlpassword -quiet -cleanup -mysqltmpdir /scratch/

Will load the mouse databases for release 59 into the mysql instance on mysqlhostname and use the directory 
/scratch/ to use as the tmp directory for mysql.
This will load the databases:- 
   mus_musculus_cdna_59_37l
   mus_musculus_core_59_37l
   mus_musculus_funcgen_59_37l
   mus_musculus_otherfeatures_59_37l
   mus_musculus_variation_59_37l
   mus_musculus_vega_59_37l

EOH

}	
