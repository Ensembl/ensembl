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

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Net qw/do_FTP/;
use Bio::EnsEMBL::Utils::Exception qw/throw/;
use Getopt::Long;

my ($db_name, $db_host, $db_user, $db_pass, $db_port, $help, $species, $group, $no_interactive);
my @user_karyotype;

my $NCBI_BASE = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All';

$db_port = 3306;
$species = 'human';
$group = 'core';

GetOptions(
  "db_name|dbname|database=s"      => \$db_name,
  "db_host|dbhost|host=s"          => \$db_host,
  "db_user|dbuser|user|username=s" => \$db_user,
  "db_pass|dbpass|pass|password=s" => \$db_pass,
  "db_port|dbport|port=s"          => \$db_port,
  "species=s"                      => \$species,
  'group=s'                        => \$group,
  'karyotype=s@'                   => \@user_karyotype,
  'no_interactive!'                => \$no_interactive,
  "h!"                             => \$help,
  "help!"                          => \$help,
);

if ($help) { &usage; exit 0; }
unless ($db_name and $db_host) { print "Insufficient arguments\n"; &usage; exit 1; }

if(@user_karyotype) {
  @user_karyotype = strings_to_array(@user_karyotype);
}

sub get_adaptor {
  my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -species => $species,
    -group   => $group,
    -dbname  => $db_name,
    -host    => $db_host,
    -user    => $db_user,
    -port    => $db_port
  );
  $dba->dbc->password($db_pass) if $db_pass;
  return $dba;
}

sub run {
  my $dba = get_adaptor();
  
  my @karyotype;
  if(@user_karyotype) {
    @karyotype = @user_karyotype;
  }
  else {
    @karyotype = @{ncbi_karyotype($dba)};
  }
  
  print_karyotype(@karyotype);
  
  my $ok_to_write = confirm('Is this correct?');
  if(!$ok_to_write) {
    my $new_karyotype = capture_user_input("Please give the correct karyotype (comma separated)");
    @karyotype = strings_to_array($new_karyotype);
    print_karyotype(@karyotype);
  }
  
  my $slices = karyotype_to_Slices($dba, @karyotype);
  my $has_karyotype = has_karyotype($slices);
  my $overwrite = 0;
  if($has_karyotype) {
    $overwrite = confirm('Karyotypes have already been assigned. Are you sure you want to overwrite');
    if(!$overwrite) {
      die "Cannot continue. Karyotypes have already been assigned to these Slices";
    }
  }
  
  write_ranks($dba, $slices);
  
  return;
}

sub accession {
  my ($dba) = @_;
  my $mc = $dba->get_MetaContainer();
  my $acc = $mc->single_value_by_key('assembly.accession', 1);
  if(!$acc) {
    print "Cannot continue. Species does not have a valid Genome Collections accession in its meta table\n";
    exit 1;
  }
  return $acc;
}

sub ncbi_karyotype {
  my ($dba) = @_;
  my $accession = accession($dba);
  print STDOUT "Fetching karyotype from NCBI using accession '$accession'\n";
  my $url = "$NCBI_BASE/${accession}.assembly.txt";
  my $content = do_FTP($url, 5, 2);
  my @lines = split(/\n/, $content);
  my @karyotype;
  foreach my $line (@lines) {
    next if $line =~ /^#/;
    chomp $line;
    my ($name, $role) = split(/\t/, $line);
    if($role eq 'assembled-molecule') {
      push(@karyotype, $name);
    }
  }
  return \@karyotype;
}

sub karyotype_to_Slices {
  my ($dba, @karyotype) = @_;
  my $sa = $dba->get_SliceAdaptor();
  my @slices;
  foreach my $name (@karyotype) {
    my $slice = $sa->fetch_by_region('toplevel', $name);
    push(@slices, $slice);
  }
  return \@slices;
}

sub confirm {
  my ($question) = @_;
  if($no_interactive) {
    return 1;
  }
  my $userinput = capture_user_input("$question [y(es)/N(o)]");
  return ($userinput =~ /^y(?:es)?$/xmsi) ? 1 : 0;
}

sub capture_user_input {
  my ($msg) = @_;
  die "Cannot continue; asking for user input but -no_interactive is on" if $no_interactive;
  print STDOUT "$msg: ";
  my $userinput =  <STDIN>;
  chomp ($userinput);
  return $userinput;
}

sub print_karyotype {
  my (@karyotype) = @_;
  print STDOUT "Karyotype to be used: ";
  print STDOUT join(q{,},@karyotype);
  print "\n";
  return;
}

sub strings_to_array {
  return split(q{,}, join(q{,}, @_));
}

sub has_karyotype {
  my ($slices) = @_;
  my $has_karyotype = 0;
  foreach my $slice (@{$slices}) {
    if($slice->has_karyotype()) {
      $has_karyotype = 1;
      last;
    }
  }
  return $has_karyotype;
}

sub write_ranks {
  my ($dba, $slices) = @_;
  my $aa = $dba->get_AttributeAdaptor();
  my $code = 'karyotype_rank';
  my $rank = 1;
  foreach my $slice (@{$slices}) {
    printf STDOUT '%s has been assigned rank %d', $slice->seq_region_name(), $rank;
    print "\n";
    $aa->remove_from_Slice($slice, $code);
    $aa->store_on_Slice($slice, $code, $rank);
    $rank++;
  }
  return;
}

sub usage {
    print "Description:

A program for writing karyotype ranks into a core-like database. If one is
not specified then we query NCBI's assembly report resource for a candidate
list. The program also allows you to specify a custom rank if required.

This module uses LWP to communicate with NCBI's FTP site. Should you have
any issues please consult its documentation on how to debug the issue. Also 
make sure you can access $NCBI_BASE

Synopsis:

  perl karyotype_rank_assigner.pl -db_name NAME -db_host HOST -db_user USR -db_pass PASS -species human -group GRP

Options:
    
    -db_name          The DB to add the rank to
    -database
    -dbname
    
    -db_host          Hostname for the DB
    -host
    -dbhost
    
    -db_user          Username for the DB
    -user
    -username
    -dbuser
    
    -db_pass          Password for the DB
    -pass
    -password
    -dbpass
    
    -db_port          Port for the DB
    -dbport
    -port
    
    -species          Name of the species; defaults to human
    -group            Name of the DB group; defaults to core
    
    -karyotype        Comma separated list of karyotypes to write. Also 
                      supports specifying this option multiple times
    
    -no_interactive   Do not prompt for confirmation
    
    -help             Print this message
";
}

run();
