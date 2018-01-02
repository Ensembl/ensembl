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


# Designed to work on BED data such as that available from UCSC or 1KG:
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/accessible_genome_masks
#

use strict;
use warnings;

use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor;
use Bio::EnsEMBL::Utils::IO qw/iterate_file/;
use POSIX qw/strftime/;
use Getopt::Long;

my ($file,$db_name,$db_host,$db_user,$db_pass,$db_port,$help,$species,$group);
my ($dna_db_name,$dna_db_host,$dna_db_user,$dna_db_pass,$dna_db_port, $dna_group);
my ($logic_name, $description, $display_label);
my $write_every = -1;
$species = "human";
$group = 'core';

GetOptions ("file=s" => \$file,
            "db_name|dbname|database=s" => \$db_name,
            "db_host|dbhost|host=s" => \$db_host,
            "db_user|dbuser|user|username=s" => \$db_user,
            "db_pass|dbpass|pass|password=s" => \$db_pass,
            "db_port|dbport|port=s" => \$db_port,
            "dna_db_name|dna_dbname|dna_database=s" => \$dna_db_name,
            "dna_db_host|dna_dbhost|dna_host=s" => \$dna_db_host,
            "dna_db_user|dna_dbuser|dna_user|dna_username=s" => \$dna_db_user,
            "dna_db_pass|dna_dbpass|dna_pass|dna_password=s" => \$dna_db_pass,
            "dna_db_port|dna_dbport|dna_port=s" => \$dna_db_port,
            "dna_group=s" => \$dna_group,
            "species=s" => \$species,
            'group=s'   => \$group,
            'logic_name=s' => \$logic_name,
            'description=s' => \$description,
            'display_label=s' => \$display_label,
            'write_every=i' => \$write_every,
            "h!"        => \$help,
            "help!"     => \$help,
);

if ($help) {&usage; exit 0;}
unless ($file and $db_name and $db_host) {print "Insufficient arguments\n"; &usage; exit 1;}
unless ($logic_name) { print "No logic name given\n"; usage(); exit 1; } 

sub get_adaptor {
  my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -species => $species,
    -group => $group,
    -dbname => $db_name,
    -host => $db_host,
    -user => $db_user,
    -port => $db_port
  );
  $dba->dbc->password($db_pass) if $db_pass;
  
  if($dna_db_name) {
    $dna_group ||= 'core';
    $dna_db_host ||= $db_host;
    $dna_db_port ||= $db_port;
    $dna_db_user ||= $db_user;
    my $dna_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
      -species => $species,
      -group => $dna_group,
      -dbname => $dna_db_name,
      -host => $dna_db_host,
      -user => $dna_db_user,
      -port => $dna_db_port
    );
    $dna_dba->dbc->password($dna_db_pass) if $dna_db_pass;
    $dba->dnadb($dna_dba);
  }
  return $dba;
}

# see bottom of file for this method call
sub run {
  my $dba = get_adaptor();
  if(! -f $file) {
    die "No file found at $file";
  }
  process_file($file, $dba);
  return;
}

sub process_file {
  my ($f, $dba) = @_;
  my $analysis = get_Analysis($dba);
  my @features;
  my $count = 0;
  my $commit_count = 0;
  iterate_file($f, sub {
    my ($line) = @_;
    if($count != 0 && $count % 2000 == 0) {
      info("Processed %s records", $count);
    }
    chomp $line;
    my $sf = line_to_SimpleFeature($line, $analysis, $dba);
    push(@features, $sf);
    $count++;
    $commit_count++;
    
    if($commit_count == $write_every) {
      _store(\@features, $dba);
      @features = ();
      $commit_count = 0;
    }
  });
  _store(\@features, $dba);
  return;
}

sub _store {
  my ($features, $dba) = @_;
  my $sfa = $dba->get_SimpleFeatureAdaptor();
  my $count = scalar(@{$features});
  if($count > 0) {
    info("Writing %d feature(s)", $count);
    $sfa->store(@{$features});
    info("Done");
  }
  return;
}

sub line_to_SimpleFeature {
  my ($line, $analysis, $dba) = @_;
  my ($chr, $start, $end, $label, $score, $ucsc_strand) = split(/\t/, $line);
  $start++; # UCSC is 0 idx start
  #if was defined & -ve then set as so. +ve is default
  my $strand = (defined $ucsc_strand && $strand eq '-') ? -1 : 1;
  my $slice = get_Slice($chr, $dba);
  my %args = (
    -start => $start,
    -end => $end,
    -analysis => $analysis,
    -slice => $slice,
    -strand => $strand,
    -display_label => $label,
  );
  $args{-SCORE} = $score if defined $score; # only add score if it was there
  my $sf = Bio::EnsEMBL::SimpleFeature->new(%args);
  return $sf;
}

my %slices;
sub get_Slice {
  my ($original, $dba) = @_;
  my $name = $original;
  $name =~ s/^chr//;
  return $slices{$name} if exists $slices{name};
  my $slice = $dba->get_SliceAdaptor()->fetch_by_region('toplevel', $name);
  if(!$slice) {
    die "Could not get a Slice from the Ensembl database for the given region '$original' or '$name' and coorindate system 'toplevel'. Check your core database";
  }
  $slices{$name} = $slice;
  return $slice;
}

sub get_Analysis {
  my ($dba) = @_;
  my $aa = $dba->get_AnalysisAdaptor();
  my $analysis = $aa->fetch_by_logic_name($logic_name);
  if(!$analysis) {
    $analysis = Bio::EnsEMBL::Analysis->new(
      -displayable => 1,
      -logic_name => $logic_name,
    );
    $analysis->description($description) if $description;
    $analysis->display_label($display_label) if $display_label;
    $aa->store($analysis);
  }
  return $analysis;
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
    print "Launching instructions:
    Run from a folder you are happy to have filled with files.

Description:

Import data from a BED file into the simple_feature table. Only supports
6 column BED files (location, name and score).

Synopsis:

  perl import_bed_simple_feature.pl -file [PATH} -db_name NAME

Options:
    
    -file               Supply the file path
    -logic_name         Analysis logic name import data against
    
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
    
    -dna_db_name        The DNA DB to use if DB does not contain coordinate systems and DNA
    -dna_database
    -dna_dbname
    
    -dna_db_host        Hostname for the DNA DB. Defaults to -host
    -dna_host
    -dna_dbhost
    
    -dna_db_user        Username for the DNA DB. Defaults to -user
    -dna_user
    -dna_username
    -dna_dbuser
    
    -dna_db_pass        Password for the DNA DB. Defaults to -pass
    -dna_pass
    -dna_password
    -dna_dbpass
    
    -dna_db_port        Port for the DNA DB. Defaults to -port
    -dna_dbport
    -dna_port
    
    -dna_group          
    
    -species        Name of the species; defaults to human
    -group          Name of the DB group; defaults to core
    -description    Analysis description; only needed if analysis is not already in the DB
    -display_label  Analysis display label for the website; only needed if analysis is not already in the DB
    
    -write_every    Write features once every N lines. Defaults to -1 (write once all records are parsed)
    -help
";
}

run();
