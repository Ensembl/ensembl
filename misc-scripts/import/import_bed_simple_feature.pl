#!/usr/bin/env perl

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

use Getopt::Long;

my ($file,$db_name,$db_host,$db_user,$db_pass,$db_port,$db_version,$help,$species,$group);
my ($logic_name, $description, $display_label);
my ($directory);
$species = "human";
$group = 'core';

GetOptions ("file=s" => \$file,
            "db_name|dbname|database=s" => \$db_name,
            "db_host|dbhost|host=s" => \$db_host,
            "db_user|dbuser|user|username=s" => \$db_user,
            "db_pass|dbpass|pass|password=s" => \$db_pass,
            "db_port|dbport|port=s" => \$db_port,
            "species=s" => \$species,
            'group=s'   => \$group,
            'logic_name=s' => \$logic_name,
            'description=s' => \$description,
            'display_label=s' => \$display_label,
            'directory=s' => \$directory,
            "h!"        => \$help,
            "help!"     => \$help,
);

if ($help) {&usage; exit 0;}
unless ($file and $db_name and $db_host) {print "Insufficient arguments\n"; &usage; exit 1;}
unless ($logic_name) { print "No logic name given\n"; usage(); exit 1; } 

my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  -species => $species,
  -group => $group,
  -dbname => $db_name,
  -host => $db_host,
  -user => $db_user,
  -pass => $db_pass,
  -port => $db_port
);

run();

sub run {
  if(! -f $file) {
    die "No file found at $file";
  }
  process_file($file);
  return;
}

sub process_file {
  my ($f) = @_;
  my $analysis = get_Analysis();
  my @features;
  my $count = 0;
  iterate_file($f, sub {
    my ($line) = @_;
    if($count != 0 && $count % 500 == 0) {
      printf STDERR "Processed %s records\n", $count;
    }
    chomp $line;
    my $sf = line_to_SimpleFeature($line, $analysis);
    push(@features, $sf);
    $count++;
  });
  my $sfa = $dba->get_SimpleFeatureAdaptor();
  print STDERR "Storing\n";
  $sfa->store(@features);
  print STDERR "Done\n";
  return;
}

sub line_to_SimpleFeature {
  my ($line, $analysis) = @_;
  my ($chr, $start, $end, $label, $score, $strand) = split(/\t/, $line);
  $start++; # UCSC is 0 idx start
  $score ||= 0;
  $strand ||= 1;
  my $slice = get_Slice($chr);
  my $sf = Bio::EnsEMBL::SimpleFeature->new(
      -start => $start,
      -end => $end,
      -score => $score,
      -analysis => $analysis,
      -slice => $slice,
      -strand => $strand,
      -display_label => $label,
  );
  return $sf;
}

my %slices;
sub get_Slice {
  my ($original) = @_;
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

sub usage {
    print "Launching instructions:
    Run from a folder you are happy to have filled with files.

Description:

Import data from a BED file into the simple_feature table. Only supports
6 column BED files (location, name and score).

Synopsis:

  perl import_bed_simple_feature.pl -file [PATH} -db_name NAME

Options:
    
    -file           Supply the file path
    -directory      Location to download to; defaults to current directory
    -logic_name     Analysis logic name import data against
    
    -db_name        The DB to add these features to
    -database
    -dbname
    
    -db_host        Hostname for the DB
    -host
    -dbhost
    
    -db_user        Username for the DB
    -user
    -username
    -dbuser
    
    -db_pass        Password for the DB
    -pass
    -password
    -dbpass
    
    -db_port        Port for the DB
    -dbport
    -port
    
    -species        Name of the species; defaults to human
    -group          Name of the DB group; defaults to core
    -description    Analysis description; only needed if analysis is not already in the DB
    -display_label  Analysis display label for the website; only needed if analysis is not already in the DB
    -help
";    
}