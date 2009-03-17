#!/usr/bin/perl


use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Cwd;
use XrefMapper::db;
use XrefMapper::SubmitMapper;
use XrefMapper::ProcessMappings;
use XrefMapper::ProcessPrioritys;
use XrefMapper::ProcessPaired;
use XrefMapper::CoreInfo;
use XrefMapper::TestMappings;
use XrefMapper::XrefLoader;
use XrefMapper::Interpro;
use XrefMapper::DisplayXrefs;
use XrefMapper::CoordinateMapper;

use vars qw(@INC);

$| = 1;

my $file;
my $dumpcheck;
my $upload = 0;
my $nofarm;
my $verbose;

print "Options: ".join(" ",@ARGV)."\n";

GetOptions ('file=s'                    => \$file,
            'dumpcheck'                 => \$dumpcheck, 
            'upload'                    => \$upload,
	    'verbose'                   => \$verbose,  
            'nofarm'                    => \$nofarm, 
            'help'                      => sub { usage(); exit(0); } );


my $mapper = XrefMapper::BasicMapper->process_file($file, $verbose);


if(defined($dumpcheck)){
  $mapper->dumpcheck("yes");
}
if(defined($nofarm)){
  $mapper->nofarm("yes");
}
if(defined($verbose)){
  $mapper->verbose(1);
}
else{
  $mapper->verbose(0);
}

# find out what stage the database is in at present.
my $status = $mapper->xref_latest_status($mapper->verbose);
print "current status is $status\n" if ($mapper->verbose);


my $submitter = XrefMapper::SubmitMapper->new($mapper);

if( $status eq "parsing_finished" or $status = "xref_fasta_dumped"){ 
  print "\nDumping xref & Ensembl sequences\n"  if ($mapper->verbose);
  $submitter->dump_seqs();
}
else{
  $submitter->no_dump_xref()
}



$status = $mapper->xref_latest_status();
if($status eq "core_fasta_dumped"){
  $submitter->build_list_and_map();
  $status =  $mapper->xref_latest_status();
}
else{

}

$status = $mapper->xref_latest_status();
if($status eq "mapping_started"){
  die "Status is $status so a job is already doing the mapping. Please wait till it has finished";
}
elsif($status eq "mapping_finished"){
  my $parser = XrefMapper::ProcessMappings->new($mapper);
  $parser->process_mappings();
}

$status = $mapper->xref_latest_status();
if($status eq "mapping_processed"){  # load core data needed and add direct xrefs to the object_xref etc tables
  my $core_info = XrefMapper::CoreInfo->new($mapper);
  $core_info->get_core_data();
}


$status = $mapper->xref_latest_status();
if($status eq "direct_xrefs_parsed"){  # process the priority xrefs.
  my $pp = XrefMapper::ProcessPrioritys->new($mapper);
  $pp->process();
}


# pair data
$status = $mapper->xref_latest_status();
if($status eq "prioritys_flagged"){  # process the inferred pairs and interpro xrefs.
  my $pp = XrefMapper::ProcessPaired->new($mapper);
  $pp->process();

  my $inter = XrefMapper::Interpro->new($mapper);
  $inter->process();
}


# Biomart test each external_db on one ensembl type only and fix it
$status = $mapper->xref_latest_status();
if($status eq "processed_pairs"){ 
  $mapper->biomart_testing();
}

# species specific processing
# i.e. for mouse and human add MGI_curated_gene and HGNC_curated_gene

$status = $mapper->xref_latest_status();
if($status eq "processed_pairs"){ 
  $mapper->official_naming();
}


# Coordinate xrefs

# tests
$status = $mapper->xref_latest_status();
if($status eq "official_naming_done" || $status eq "tests_started" || $status eq "tests_failed" ){  # process the priority xrefs.
  my $tester = XrefMapper::TestMappings->new($mapper);
  if($tester->unlinked_entries){
    die "Problems found so will not load core database\n";
  }
  $tester->direct_stable_id_check(); # Will only give warnings
  $tester->entry_number_check();     # Will only give warnings
  $tester->name_change_check();      # Will only give warnings
}

# load into core database
$status = $mapper->xref_latest_status();
if($status eq "tests_finished" and $upload){

  my $coord = XrefMapper::CoordinateMapper->new($mapper);
  $coord->run_coordinatemapping($upload);
}

$status = $mapper->xref_latest_status();
if($status eq "coordinate_xref_finished" and $upload){

  my $loader = XrefMapper::XrefLoader->new($mapper);
  $loader->update();
  
}

# generate display_xrefs and descriptions for gene and transcripts. 

$status = $mapper->xref_latest_status();
if(($status eq "core_loaded" or $status eq "display_xref_done") and $upload){

  my $display = XrefMapper::DisplayXrefs->new($mapper);
  $display->genes_and_transcripts_attributes_set();

}

sub usage {

  print << "EOF";
  xref_mapper.pl -file {config_file} -verbose -upload -nofarm 

  -file             Name of configuration file.
                    For more info on this file see end of help.

  -dumpcheck        Check if the fasta files already exist.
                    If so do not redump them.

  -upload           Write data from xref database to core database.

  -verbose          Give information about progress and possible warnings.
                    (Very much recomended)

  -nofarm           Run the exonerate jobs locally and not on the compute farm.


Below is an example of the configuration file
####################################################
xref
host=host1
port=3306
dbname=hedgehog_xref_54
user=user1
password=pass1
dir=./xref

species=erinaceus_europaeus
host=host2
port=3306
dbname=erinaceus_europaeus_core_54_1e
user=user2
password=pass2
dir=./ensembl
####################################################
host1 can be the same as host2.
user1 can be the same as user2 but user2 must have write access and user1 must have at least read access.
The directorys set by dir= must already exist.



EOF

}




























