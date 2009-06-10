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
my $partupdate;
my $notverbose ;
my $reset_to_mapping_finished;
my $reset_to_parsing_finished;
my $resubmit_failed_jobs;


my $options = join(" ",@ARGV);
print "Options: ".join(" ",@ARGV)."\n";

GetOptions ('file=s'                    => \$file,
            'dumpcheck'                 => \$dumpcheck, 
            'upload'                    => \$upload,
	    'notverbose'                => \$notverbose,  
            'nofarm'                    => \$nofarm, 
            'partupdate'                => \$partupdate,
            'reset_to_mapping_finished' => \$reset_to_mapping_finished,
            'reset_to_parsing_finished' => \$reset_to_parsing_finished,
            'resubmit_failed_jobs'      => \$resubmit_failed_jobs,
            'help'                      => sub { usage(); exit(0); } );


my $mapper = XrefMapper::BasicMapper->process_file($file, !$notverbose);

$mapper->add_meta_pair("mapper options",$options);


if(defined($dumpcheck)){
  $mapper->dumpcheck("yes");
}
if(defined($nofarm)){
  $mapper->nofarm("yes");
}
if(!defined($notverbose)){
    print "running in verbose mode\n";
  $mapper->verbose(1);
}
else{
  print "running in quite mode.\n";
  $mapper->verbose(0);
}

# find out what stage the database is in at present.
my $status = $mapper->xref_latest_status(0);
print "current status is $status\n" if ($mapper->verbose);



if(defined($reset_to_mapping_finished)){
  $mapper->revert_to_mapping_finished();
  exit();
}

if(defined($reset_to_parsing_finished)){
  $mapper->revert_to_parsing_finished();
  exit();
}


my $submitter = XrefMapper::SubmitMapper->new($mapper);

if(defined($resubmit_failed_jobs)){
  print STDERR "Resubmitting failed jobs\n";
  $submitter->fix_mappings();
}

if( $status eq "parsing_finished" or $status eq "xref_fasta_dumped"){ 
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
  die "Status is $status so a job is already doing the mapping. Please wait till it has finished\n If this is not the case then you will have to manually change the status and try again\n";
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

my $fullmode;
if(defined($partupdate)){
  print "partupdate is set to $partupdate over ruling the precalulated value\n" if($mapper->verbose);
  $fullmode = 0;
}
else{
  if($mapper->get_meta_value("fullmode") eq "yes"){
    $fullmode = 1;
  }
  elsif($mapper->get_meta_value("fullmode") eq "no"){
    $fullmode = 0;
  }
  else{
    print "WARNING: No value for fullmode in meta table using fullmode anyway\n";
    $fullmode = 1;
  }
}


if(($status eq "core_loaded" or $status eq "display_xref_done") and $upload){

  
  my $display = XrefMapper::DisplayXrefs->new($mapper);
  $display->genes_and_transcripts_attributes_set($fullmode);
  
}
#}
#else{
#  my $display = XrefMapper::DisplayXrefs->new($mapper);
#  $display->pump_up_the_jam();
#  #$display->while_your_feet_are_stomping();
#  #$display->set_status(); # set KNOWN,NOVEL etc 

#}



sub usage {

  print << "EOF";
  xref_mapper.pl -file {config_file} -verbose -upload -nofarm 

  -file             Name of configuration file.
                    For more info on this file see end of help.

  -dumpcheck        Check if the fasta files already exist.
                    If so do not redump them.

  -upload           Write data from xref database to core database.

  -notverbose       Do not give information about progress and possible warnings.

  -nofarm           Run the exonerate jobs locally and not on the compute farm.

  -partupdate       Not all xrefs have been updated hence to get the gene descriptions
                    etc we need to work out these via the core database.
                    (Note this is much slower, but has to be done if you are only updating 
                     a few xref sources) 
                    By default this is calulated from the parsing options used.
                    ONLY set if you know what the consequences will be!!

  -reset_to_mapping_finished
                    Reset the status of the database (cleaning up the tables) to be equivalent
                    of the mapping having been finished and the map files are ready to be processed.

  -reset_to_parsing_finished
                    Reset the status of the database (cleaning up the tables) to be equivalent
                    of the parsing having just been done.

  -resubmit_failed_jobs
                    This will find the failed mapping jobs and rerun these and then contiune as normal.
 
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




























