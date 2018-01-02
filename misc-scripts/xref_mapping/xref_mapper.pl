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

use Getopt::Long qw(:config pass_through);
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
use XrefMapper::DisplayXrefs;
use XrefMapper::CoordinateMapper;
use XrefMapper::UniParcMapper;
use XrefMapper::RNACentralMapper;
use XrefMapper::OfficialNaming;
use XrefMapper::DirectXrefs;

use vars qw(@INC);

$| = 1;

my $file;
my $dumpcheck;
my $upload = 0;
my $nofarm;
my $notverbose ;
my $reset_to_mapping_finished;
my $reset_to_parsing_finished;
my $resubmit_failed_jobs;
my $recalc_display_xrefs;
my $no_coord_mapping = 0;

my $options = join(" ",@ARGV);
print "Options: ".join(" ",@ARGV)."\n";

my $ret = Getopt::Long::GetOptions ('file=s'          => \$file,
            'dumpcheck'                 => \$dumpcheck, 
            'upload'                    => \$upload,
	    'notverbose'                => \$notverbose,  
            'nofarm'                    => \$nofarm, 
            'reset_to_mapping_finished' => \$reset_to_mapping_finished,
            'reset_to_parsing_finished' => \$reset_to_parsing_finished,
            'resubmit_failed_jobs'      => \$resubmit_failed_jobs,
            'recalc_display_xrefs'    => \$recalc_display_xrefs,
            'no_coord_mapping'          => \$no_coord_mapping,
            'help'                      => sub { usage(); exit(0); } );


if($ARGV[0]){
  print STDERR "Unknown command line arguments:-\n";
  foreach my $a (@ARGV){
    print STDERR "\t".$a."\n";
  }
  print STDERR "Stopping script. Please fix the command line.\n";
  print STDERR "use -help for full list of command line options.\n";;
  exit(1);
}


if(!defined($file)){
  usage();
  exit;
}

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


if($recalc_display_xrefs){
  my $display = XrefMapper::DisplayXrefs->new($mapper);
  $display->remove_source_priorities(); 
  $display->genes_and_transcripts_attributes_set();
  exit();
}


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

if( $status eq "parsing_finished"){
  $mapper->get_alt_alleles();
}

$status = $mapper->xref_latest_status();
if( $status eq "alt_alleles_added" or $status eq "xref_fasta_dumped"){ 
  print "\nDumping xref & Ensembl sequences\n"  if ($mapper->verbose);
  $submitter->dump_seqs();
}
else{
  $submitter->no_dump_xref()
}

$status = $mapper->xref_latest_status();
if($status eq "core_fasta_dumped"){# load core data needed
  my $core_info = XrefMapper::CoreInfo->new($mapper);
  $core_info->get_core_data();
}

$status = $mapper->xref_latest_status();
if($status eq "core_data_loaded") {
  $submitter->build_list_and_map();
  $status =  $mapper->xref_latest_status();
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
if($status eq "mapping_processed"){ # add direct xrefs to the object_xref etctables
  my $direct_mappings = XrefMapper::DirectXrefs->new($mapper);
  $direct_mappings->process();
}

$status = $mapper->xref_latest_status();
if($status eq "direct_xrefs_parsed"){  # process the priority xrefs.
  my $pp = XrefMapper::ProcessPrioritys->new($mapper);
  $pp->process();
}


# pair data
$status = $mapper->xref_latest_status();
if($status eq "prioritys_flagged"){  # process the inferred pairs
  my $pp = XrefMapper::ProcessPaired->new($mapper);
  $pp->process();

}


# Biomart test each external_db on one ensembl type only and fix it
$status = $mapper->xref_latest_status();
if($status eq "processed_pairs"){ 
  $mapper->biomart_testing();
}


$status = $mapper->xref_latest_status();
if($status eq "biomart_test_finished"){ 
  $mapper->source_defined_move();
}

$status = $mapper->xref_latest_status();
if($status eq "source_level_move_finished"){ 
    $mapper->process_alt_alleles();
}



# species specific processing
# i.e. for mouse and human add MGI_curated_gene and HGNC_curated_gene

$status = $mapper->xref_latest_status();
if($status eq "alt_alleles_processed"){ 
  my $official_naming = XrefMapper::OfficialNaming->new($mapper);
  $official_naming->run();
}


#UniParc and RNACentral checksum mapping
#Allow for reruns if required and let the mapper decide if it can run or not
$status = $mapper->xref_latest_status();
if($status eq 'official_naming_done') {
  my $checksum_mapper = XrefMapper::UniParcMapper->new($mapper);
  $checksum_mapper->process();
  $checksum_mapper = XrefMapper::RNACentralMapper->new($mapper);
  $checksum_mapper->process();
}

# Coordinate xrefs

# tests
$status = $mapper->xref_latest_status();
if($status eq "official_naming_done"  || $status eq "checksum_xrefs_finished" || $status eq "tests_started" || $status eq "tests_failed" ){
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
if($status eq "tests_finished" and $upload and !$no_coord_mapping){
  my $coord = XrefMapper::CoordinateMapper->new($mapper);
  $coord->run_coordinatemapping($upload);
}

$status = $mapper->xref_latest_status();
if(($status eq "coordinate_xref_finished" and $upload) or ($no_coord_mapping and $status eq "tests_finished" and $upload) ){
  my $loader = XrefMapper::XrefLoader->new($mapper);
  $loader->update();
  
}

# generate display_xrefs and descriptions for gene and transcripts. 

$status = $mapper->xref_latest_status();

#my $fullmode;
#if(defined($partupdate)){
#  print "partupdate is set to $partupdate over ruling the precalulated value\n" if($mapper->verbose);
#  $fullmode = 0;
#}
#else{
#  if($mapper->get_meta_value("fullmode") eq "yes"){
#    $fullmode = 1;
#  }
#  elsif($mapper->get_meta_value("fullmode") eq "no"){
#    $fullmode = 0;
#  }
#  else{
#    print "WARNING: No value for fullmode in meta table using fullmode anyway\n";
#    $fullmode = 1;
#  }
#}


if(($status eq "core_loaded" or $status eq "display_xref_done") and $upload){

  
  my $display = XrefMapper::DisplayXrefs->new($mapper);
  $display->genes_and_transcripts_attributes_set();
  
}

if(!defined($notverbose)){
  print "xref_mapper.pl FINISHED NORMALLY\n";
}

sub usage {

  print << "EOF";
  xref_mapper.pl -file {config_file} -upload -nofarm 

  -file             Name of configuration file.
                    For more info on this file see end of help.

  -dumpcheck        Check if the fasta files already exist.
                    If so do not redump them.

  -upload           Write data from xref database to core database.

  -notverbose       Do not give information about progress and possible warnings.

  -nofarm           Run the exonerate jobs locally and not on the compute farm.

  -reset_to_mapping_finished
                    Reset the status of the database (cleaning up the tables) to be 
                    equivalent of the mapping having been finished and the map files 
                    are ready to be processed.

  -reset_to_parsing_finished
                    Reset the status of the database (cleaning up the tables) to be 
                    equivalent of the parsing having just been done.

  -resubmit_failed_jobs
                    This will find the failed mapping jobs and rerun these and then 
                    contiune as normal.

  -recalc_display_xrefs
                    Only recalculate the display xrefs and genes descriptions. 
                    

                  
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

farm
queue=long
exonerate=/software/ensembl/bin/exonerate-1.4.0
####################################################
host1 can be the same as host2.
user1 can be the same as user2 but user2 must have write access and user1 must have at 
least read access. The directorys set by dir= must already exist.



EOF

}




























