#!/usr/local/ensembl/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor;
use Bio::EnsEMBL::AssemblyExceptionFeature;

$| = 1;

my $host         = '';
my $user         = '';
my $pass         = '';
my $port         = '';
my $dbname       = '';
my @patch_types = ('PATCH_FIX','PATCH_NOVEL');
my $align_by_component;
my $align_non_ident;
my $fix_overlap;
my $new_align_session;
my $check_repeat_masked;
my $out_dir;
my $perl;
my $pipe_path;
my $script;
my $assembly;

&GetOptions( 'dbhost:s'       => \$host,
             'dbuser:s'       => \$user,
             'dbname:s'       => \$dbname,
             'dbpass:s'       => \$pass,
             'dbport:n'       => \$port,
             'align_by_component'  => \$align_by_component,
             'align_non_ident'     => \$align_non_ident,
             'new_align_session'   => \$new_align_session,
             'check_repeat_masked' => \$check_repeat_masked,
             'fix_overlap'         => \$fix_overlap,
             'out_dir:s'           => \$out_dir,
             'perl:s'              => \$perl,
             'pipeline_path:s'     => \$pipe_path,
             'assembly:s'          => \$assembly);

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $host,
                                             -user   => $user,
                                             -pass   => $pass,
                                             -port   => $port,
                                             -dbname => $dbname );

if((!($align_by_component or $align_non_ident or $fix_overlap or $new_align_session or $check_repeat_masked))){
  throw("Must specify one of -align_by_component, -align_non_ident, -new_align_session, -check_repeat_masked or -fix_overlap.\n");
}
if((!$out_dir) or (!(-d $out_dir))){
  throw("Must specify a valid output directory.\n");
}
if(!$perl){
  throw("Must specify perl path.\n");
}
if(!$pipe_path){
  throw("Must specify path to ensembl-pipeline.\n");
}
if(!$assembly){
  throw("Must specify assembly i.e. GRCh37.\n");
}
if($align_by_component){
  $script = $pipe_path."/scripts/Finished/assembly/align_by_component_identity.pl";
  if(! -e $script){
    throw("Can't find script: ".$script."\n");
  }
}
elsif($align_non_ident){
  $script = $pipe_path."/scripts/Finished/assembly/align_nonident.pl";
  if(! -e $script){
    throw("Can't find script: ".$script."\n");
  }
}
elsif($new_align_session){
  $script = $pipe_path."/scripts/Finished/assembly/new_align_session.pl";
  if(! -e $script){
    throw("Can't find script: ".$script."\n");
  }
}
elsif($check_repeat_masked){
  $script = $pipe_path."/scripts/Finished/assembly/check_repeat_masked.pl";
  if(! -e $script){
    throw("Can't find script: ".$script."\n");
  }
}
elsif($fix_overlap){
  $script = $pipe_path."/scripts/Finished/assembly/fix_overlaps.pl";
  if(! -e $script){
    throw("Can't find script: ".$script."\n");
  }
}

#get the patches
print "Getting patches...\n";
my $asm_exc_adaptor = $db->get_AssemblyExceptionFeatureAdaptor();
my @exceptions = @{$asm_exc_adaptor->fetch_all()};
my @patches;
EXC: foreach my $exc (@exceptions){
  foreach my $type (@patch_types){
    if($exc->type() =~ m/$type/){
      push(@patches, $exc);
      next EXC;
    }
  }
}
#Assuming that AssemblyExceptionFeatureAdaptor's fetch_all will always return two 
#entries for each patch and that those two entries are adjacent
my $num_patches = scalar(@patches)/2;
print "Have ".$num_patches." patches.\n";

my $command = "";
my @ref_chroms;
my @patch_chroms;

#for each patch
for (my $i = 0;$i < $num_patches;$i++){

  #get the two slices
  my $ref_slice;
  my $patch_slice;
  for(my $j = 0;$j < 2; $j++){
    my $exc = pop(@patches);
    #if this is the ref version
    if($exc->type =~ m/REF/){
      #alt is only the patch slice
      $patch_slice = $exc->alternate_slice();
    }
    else{
      #alt is replaced region of ref
      $ref_slice = $exc->alternate_slice();
    }    
  }
  if(!($patch_slice and $ref_slice)){
    throw("Something is wrong, the patch and ref slices were not set correctly.\n");
  }

  print "Reached ".$patch_slice->name()." and ".$ref_slice->name()."\n";

  #running for all patches together below
  if($align_by_component or $fix_overlap or $new_align_session or $check_repeat_masked){
    push @ref_chroms, $ref_slice->seq_region_name;
    push @patch_chroms, $patch_slice->seq_region_name;
  }
  elsif($align_non_ident){
    #run for each patch
    $command = "bsub -R'select[mem>4000] rusage[mem=4000]' -M4000000 -o ".$out_dir.$patch_slice->seq_region_name."_nonident.out -e ".
      $out_dir.$patch_slice->seq_region_name."_nonident.err ".
      $perl." ".$script." --dbhost ".$host." --dbport ".$port." --dbuser ".$user." --dbpass ".$pass." --dbname ".$dbname.
      " --assembly ".$assembly." --altdbname ".$dbname." --altassembly ".$assembly." --chromosomes ".$ref_slice->seq_region_name.
      " --altchromosomes ".$patch_slice->seq_region_name." --logpath ".$out_dir." --logfile ".$patch_slice->seq_region_name."_nonident.log --force-stage";

    if ( ($patch_slice->seq_region_name eq 'HG1472_PATCH') or ($patch_slice->seq_region_name eq 'HG858_PATCH') ) {
      $command = "bsub -R'select[mem>10000] rusage[mem=10000]' -M10000000 -o ".$out_dir.$patch_slice->seq_region_name."_nonident.out -e ".
        $out_dir.$patch_slice->seq_region_name."_nonident.err ".
        $perl." ".$script." --dbhost ".$host." --dbport ".$port." --dbuser ".$user." --dbpass ".$pass." --dbname ".$dbname.
        " --assembly ".$assembly." --altdbname ".$dbname." --altassembly ".$assembly." --chromosomes ".$ref_slice->seq_region_name.
        " --altchromosomes ".$patch_slice->seq_region_name." --logpath ".$out_dir." --logfile ".$patch_slice->seq_region_name."_nonident.log --force-stage";
    }

    print $command."\n";
    my $cmd_rtn = `$command`;
    print $cmd_rtn."\n";
  }

}

if($align_by_component){

  my $ref_chromosomes = join ",", @ref_chroms;
  my $patch_chromosomes = join ",", @patch_chroms;

  $command = "bsub -R\"select[mem>3800] rusage[mem=3800]\" -M3800000 -o ".$out_dir."component_ident.out -e ".$out_dir."component_ident.err ".
    $perl." ".$script." --host ".$host." --port ".$port." --user ".$user." --pass ".$pass." --dbname ".$dbname." --altdbname ".$dbname.
    " --assembly ".$assembly." --altassembly ".$assembly." --chromosomes ".$ref_chromosomes.
    " --altchromosomes ".$patch_chromosomes." --logpath ".$out_dir." --logfile component_ident.log";

  print $command."\n";

  my $cmd_rtn = `$command`;
  print $cmd_rtn."\n";
}
elsif($new_align_session){
  my $ref_chromosomes = join ",", @ref_chroms;
  my $patch_chromosomes = join ",", @patch_chroms;

  $command = "bsub -o ".$out_dir."new_align_session.out -e ".$out_dir."new_align_session.err ".
    $perl." ".$script." --host ".$host." --port ".$port." --user ".$user." --pass ".$pass." --dbname ".$dbname." --altdbname ".$dbname.
    " --assembly ".$assembly." --altassembly ".$assembly." --chromosomes ".$ref_chromosomes.
    " --altchromosomes ".$patch_chromosomes." --logpath ".$out_dir." --logfile new_align_session.log -author patch -nolog";

  print $command."\n";

  my $cmd_rtn = `$command`;
  print $cmd_rtn."\n";

}
elsif($fix_overlap){

  my $ref_chromosomes = join ",", @ref_chroms;
  my $patch_chromosomes = join ",", @patch_chroms;

  $command = "bsub -o ".$out_dir."fix_overlap.out -e ".$out_dir."fix_overlap.err ".
    $perl." ".$script." --dbhost ".$host." --dbport ".$port." --dbuser ".$user." --dbname ".$dbname." --altdbname ".$dbname.
    " --dbpass ".$pass." --assembly ".$assembly." --altassembly ".$assembly." --chromosomes ".$ref_chromosomes." --altchromosomes ".
    $patch_chromosomes." --logpath ".$out_dir." --logfile fix_overlap.log --force-stage"; 

  print $command."\n";

  my $cmd_rtn = `$command`;
  print $cmd_rtn."\n";
}
print "Completed\n";
