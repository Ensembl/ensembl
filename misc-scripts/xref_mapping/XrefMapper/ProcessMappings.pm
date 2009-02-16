package XrefMapper::ProcessMappings;

use vars '@ISA';
@ISA = qw{ XrefMapper::BasicMapper };

use strict;
use warnings;
use XrefMapper::BasicMapper;

use Cwd;
use DBI;
use File::Basename;
use IPC::Open3;

##################################################################
# JOB 2
##################################################################

# process exonerate mapping file (include checks)

# Save all data to object_xref. Plus remember to add dependents

# process priority xrefs to leave those excluded flagged as such

# Do tests on xref database wrt to core before saving data.

sub new {
  my($class, $mapper) = @_;

  my $self ={};
  bless $self,$class;
  $self->core($mapper->core);
  $self->xref($mapper->xref);
  return $self;
}

sub process_mappings {
  my ($self) = @_;

  # get the jobs from the mapping table
  # for i =1 i < jobnum{
  #   if( Not parsed already see mapping_jobs){
  #     check the .err, .out and .map files in that order.
  #     put data into object_xref, identity_xref, go_xref etc... 
  #     add data to mapping_job table
  #   }
  # }



  my %query_cutoff;
  my %target_cutoff;
  my ($job_id, $percent_query_cutoff, $percent_target_cutoff);
  my $sth = $self->xref->dbc->prepare("select job_id, percent_query_cutoff, percent_target_cutoff from mapping");
  $sth->execute();
  $sth->bind_columns(\$job_id, \$percent_query_cutoff, \$percent_target_cutoff);

  while($sth->fetch){
    $query_cutoff{$job_id} = $percent_query_cutoff;    
    $target_cutoff{$job_id} = $percent_target_cutoff;
  }
  $sth->finish;

  my ($root_dir, $map, $status, $out, $err, $array_number); 
  my ($map_file, $out_file, $err_file);
  my $map_sth = $self->xref->dbc->prepare("select root_dir, map_file, status, out_file, err_file, array_number, job_id from mapping_jobs");
  $map_sth->execute();
  $map_sth->bind_columns(\$root_dir, \$map, \$status, \$out, \$err, \$array_number, \$job_id);
  my $already_processed_count = 0;
  my $processed_count = 0;
  my $error_count = 0;

  my $stat_sth = $self->xref->dbc->prepare("update mapping_jobs set status = ? where job_id = ? and array_number = ?");

  while($map_sth->fetch()){
    my $err_file = $root_dir."/".$err;
    my $out_file = $root_dir."/".$out;
    my $map_file = $root_dir."/".$map;
    if($status eq "SUCCESS"){
      $already_processed_count++;
    }
    else{
      if(-s $err_file){
	$error_count++;
	print "Problem $err_file is non zero\n";
	if(open(ERR,"<$err_file")){
	  while(<ERR>){
	    print "#".$_;
	  }
	}
	else{
	  print "No file exists $err_file???\n Resubmit this job\n";
	}
	if($status eq "SUBMITTED"){
	  $stat_sth->execute('FAILED',$job_id, $array_number);
	}
      }
      else{ #err file checks out so process the mapping file.
	if(-e $map_file){
	  if($self->process_map_file($map_file, $query_cutoff{$job_id}, $target_cutoff{$job_id}, $job_id, $array_number ) >= 0){
	    $processed_count++;
	    $stat_sth->execute('SUCCESS',$job_id, $array_number);
	  }
	  else{
	    $error_count++;
	    $stat_sth->execute('FAILED',$job_id, $array_number);	    
	  }	
	}
	else{
	  $error_count++;
	  print "Could not open file $map_file???\n Resubmit this job\n";
	  $stat_sth->execute('FAILED',$job_id, $array_number);
	}
      }      
    }
  }
  $map_sth->finish;
  $stat_sth->finish;

  print "already processed = $already_processed_count, processed = $processed_count, errors = $error_count\n"; 

  if(!$error_count){
    my $sth = $self->xref->dbc->prepare("insert into process_status (status, date) values('mapping_processed',now())");
    $sth->execute();
    $sth->finish;
  }

}



#return number of lines parsed if succesfull. -1 for fail
sub process_map_file{
  my ($self, $map_file, $query_cutoff, $target_cutoff, $job_id, $array_number) = @_;
  my $ret = 1;


  my $ensembl_type = "Translation";
  if($map_file =~ /_dna_/){
    $ensembl_type = "Transcript";
  }

  if(!open(MAP ,"<$map_file")){
    print "Could not open file $map_file\n Resubmit this job??\n";
    return -1;
  }

  my $total_lines = 0;
  my $root_dir = $self->core->dir;

  my $ins_go_sth = $self->xref->dbc->prepare("insert ignore into go_xref (object_xref_id, linkage_type, source_xref_id) values(?,?,?)");
  my $dep_sth    = $self->xref->dbc->prepare("select dependent_xref_id, linkage_annotation from dependent_xref where master_xref_id = ?");
  my $start_sth  = $self->xref->dbc->prepare("update mapping_jobs set object_xref_start = ? where job_id = ? and array_number = ?");
  my $end_sth    = $self->xref->dbc->prepare("update mapping_jobs set object_xref_end = ? where job_id = ? and array_number = ?");
  my $update_dependent_xref_sth = $self->xref->dbc->prepare("update dependent_xref set object_xref_id = ? where master_xref_id = ? and dependent_xref_id =?");

  my $object_xref_id;
  my $sth = $self->xref->dbc->prepare("select max(object_xref_id) from object_xref");
  $sth->execute();
  $sth->bind_columns(\$object_xref_id);
  $sth->fetch();
  $sth->finish;

  my $object_xref_sth = $self->xref->dbc->prepare("insert into object_xref (object_xref_id, ensembl_id,ensembl_object_type, xref_id, linkage_type, ox_status ) values (?, ?, ?, ?, ?, ?)");
  local $object_xref_sth->{RaiseError}; #catch duplicates
  local $object_xref_sth->{PrintError}; # cut down on error messages

  my $identity_xref_sth = $self->xref->dbc->prepare("insert into identity_xref (object_xref_id, query_identity, target_identity, hit_start, hit_end, translation_start, translation_end, cigar_line, score ) values (?, ?, ?, ?, ?, ?, ?, ?, ?)");
 
  $start_sth->execute(($object_xref_id+1),$job_id, $array_number);
  while(<MAP>){
    $total_lines++;
    chomp();
    my ($label, $query_id, $target_id, $identity, $query_length, $target_length, $query_start, $query_end, $target_start, $target_end, $cigar_line, $score) = split(/:/, $_);


    if(!defined($score)){
      $end_sth->execute(($object_xref_id),$job_id, $array_number);
      die "No score on line. Possible file corruption\n$_\n";      
    }

    # calculate percentage identities
    my $query_identity = int (100 * $identity / $query_length);
    my $target_identity = int (100 * $identity / $target_length);

    my $status = "DUMP_OUT";
    if($query_identity < $query_cutoff and $target_identity < $target_cutoff){
      $status = "FAILED_CUTOFF";
    }		


    $object_xref_id++;
    $object_xref_sth->execute($object_xref_id, $target_id, $ensembl_type, $query_id, 'SEQUENCE_MATCH', $status) ;
    if($object_xref_sth->err){
      my $err = $object_xref_sth->errstr;
      if($err =~ /Duplicate/){
	# can get match from same xref and ensembl entity e.g.
	# ensembl/ExonerateGappedBest1_dna_569.map:xref:934818:155760:54:1617:9648:73:12:3456:3517: M 61:242
	# ensembl/ExonerateGappedBest1_dna_569.map:xref:934818:151735:58:1617:10624:73:6:5329:5397: M 48 D 1 M 19:242
	next;
      }
      else{
	$end_sth->execute(($object_xref_id),$job_id, $array_number);
	die "Problem loading error is $err\n";
      } 
    }  



    $cigar_line =~ s/ //g;
    $cigar_line =~ s/([MDI])(\d+)/$2$1/ig;


    if(!$identity_xref_sth->execute($object_xref_id, $query_identity, $target_identity, $query_start+1, $query_end, $target_start+1, $target_end, $cigar_line, $score)){
      $end_sth->execute(($object_xref_id),$job_id, $array_number);
      die "Problem loading identity_xref";
    }

     my @master_xref_ids;
     push @master_xref_ids, $query_id;
     while(my $master_xref_id = pop(@master_xref_ids)){
       my ($dep_xref_id, $link);
       $dep_sth->execute($master_xref_id);
       $dep_sth->bind_columns(\$dep_xref_id, \$link);
       while($dep_sth->fetch){
         $object_xref_id++;
         $object_xref_sth->execute($object_xref_id, $target_id, $ensembl_type, $dep_xref_id, 'DEPENDENT', $status);
	 if($object_xref_sth->err){
	   my $err = $object_xref_sth->errstr;
	   if($err =~ /Duplicate/){
#	     $duplicate_dependent_count++;
	     next;
	   }
	   else{
	     $end_sth->execute($object_xref_id,$job_id, $array_number);
	     die "Problem loading error is $err\n";
	   } 
	 }
	 $update_dependent_xref_sth->execute($object_xref_id, $master_xref_id, $dep_xref_id);
	 if($object_xref_sth->err){
	   print "WARNING: Should not reach here??? object_xref_id = $object_xref_id\n";
	 }
	 push @master_xref_ids, $dep_xref_id; # get the dependent, dependents just in case
	 if(defined($link) and $link ne ""){ # we have a go term linkage type
           $ins_go_sth->execute($object_xref_id, $link, $master_xref_id);
         }
       }
     }

  }	
  close MAP;
  $end_sth->execute($object_xref_id,$job_id, $array_number);
  $start_sth->finish;
  $end_sth->finish;
  $dep_sth->finish;
  $ins_go_sth->finish;
  $object_xref_sth->finish;
  $identity_xref_sth->finish;

  return $total_lines;
}




1;
