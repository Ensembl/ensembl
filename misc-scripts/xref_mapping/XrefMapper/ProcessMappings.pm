=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package XrefMapper::ProcessMappings;
use strict;

use vars '@ISA';
@ISA = qw{ XrefMapper::BasicMapper };

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
  $self->verbose($mapper->verbose);
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
  my $dbi = $self->xref->dbc();
  my $sth = $dbi->prepare("select job_id, percent_query_cutoff, percent_target_cutoff from mapping");
  $sth->execute();
  $sth->bind_columns(\$job_id, \$percent_query_cutoff, \$percent_target_cutoff);
  my $object_xref_id;

  while($sth->fetch){
    $query_cutoff{$job_id} = $percent_query_cutoff;    
    $target_cutoff{$job_id} = $percent_target_cutoff;
  }
  $sth->finish;

  my ($root_dir, $map, $status, $out, $err, $array_number); 
  my ($map_file, $out_file, $err_file);
  my $map_sth = $dbi->prepare("select root_dir, map_file, status, out_file, err_file, array_number, job_id from mapping_jobs");
  $map_sth->execute();
  $map_sth->bind_columns(\$root_dir, \$map, \$status, \$out, \$err, \$array_number, \$job_id);
  my $already_processed_count = 0;
  my $processed_count = 0;
  my $error_count = 0;
  my $empty_count = 0;

  my $stat_sth = $dbi->prepare("update mapping_jobs set status = ? where job_id = ? and array_number = ?");

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
	print STDERR "Problem $err_file is non zero\n";
	if(open my $eh ,"<", $err_file){
	  while(<$eh>){
	    print STDERR "#".$_;
	  }
	  close $eh;
	}
	else{
	  print STDERR "No file exists $err_file???\n Resubmit this job\n";
	}
	if($status eq "SUBMITTED"){
	  $stat_sth->execute('FAILED',$job_id, $array_number);
	}
      }
      else{ #err file checks out so process the mapping file.
	if(-e $map_file){
	  my $count = $self->process_map_file($map_file, $query_cutoff{$job_id}, $target_cutoff{$job_id}, $job_id, $array_number, $dbi);
	  if( $count > 0){
	    $processed_count++;
	    $stat_sth->execute('SUCCESS',$job_id, $array_number);
	  }
	  elsif($count ==0){
	    print STDERR "WARNING $map_file was empty could be okay but if there are alot of these it is probably a problem\n";
	    $processed_count++;
	    $empty_count++;
 	    $stat_sth->execute('SUCCESS',$job_id, $array_number);
	  }	
	  else{
	    $error_count++;
	    $stat_sth->execute('FAILED',$job_id, $array_number);
	  }	
	}
	else{
	  $error_count++;
	  print STDERR "Could not open file $map_file???\n Resubmit this job\n";
	  $stat_sth->execute('FAILED',$job_id, $array_number);
	}
      }      
    }
  }
  $map_sth->finish;
  $stat_sth->finish;

  print "already processed = $already_processed_count, processed = $processed_count, errors = $error_count, empty = $empty_count\n" if($self->verbose); 

  if(!$error_count){
    my $sth = $dbi->prepare("insert into process_status (status, date) values('mapping_processed',now())");
    $sth->execute();
    $sth->finish;
  }

}



#return number of lines parsed if succesfull. -1 for fail
sub process_map_file{
  my ($self, $map_file, $query_cutoff, $target_cutoff, $job_id, $array_number, $dbi) = @_;
  my $ret = 1;


 
  my $ensembl_type = "Translation";
  if($map_file =~ /dna_/){
    $ensembl_type = "Transcript";
  }

  my $mh;
  if(!(open $mh ,"<",$map_file) ){
    print STDERR "Could not open file $map_file\n Resubmit this job??\n";
    return -1;
  }

  my $total_lines = 0;
  my $root_dir = $self->core->dir;

  my $ins_go_sth = $dbi->prepare("insert ignore into go_xref (object_xref_id, linkage_type, source_xref_id) values(?,?,?)");
  my $dep_sth    = $dbi->prepare("select dependent_xref_id, linkage_annotation from dependent_xref where master_xref_id = ?");
  my $start_sth  = $dbi->prepare("update mapping_jobs set object_xref_start = ? where job_id = ? and array_number = ?");
  my $end_sth    = $dbi->prepare("update mapping_jobs set object_xref_end = ? where job_id = ? and array_number = ?");
#  my $update_dependent_xref_sth = $dbi->prepare("update dependent_xref set object_xref_id = ? where master_xref_id = ? and dependent_xref_id =?");

  my $object_xref_id;
  my $sth = $dbi->prepare("select max(object_xref_id) from object_xref");
  $sth->execute();
  $sth->bind_columns(\$object_xref_id);
  $sth->fetch();
  $sth->finish;
  if(!defined($object_xref_id)){
    $object_xref_id = 0;
  }
  

  my $object_xref_sth = $dbi->prepare("insert into object_xref (ensembl_id,ensembl_object_type, xref_id, linkage_type, ox_status ) values (?, ?, ?, ?, ?)");
  my $object_xref_sth2 = $dbi->prepare("insert into object_xref (ensembl_id,ensembl_object_type, xref_id, linkage_type, ox_status, master_xref_id ) values (?, ?, ?, ?, ?, ?)");
  my $get_object_xref_id_sth = $dbi->prepare("select object_xref_id from object_xref where ensembl_id = ? and ensembl_object_type = ? and xref_id = ? and linkage_type = ? and ox_status = ?");
  my $get_object_xref_id_master_sth = $dbi->prepare("select object_xref_id from object_xref where ensembl_id = ? and ensembl_object_type = ? and xref_id = ? and linkage_type = ? and ox_status = ? and master_xref_id = ?");
  local $object_xref_sth->{RaiseError}; #catch duplicates
  local $object_xref_sth->{PrintError}; # cut down on error messages
  local $object_xref_sth2->{RaiseError}; #catch duplicates
  local $object_xref_sth2->{PrintError}; # cut down on error messages

  my $identity_xref_sth = $dbi->prepare("insert ignore into identity_xref (object_xref_id, query_identity, target_identity, hit_start, hit_end, translation_start, translation_end, cigar_line, score ) values (?, ?, ?, ?, ?, ?, ?, ?, ?)");

  my $ins_dep_ix_sth = $dbi->prepare("insert ignore into identity_xref (object_xref_id, query_identity, target_identity) values(?, ?, ?)");

  my $source_name_sth = $dbi->prepare("select s.name from xref x join source s using(source_id) where x.xref_id = ?");
  my $biotype_sth = $dbi->prepare("select biotype from transcript_stable_id where internal_id = ?");

  my $last_query_id = 0;
  my $best_match_found = 0;
  my $best_identity = 0;
  my $best_score = 0;
  
  my %mRNA_biotypes = ( 'protein_coding' => 1,
			'TR_C_gene' => 1,
			'IG_V_gene' => 1,
			'nonsense_mediated_decay' => 1,
			'polymorphic_pseudogene' => 1);

  my $first = 1;
  while(<$mh>){
    my $load_object_xref = 0;
    $total_lines++;
    chomp();
    my ($label, $query_id, $target_id, $identity, $query_length, $target_length, $query_start, $query_end, $target_start, $target_end, $cigar_line, $score) = split(/:/, $_);

    if ($last_query_id != $query_id) {
	$best_match_found = 0;
	$best_score = 0;
	$best_identity = 0;
    } else {

	#ignore mappings with worse identity or score if we already found a good mapping
	if ( ($identity < $best_identity || $score < $best_score) && $best_match_found) {
	    next;
	}

    }

    if ($ensembl_type eq "Translation"){
	$load_object_xref = 1;
	
    } else {

	    #check if source name is RefSeq_ncRNA or RefSeq_mRNA
	    #if yes check biotype, if ok store object xref
	    $source_name_sth->execute($query_id);
	    my ($source_name)  = $source_name_sth->fetchrow_array;

	    if ($source_name && ($source_name =~ /^RefSeq_(m|nc)RNA/ || $source_name =~ /^miRBase/ || $source_name =~ /^RFAM/)) { 

		#make sure mRNA xrefs are matched to protein_coding biotype only
		$biotype_sth->execute($target_id);
		my ($biotype) = $biotype_sth->fetchrow_array; 
	   
		if ($source_name =~ /^RefSeq_mRNA/ && exists($mRNA_biotypes{$biotype})) {
		    $load_object_xref = 1;
		}
		if ($source_name =~ /^RefSeq_ncRNA/ && !exists($mRNA_biotypes{$biotype}) ) {
		    $load_object_xref = 1;	
		}	
                if (($source_name =~ /miRBase/ || $source_name =~ /^RFAM/) && $biotype =~ /RNA/) {
                    $load_object_xref = 1;
                }
	    } else {
		$load_object_xref = 1;
	    }
	
    }

    $last_query_id = $query_id;

    if ($score > $best_score || $identity > $best_identity) {
	$best_score = $score;
	$best_identity = $identity;
    }

    if (!$load_object_xref) {
	next;
    } else {
	$best_match_found = 1;
    }

    if(!defined($score)){
      $end_sth->execute(($object_xref_id),$job_id, $array_number);
      die "No score on line. Possible file corruption\n$_\n";      
    }

    # calculate percentage identities
    my $query_identity = int (100 * $identity / $query_length);
    my $target_identity = int (100 * $identity / $target_length);

    my $status = "DUMP_OUT";
# Only keep alignments where both sequences match cutoff
    if($query_identity < $query_cutoff or $target_identity < $target_cutoff){
      $status = "FAILED_CUTOFF";
    }		

    $object_xref_sth->execute($target_id, $ensembl_type, $query_id, 'SEQUENCE_MATCH', $status) ;
    $get_object_xref_id_sth->execute($target_id, $ensembl_type, $query_id, 'SEQUENCE_MATCH', $status);
    $object_xref_id = ($get_object_xref_id_sth->fetchrow_array())[0];
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
    if($first){
      $start_sth->execute($object_xref_id,$job_id, $array_number);
      $first = 0;
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
         $object_xref_sth2->execute($target_id, $ensembl_type, $dep_xref_id, 'DEPENDENT', $status, $master_xref_id);
         $get_object_xref_id_master_sth->execute($target_id, $ensembl_type, $dep_xref_id, 'DEPENDENT', $status, $master_xref_id);
         $object_xref_id = ($get_object_xref_id_master_sth->fetchrow_array())[0];
	 if($object_xref_sth2->err){
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
	 if($object_xref_sth2->err){
	   print STDERR "WARNING: Should not reach here??? object_xref_id = $object_xref_id\n";
	 }

	 $ins_dep_ix_sth->execute($object_xref_id, $query_identity, $target_identity);
# now store in object_xref	 $update_dependent_xref_sth->execute($object_xref_id, $master_xref_id, $dep_xref_id);

	 push @master_xref_ids, $dep_xref_id; # get the dependent, dependents just in case
	 if(defined($link) and $link ne ""){ # we have a go term linkage type
           $ins_go_sth->execute($object_xref_id, $link, $master_xref_id);
         }
       }
     }

  }	
  close $mh;
  $end_sth->execute($object_xref_id,$job_id, $array_number);
  $start_sth->finish;
  $end_sth->finish;
  $dep_sth->finish;
  $ins_go_sth->finish;
  $object_xref_sth->finish;
  $identity_xref_sth->finish;
  $ins_dep_ix_sth->finish;
  $source_name_sth->finish;
  $biotype_sth->finish;

  return $total_lines;
}




1;
