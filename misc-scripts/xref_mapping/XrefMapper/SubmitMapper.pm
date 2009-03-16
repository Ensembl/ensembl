package XrefMapper::SubmitMapper;

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
# JOB 1 Do exonerate jobs and get core info
##################################################################

# Also load database with connection details of the core database for use later in JOBS 

# One get all info from core that is needed. :- stable_id's etc.

# Process the direct xrefs by putting them in the object_xref table

# dump the fasta files

# submit exonerate jobs (fill tables mapping, mapping_jobs) if not already done

# submit coordinate xrefs if not processed


sub new {
  my($class, $mapper) = @_;

  my $self ={};
  bless $self,$class;
  $self->core($mapper->core);
  $self->xref($mapper->xref);
  $self->mapper($mapper);
  return $self;
}


sub mapper{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_mapper} = $arg );
  return $self->{_mapper};
}


sub store_core_database_details{
  my ($self, $port, $user, $pass, $dbname, $dir);

    
}





##############################################################################
# dump fasta files code
##############################################################################
=head2 dump_seqs

  Arg[1]: xref object which holds info needed for the dump of xref

  Description: Dumps out the files for the mapping. Xref object should hold
              the value of the databases and source to be used.
  Returntype : none
  Exceptions : will die if species not known or an error occurs while
             : trying to write to files.
  Caller     : general

=cut



sub dump_seqs{

  my ($self, $location) = @_;

  $self->core->dbc->disconnect_if_idle(1);
  $self->core->dbc->disconnect_when_inactive(1);

  $self->dump_xref();

  $self->core->dbc->disconnect_if_idle(0);
  $self->core->dbc->disconnect_when_inactive(0);




  $self->xref->dbc->disconnect_when_inactive(1);
  $self->xref->dbc->disconnect_if_idle(1);

  $self->dump_ensembl($location);

  $self->xref->dbc->disconnect_if_idle(0);
  $self->xref->dbc->disconnect_when_inactive(0);

}

sub no_dump_xref {
  my ($self) = @_;

  my @method=();
  my @lists =@{$self->get_set_lists()};

  my $i=0;
  my $k = 0;
  foreach my $list (@lists){
    $method[$k++] = shift @$list;
  }
  $self->method(\@method);

  $self->core->dna_file($self->core->dir."/".$self->core->species."_dna.fasta");
  $self->core->protein_file($self->core->dir."/".$self->core->species."_protein.fasta");
}

=head2 dump_xref

  Arg[1]: xref object which holds info on method and files.

  Description: Dumps the Xref data as fasta file(s)
  Returntype : none
  Exceptions : none
  Caller     : dump_seqs

=cut

sub dump_xref{
  my ($self) = @_;

  my $xref =$self->xref();
  if(!defined($xref->dir())){
    if(defined($self->dir)){
      $xref->species($self->dir);
      $self->species_id($self->get_id_from_species_name($self->species));
    }
    else{
      $xref->dir(".");
    }
  }

  my @method=();
  my @lists =@{$self->get_set_lists()};

  my $k = 0;
  foreach my $list (@lists){
    $method[$k++] = shift @$list;
  }
  $self->method(\@method);
  
  my $i=0;
  if(defined($self->dumpcheck())){
    my $skip = 1;
    foreach my $list (@lists){
      if(!-e $xref->dir()."/xref_".$i."_dna.fasta"){
        $skip = 0;
      }
      if(!-e $xref->dir()."/xref_".$i."_peptide.fasta"){
        $skip = 0;
      }
      $i++;
    }
    if($skip){
      print "Xref fasta files found and will be used (No new dumping)\n" if($self->verbose);
      return;
    }
  }

  print "Dumping Xref fasta files\n" if($self->verbose());
  for my $sequence_type ('dna', 'peptide') {

    my $filename = $xref->dir() . "/xref_0_" . $sequence_type . ".fasta";
    open(XREF_DUMP,">$filename") || die "Could not open $filename";
    
    my $sql = "SELECT p.xref_id, p.sequence, x.species_id , x.source_id ";
    $sql   .= "  FROM primary_xref p, xref x ";
    $sql   .= "  WHERE p.xref_id = x.xref_id AND ";
    $sql   .= "        p.sequence_type ='$sequence_type' ";
    
    my $sth = $xref->dbc->prepare($sql);
    $sth->execute();
    while(my @row = $sth->fetchrow_array()){
      
      $row[1] =~ s/(.{60})/$1\n/g;
      print XREF_DUMP ">".$row[0]."\n".$row[1]."\n";
      
    }
    
    close(XREF_DUMP);
    $sth->finish();
    
  }
  my $sth = $xref->dbc->prepare("insert into process_status (status, date) values('xref_fasta_dumped',now())");
  $sth->execute();
  $sth->finish;
  
  
  return;

}


=head2 dump_ensembl

  Description: Dumps the ensembl data to a file in fasta format.
  Returntype : none
  Exceptions : none
  Caller     : dump_seqs

=cut

sub dump_ensembl{
  my ($self, $location) = @_;

  $self->fetch_and_dump_seq($location);

}

=head2 fetch_and_dump_seq

  Description: Dumps the ensembl data to a file in fasta format.
  Returntype : none
  Exceptions : wil die if the are errors in db connection or file creation.
  Caller     : dump_ensembl

=cut

sub fetch_and_dump_seq{
  my ($self) = @_;

  my $ensembl = $self->core;
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-dbconn => $ensembl->dbc);

  #
  # store ensembl dna file name and open it
  #
  if(!defined($ensembl->dir())){
    $ensembl->dir(".");
  }
  $ensembl->dna_file($ensembl->dir."/".$ensembl->species."_dna.fasta");


  #
  # store ensembl protein file name and open it
  #
  $ensembl->protein_file($ensembl->dir."/".$ensembl->species."_protein.fasta");


  if(defined($self->dumpcheck()) and -e $ensembl->protein_file() and -e $ensembl->dna_file()){
    my $sth = $self->xref->dbc->prepare("insert into process_status (status, date) values('core_fasta_dumped',now())");
    $sth->execute();    
    print "Ensembl Fasta files found (no new dumping)\n" if($self->verbose());
    return;
  }

  print "Dumping Ensembl Fasta files\n" if($self->verbose());

  open(DNA,">".$ensembl->dna_file())
    || die("Could not open dna file for writing: ".$ensembl->dna_file."\n");

  open(PEP,">".$ensembl->protein_file())
    || die("Could not open protein file for writing: ".$ensembl->protein_file."\n");

  my $gene_adaptor = $db->get_GeneAdaptor();


  # fetch by location, or everything if not defined

  my @genes;
  my $constraint;


# TEST PURPOSES ONLY#################################################
#####################################################################
  @genes = @{$gene_adaptor->fetch_all()};

#  push @genes, $gene_adaptor->fetch_by_stable_id("ENSG00000139618");
#####################################################################
  my $max = undef;

  my $i =0;
  my $rna = 0;
  foreach my $gene (@genes){
    next if $gene->biotype eq 'J_segment';
    next if $gene->biotype eq 'D_segment';

    foreach my $transcript (@{$gene->get_all_Transcripts()}) {
      $i++;
      my $seq = $transcript->spliced_seq();
      $seq =~ s/(.{60})/$1\n/g;
      print DNA ">" . $transcript->dbID() . "\n" .$seq."\n";
      my $trans = $transcript->translation();
      my $translation = $transcript->translate();

      if(defined($translation)){
        my $pep_seq = $translation->seq();
        $pep_seq =~ s/(.{60})/$1\n/g;
        print PEP ">".$trans->dbID()."\n".$pep_seq."\n";
      }
    }

     last if(defined($max) and $i > $max);

  }
  close DNA;
  close PEP;
  my $sth = $self->xref->dbc->prepare("insert into process_status (status, date) values('core_fasta_dumped',now())");
  $sth->execute();
  $sth->finish;

}

=head2 get_set_lists

  Description: specifies the list of databases and source to be used in the
             : generation of one or more data sets.
  Returntype : list of lists
  Example    : my @lists =@{$self->get_set_lists()};
  Exceptions : none
  Caller     : dump_xref

=cut

sub get_set_lists{
  my ($self) = @_;

  return [["ExonerateGappedBest1", ["*","*"]]];

}




###################################################################################################
# exonerate subs
###################################################################################################
=head2 jobcount

  Arg [1]    : (optional)
  Example    : $mapper->jobcount(1004);
  Description: Getter / Setter for number of jobs submitted.
  Returntype : scalar
  Exceptions : none

=cut

sub jobcount {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_jobcount} = $arg );
  return $self->{_jobcount};
}


sub method{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_method} = $arg );
  return $self->{_method};
}

=head2 build_list_and_map

  Arg[1]: xref object which holds info on method and files.

  Description: runs the mapping of the list of files with species methods
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub build_list_and_map {

  my ($self) = @_;

  my @list=();

  my $i = 0;

  foreach my $method (@{$self->method()}){
    my @dna=();
    my $q_dna_file = $self->xref->dir."/xref_".$i."_dna.fasta";
    if (-e $q_dna_file and -s $q_dna_file) {
      push @dna, $method;
      push @dna, $q_dna_file;
      push @dna, $self->core->dna_file();
      push @list, \@dna;
    }

    my @pep=();
    my $q_pep_file =  $self->xref->dir."/xref_".$i."_peptide.fasta";
    if (-e $q_pep_file and -s $q_pep_file) {
      push @pep, $method;
      push @pep, $self->xref->dir."/xref_".$i."_peptide.fasta";
      push @pep, $self->core->protein_file();
      push @list, \@pep;
    }
    $i++;
  }
  $self->run_mapping(\@list);

}

=head2 run_mapping

  Arg[1]     : List of lists of (method, query, target)
  Arg[2]     :
  Example    : none
  Description: Create and submit mapping jobs to LSF, and wait for them to finish.
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub run_mapping {

  my ($self, $lists) = @_;

  # delete old output files in target directory if we're going to produce new ones

  my $dir = $self->core->dir();
  print "Deleting out, err and map files from output dir: $dir\n" if($self->verbose());
  unlink (<$dir/*.map $dir/*.out $dir/*.err>);

  $self->remove_all_old_output_files();
  #disconnect so that we can then reconnect after the long mapping bit.
  $self->core->dbc->disconnect_if_idle(1);
  $self->xref->dbc->disconnect_if_idle(1);
  $self->core->dbc->disconnect_when_inactive(1);
  $self->xref->dbc->disconnect_when_inactive(1);

  # foreach method, submit the appropriate job & keep track of the job name
  # note we check if use_existing_mappings is set here, not earlier, as we
  # still need to instantiate the method object in order to fill
  # method_query_threshold and method_target_threshold

  my @job_names;
  my @running_methods;
  foreach my $list (@$lists){

    my ($method, $queryfile ,$targetfile)  =  @$list;

    my $obj_name = "XrefMapper::Methods::$method";
    # check that the appropriate object exists
    eval "require $obj_name";
    if($@) {

      warn("Could not find object $obj_name corresponding to mapping method $method, skipping\n$@");

    } else {

      my $obj = $obj_name->new();
 
      my $job_name = $obj->run($queryfile, $targetfile, $self);
      push @job_names, $job_name;
      push @running_methods, $obj;
      sleep 1; # make sure unique names really are unique
      
      $self->jobcount($self->jobcount+$obj->jobcount);
    }
  } # foreach method

  # submit depend job to wait for all mapping jobs
  foreach my $method( @running_methods ){
    # Submit all method-specific depend jobs
    if( $method->can('submit_depend_job') ){
      $method->submit_depend_job;
    }
  }
  # Submit generic depend job. Defaults to LSF
  $self->submit_depend_job($self->core->dir, @job_names);
  $self->core->dbc->disconnect_if_idle(0);
  $self->xref->dbc->disconnect_if_idle(0);
  $self->core->dbc->disconnect_when_inactive(0);
  $self->xref->dbc->disconnect_when_inactive(0);

  $self->check_err($self->core->dir);

} # run_mapping

sub nofarm{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_nofarm} = $arg );
  return $self->{_nofarm};
}

sub check_err {

  my ($self, $dir) = @_;

  foreach my $err (glob("$dir/*.err")) {

    print STDERR "\n\n*** Warning: $err has non-zero size; may indicate".
      " problems with exonerate run\n\n\n" if (-s $err);

  }
}


=head2 submit_depend_job

  Arg[1]     : List of job names.
  Arg[2]     :
  Example    : none
  Description: Submit an LSF job that waits for other jobs to finish.
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub submit_depend_job {

  my ($self, $root_dir, @job_names) = @_;


  if(defined($self->nofarm)){
    return;
  }

  # Submit a job that does nothing but wait on the main jobs to
  # finish. This job is submitted interactively so the exec does not
  # return until everything is finished.

  # build up the bsub command; first part
  my @depend_bsub = ('bsub', '-K');

  # build -w 'ended(job1) && ended(job2)' clause
  my $ended_str = "-w ";
  my $i = 0;
  foreach my $job (@job_names) {
    $ended_str .= "ended($job)";
    $ended_str .= " && " if ($i < $#job_names);
    $i++;
  }

  push @depend_bsub, $ended_str;

  # rest of command
  push @depend_bsub, ('-q', 'small', '-o', "$root_dir/depend.out", '-e', "$root_dir/depend.err");


  my $jobid = 0;

  eval {
    my $pid;
    my $reader;

    local *BSUB;
    local *BSUB_READER;

    if ( ( $reader = open( BSUB_READER, '-|' ) ) ) {
      while (<BSUB_READER>) {
        if (/^Job <(\d+)> is submitted/) {
          $jobid = $1;
          print "LSF job ID for depend job: $jobid\n" if($self->verbose);
        }
      }
      close(BSUB_READER);
    } else {
      die("Could not fork : $!\n") unless ( defined($reader) );
      open( STDERR, ">&STDOUT" );
      if ( ( $pid = open( BSUB, '|-' ) ) ) {
        print BSUB "/bin/true\n";
        close BSUB;
        if ( $? != 0 ) {
          die(   "bsub exited with non-zero status ($?) "
               . "- job not submitted\n" );
        }
      } else {
        if ( defined($pid) ) {
          exec(@depend_bsub);
          die("Could not exec bsub : $!\n");
        } else {
          die("Could not fork : $!\n");
        }
      }
      exit(0);
    }
  };

  if ($@) {
    # Something went wrong
    warn("Job submission failed:\n$@\n");
  }
  else{
    my $sth = $self->xref->dbc->prepare("insert into process_status (status, date) values('mapping_finished',now())");
    $sth->execute();
    $sth->finish;
  }	
}

sub remove_all_old_output_files{
  my ($self) =@_;

  my $dir = $self->core->dir();

  print "Deleting txt and sql files from output dir: $dir\n" if($self->verbose);
  unlink(<$dir/*.txt $dir/*.sql>);
#  $self->cleanup_projections_file();  # now to be done when we load core.
}



1;
