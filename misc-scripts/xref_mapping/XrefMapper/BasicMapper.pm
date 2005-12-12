package XrefMapper::BasicMapper;

use strict;
use Cwd;
use DBI;
use File::Basename;
use IPC::Open3;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Translation;
use XrefMapper::db;
use vars '@ISA';
use strict;

=head1 NAME

XrefMapper::BasicMapper

=head1 DESCIPTION

This is the basic mapper routine. It will create the necessary fasta files for
both the xref and ensembl sequences. These will then be matched using exonerate
and the results written to another file. By creating a <species>.pm file and 
inheriting from this base class different matching routines, parameters, data 
sets etc can be set.

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk


=cut

# Hashes to hold method-specific thresholds
my %method_query_threshold;
my %method_target_threshold;

# Various useful variables.
my %translation_to_transcript;
my %transcript_to_translation;
my %genes_to_transcripts;
my %xref_to_source;
my %object_xref_mappings;
my %object_xref_identities;
my %xref_descriptions;
my %xref_accessions;
my %xref_labels;
my %source_to_external_db;
my %xrefs_written;
my %object_xrefs_written;
my %failed_xref_mappings;
my %updated_source;

=head2 new

  Description: Constructor for BasicMapper.
  Returntype : BasicMapper
  Exceptions : none
  Caller     : general

=cut

sub new{
  my($class, @args) = @_;

  my $self ={};
  bless $self,$class;
  $self->jobcount(0);

  return $self;
}


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

  $self->dump_xref();
  $self->dump_ensembl($location);

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
    push @dna, $method;
    push @dna, $self->xref->dir."/xref_".$i."_dna.fasta";
    push @dna, $self->core->dna_file();
    push @list, \@dna;
    my @pep=();
    push @pep, $method;
    push @pep, $self->xref->dir."/xref_".$i."_peptide.fasta";
    push @pep, $self->core->protein_file();
    push @list, \@pep;
    $i++;
  }

  $self->run_mapping(\@list);

}


=head2 get_species_id_from_species_name

  Arg[1]: species name

  Description: get the species_id from the database for the named database.
  Example    : my $id = get_species_id_from_species_name('homo_sapiens');
  Returntype : int (species_id)
  Exceptions : will die if species does not exist in given xref database.
  Caller     : general

=cut

sub get_species_id_from_species_name{
  my ($self,$species) = @_;


  my $sql = "select species_id from species where name = '".$species."'";
  my $sth = $self->dbc->prepare($sql);
  $sth->execute();
  my @row = $sth->fetchrow_array();
  my $species_id;
  if (@row) {
    $species_id = $row[0];
  } else {
    print STDERR "Couldn't get ID for species ".$species."\n";
    print STDERR "It must be one of :-\n";
    $sql = "select name from species";
    $sth = $self->dbc->prepare($sql);
    $sth->execute();
    while(my @row = $sth->fetchrow_array()){
      print STDERR $row[0]."\n";
    }
    die("Please try again :-)\n");
  }
  $sth->finish();

  return $species_id;
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

=head2 get_source_id_from_source_name

  Arg[1]: source name

  Description: get the source_id from the database for the named source.
  Example    : my $id = get_source_id_from_source_name('RefSeq');
  Returntype : int (source_id)
  Exceptions : will die if source does not exist in given xref database.
  Caller     : general

=cut

sub get_source_id_from_source_name{
  my ($self, $source) = @_;
  my $source_id;
  
  my $sql = "select source_id from source where name = '".$source."'";
  my $sth = $self->dbc->prepare($sql);
  $sth->execute();
  my @row = $sth->fetchrow_array();
  if (defined $row[0] and $row[0] ne '') {
    $source_id = $row[0];
  } else {
    print STDERR "Couldn't get ID for source ".$source."\n";
    print STDERR "It must be one of :-\n";
    $sql = "select name from source";
    $sth = $self->dbc->prepare($sql);
    $sth->execute();
    while(my @row = $sth->fetchrow_array()){
      print STDERR $row[0]."\n";
    }
    die("Please try again :-)\n");
  }
  $sth->finish();

  return $source_id;
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
    }
    else{
      $xref->dir(".");
    }
  }
  
  my @method=();
  
  my @lists =@{$self->get_set_lists()};
  
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
      my $k = 0;
      foreach my $list (@lists){
	$method[$k++] = shift @$list;
      }
      $self->method(\@method);
      return;
    }
  }

  $i=0;
  foreach my $list (@lists){
    $method[$i] = shift @$list;
    my $j = 0;
    my @source_id=();
    my @species_id=();
    foreach my $element (@$list){
      while(my $species = shift(@$element)){
	if($species ne "*"){
	  $species_id[$j] = get_species_id_from_species_name($xref,$species);
	}
	else{
	  $species_id[$j] = -1;
	}
	my $source = shift(@$element);
	if($source ne "*"){
	  $source_id[$j] = get_source_id_from_source_name($xref,$source);
	}
	else{
	  $source_id[$j] = -1;
	}
	$j++;
      }
    }
    #method data fully defined now
    $self->dump_subset($xref,\@species_id,\@source_id,$i);    
    $i++;
  }
  
  $self->method(\@method);

  return;
  
}

=head2 dump_subset

  Arg[1]: xref object which holds info on files.
  Arg[2]: list of species to use.
  Arg[3]: list of sources to use.
  Arg[4]: index to be used in file creation.
  
  Description: Dumps the Xref data for one set of species/databases
  Returntype : none
  Exceptions : none
  Caller     : dump_xref

=cut


sub dump_subset{

  my ($self,$xref,$rspecies_id,$rsource_id,$index) = @_;

  # generate or condition list for species and sources
  my $final_clause;
  my $use_all = 0;
  my @or_list;
  for (my $j = 0; $j < scalar(@$rspecies_id); $j++){
    my @condition;
    if($$rspecies_id[$j] > 0){
      push @condition, "x.species_id=" . $$rspecies_id[$j];
   }
    if($$rsource_id[$j] > 0){
      push @condition, "x.source_id=" . $$rsource_id[$j];
    }

    # note if both source and species are * (-1) there's no need for a final clause

    if ( !@condition ) {
      $use_all = 1;
      last;
    }

    push @or_list, join (" AND ", @condition);

  }

  $final_clause = " AND ((" . join(") OR (", @or_list) . "))" unless ($use_all) ;


  for my $sequence_type ('dna', 'peptide') {

    my $filename = $xref->dir() . "/xref_" . $index . "_" . $sequence_type . ".fasta";
    open(XREF_DUMP,">$filename") || die "Could not open $filename";

    my $sql = "SELECT p.xref_id, p.sequence, x.species_id , x.source_id ";
    $sql   .= "  FROM primary_xref p, xref x ";
    $sql   .= "  WHERE p.xref_id = x.xref_id AND ";
    $sql   .= "        p.sequence_type ='$sequence_type' ";
    $sql   .= $final_clause;

    if(defined($self->maxdump())){
      $sql .= " LIMIT ".$self->maxdump()." ";
    }

    my $sth = $xref->dbc->prepare($sql);
    $sth->execute();
    while(my @row = $sth->fetchrow_array()){

      $row[1] =~ s/(.{60})/$1\n/g;
      print XREF_DUMP ">".$row[0]."\n".$row[1]."\n";

    }

    close(XREF_DUMP);
    $sth->finish();

  }

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
  my ($self, $location) = @_;

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
    return;
  }
  open(DNA,">".$ensembl->dna_file()) 
    || die("Could not open dna file for writing: ".$ensembl->dna_file."\n");

  open(PEP,">".$ensembl->protein_file()) 
    || die("Could not open protein file for writing: ".$ensembl->protein_file."\n");

  my $gene_adaptor = $db->get_GeneAdaptor();


  # fetch by location, or everything if not defined
  my @genes;
  if ($location) {

    my $slice_adaptor = $db->get_SliceAdaptor();
    my $slice = $slice_adaptor->fetch_by_name($location);
    @genes = @{$gene_adaptor->fetch_all_by_Slice($slice)};

  } else {

    @genes = @{$gene_adaptor->fetch_all()};

  }

  my $max = undef;
  if(defined($self->maxdump())){
    $max = $self->maxdump();
  }
  my $i =0;
  foreach my $gene (@genes){
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

}



###
# Getter/Setter methods
###





=head2 method
 
  Arg [1]    : (optional) list reference $arg
               reference to a list of method names 
  Example    : my @methods = @{$self->method()};
  Description: Getter / Setter for the methods 
  Returntype : list
  Exceptions : none

=cut


sub method{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_method} = $arg );
  return $self->{_method};
}

=head2 core
 
  Arg [1]    : (optional) 
  Example    : $mapper->core($new_core);
  Description: Getter / Setter for the core. 
               info for the ensembl core database. 
  Returntype : XrefMapper::db
  Exceptions : none

=cut

sub core{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_core} = $arg );
  return $self->{_core};
}


=head2 dumpcheck
 
  Arg [1]    : (optional) 
  Example    : $mapper->dumpcheck("yes");
  Description: Getter / Setter for dumpcheck. 
               If set the mapper will not dump fasta files 
               if they exist already. 
  Returntype : scalar
  Exceptions : none

=cut

sub dumpcheck {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_dumpcheck} = $arg );
  return $self->{_dumpcheck};
}


=head2 maxdump
 
  Arg [1]    : (optional) 
  Example    : $mapper->maxdump(10);
  Description: Getter / Setter for maxdump. 
               If set the mapper will only dump that number of 
               sequences into the fasta files. (Mainly used for testing).
  Returntype : scalar
  Exceptions : none

=cut

sub maxdump {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_maxdump} = $arg );
  return $self->{_maxdump};
}



=head2 use_existing_mappings

  Arg [1]    : (optional) 
  Example    : $mapper->use_existing_mappings("yes");
  Description: Getter / Setter for use_existing_mappings. 
               If set the mapper will not redo the mapping
               but parse the existing .map files.
  Returntype : scalar
  Exceptions : none

=cut

sub use_existing_mappings {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_use_existing_mappings} = $arg );
  return $self->{_use_existing_mappings};
}


=head2 xref
 
  Arg [1]    : (optional) 
  Example    : $mapper->core($new_core);
  Description: Getter / Setter for the core. 
               info for the xref database. 
  Returntype : XrefMapper::db
  Exceptions : none

=cut

sub xref{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_xref} = $arg );
  return $self->{_xref};
}

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
  if (!defined($self->use_existing_mappings())) {
    print "Deleting out err and map files from output dir\n";
    my $dir = $self->core->dir();
    unlink (<$dir/*.map $dir/*.out $dir/*.err>);
  }
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

  foreach my $list (@$lists){

    my ($method, $queryfile ,$targetfile)  =  @$list;

    my $obj_name = "XrefMapper::Methods::$method";
    # check that the appropriate object exists
    eval "require $obj_name";
    if($@) {

      warn("Could not find object $obj_name corresponding to mapping method $method, skipping\n$@");

    } else {

      my $obj = $obj_name->new();
      $method_query_threshold{$method} = $obj->query_identity_threshold();
      $method_target_threshold{$method} = $obj->target_identity_threshold();


      if (!defined($self->use_existing_mappings)) {
	my $job_name = $obj->run($queryfile, $targetfile, $self->core->dir());
	push @job_names, $job_name;
	sleep 1; # make sure unique names really are unique
      }
      $self->jobcount($self->jobcount+$obj->jobcount);
    }

  } # foreach method


  if (!defined($self->use_existing_mappings)) {
    # submit depend job to wait for all mapping jobs

    submit_depend_job($self->core->dir, @job_names);
  }
  $self->check_err($self->core->dir); 
} # run_mapping


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

  my ($root_dir, @job_names) = @_;

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

  #print "##depend bsub:\n" . join (" ", @depend_bsub) . "\n";

  my $jobid = 0;

  eval {
    my $pid;
    my $reader;

    local *BSUB;
    local *BSUB_READER;

    if (($reader = open(BSUB_READER, '-|'))) {
      while (<BSUB_READER>) {
	if (/^Job <(\d+)> is submitted/) {
	  $jobid = $1;
	  print "LSF job ID for depend job: $jobid\n"
	}
      }
      close(BSUB_READER);
    } else {
      die("Could not fork : $!\n") unless (defined($reader));
      open(STDERR, ">&STDOUT");
      if (($pid = open(BSUB, '|-'))) {
	
	print BSUB "/bin/true\n";
	close BSUB;
	if ($? != 0) {
	  die("bsub exited with non-zero status ($?) - job not submitted\n");
	}
      } else {
	if (defined($pid)) {
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

}

=head2 parse_mappings

  Example    : none
  Description: Parse exonerate output files and build files for loading into target db tables.
  Returntype : List of strings
  Exceptions : none
  Caller     : general

=cut

sub parse_mappings {

  my ($self, $notriage) = @_;

  my $ensembl = $self->core;
  my $xref = $self->xref;
  my $dir = $ensembl->dir();

  # incase timed out, force reconnection
  $ensembl->dbc->disconnect_if_idle(0);
  $ensembl->dbc->disconnect_when_inactive(0);
  $ensembl->dbc->connect();

  $xref->dbc->disconnect_if_idle(0);
  $xref->dbc->disconnect_when_inactive(0);
  $xref->dbc->connect();

  # cache xref id->label info; useful for debugging
#  $self->cache_xref_labels($xref->dbc);

  # get current max object_xref_id
  my $row = @{$ensembl->dbc->db_handle->selectall_arrayref("SELECT MAX(object_xref_id) FROM object_xref")}[0];
  my $max_object_xref_id = @{$row}[0];
  if (!defined $max_object_xref_id) {
    print "No existing object_xref_ids, will start from 1\n";
    $max_object_xref_id = 0;
  } else {
    print "Maximum existing object_xref_id = $max_object_xref_id\n";
  }
  my $object_xref_id_offset = $max_object_xref_id + 1;
  my $object_xref_id = $object_xref_id_offset;

  $row = @{$ensembl->dbc->db_handle->selectall_arrayref("SELECT MAX(xref_id) FROM xref")}[0];
  my $max_xref_id = @$row[0];
  if (!defined $max_xref_id) {
    print "No existing xref_ids, will start from 1\n";
    $max_xref_id = -1; # so that generated xref_ids will be the same as in xref db
  } else {
    print "Maximum existing xref_id = $max_xref_id\n";
  }

  my $xref_id_offset = $max_xref_id + 1;

  $self->xref_id_offset($xref_id_offset); #store 

  # files to write table data to
  open (OBJECT_XREF,   ">$dir/object_xref.txt");
  open (IDENTITY_XREF, ">$dir/identity_xref.txt");

  my $total_lines = 0;
  my $last_lines = 0;
  my $total_files = 0;

  # keep a (unique) list of xref IDs that need to be written out to file as well
  # this is a hash of hashes, keyed on xref id that relates xrefs to e! objects (may be 1-many)
  my %primary_xref_ids = ();


  # and a list of mappings of ensembl objects to xrefs
  # (primary now, dependent added in dump_core_xrefs)
  # this is required for display_xref generation later
  # format:
  #   key: ensembl object type:ensembl object id
  #   value: list of xref_id (with offset)
  # Note %object_xref_mappings is global


  my @dna_check=();
  my @pep_check=();

  foreach my $file (glob("$dir/*.map")) {

    #print "Parsing results from " . basename($file) .  "\n";
    open(FILE, $file);
    $total_files++;

    # files are named Method_(dna|peptide)_N.map
#    my $type = get_ensembl_object_type($file);
#
#    my $method = get_method($file);
    my ($method, $type, $part) = get_parts($file);
    
    if($type =~ 'Translation'){
      $pep_check[$part] = $part;
    }
    elsif($type =~ 'Transcript'){
      $dna_check[$part] =  $part;
    }
    else{
      die "unknown type $type\n";
    }
    # get or create the appropriate analysis ID
    # XXX restore when using writeable database
    my $analysis_id = $self->get_analysis_id($type);
    #    my $analysis_id = 999;

    while (<FILE>) {

      $total_lines++;
      chomp();
      my ($label, $query_id, $target_id, $identity, $query_length, $target_length, $query_start, $query_end, $target_start, $target_end, $cigar_line, $score) = split(/:/, $_);
      $cigar_line =~ s/ //g;

      # calculate percentage identities
      my $query_identity = int (100 * $identity / $query_length);
      my $target_identity = int (100 * $identity / $target_length);

      # only take mappings where there is a good match on one or both sequences
      if ($query_identity  < $method_query_threshold{$method} &&
	  $target_identity < $method_target_threshold{$method}){
	my $reason = $target_id."|".$type."|".$query_identity."|".$target_identity."|";
	   $reason .= $method_query_threshold{$method}."|". $method_target_threshold{$method};
	$failed_xref_mappings{$query_id} = $reason;
	next;
      }

      # note we add on $xref_id_offset to avoid clashes
      print OBJECT_XREF "$object_xref_id\t$target_id\t$type\t" . ($query_id+$xref_id_offset) . "\n";
      print IDENTITY_XREF join("\t", ($object_xref_id, $query_identity, $target_identity, $query_start+1, $query_end, $target_start+1, $target_end, $cigar_line, $score, "\\N", $analysis_id)) . "\n";

      # TODO - evalue?
      $object_xref_id++;

      # store mapping for later - note NON-OFFSET xref_id is used
      my $key = $type . "|" . $target_id;
      my $xref_id = $query_id;
      push @{$object_xref_mappings{$key}}, $xref_id;

      # store query & target identities
      # Note this is a hash (type|object id) of hashes (xref id) of hashes ("query_identity" or "target_identity")
      $object_xref_identities{$key}->{$xref_id}->{"query_identity"} = $query_identity;
      $object_xref_identities{$key}->{$xref_id}->{"target_identity"} = $target_identity;

      # note the NON-OFFSET xref_id is stored here as the values are used in
      # a query against the original xref database
      $primary_xref_ids{$query_id}{$target_id."|".$type} = $target_id."|".$type;

    }

    close(FILE);
    #print "After $file, lines read increased by " . ($total_lines-$last_lines) . "\n";
    $last_lines = $total_lines;
  }

  close(IDENTITY_XREF);
  close(OBJECT_XREF);

  print "Read $total_lines lines from $total_files exonerate output files\n";

  if($self->jobcount() != $total_files){
    print $dna_check[-1]." dna map files\n";
    print $pep_check[-1]." peptide map files\n";
    my $test_failed = 0;
    for(my $i=1; $i < $dna_check[-1]; $i++){
      if($dna_check[$i] != $i){
	print "DNA $i file not found\n"; 
	$test_failed = 1;
      }
    }
    for(my $i=1; $i < $pep_check[-1]; $i++){
      if($pep_check[$i] != $i){
	print "PEPTIDE $i file not found\n"; 
	$test_failed = 1;
      }
    }
    if($test_failed){
      die "Missing Files aborting run\n";
    }
    if(!defined($self->use_existing_mappings())){
      print  "There should be ".$self->jobcount()." files. Please check\n";
      print  "As this is the number of jobs submitted\n";
      die "There should be ".$self->jobcount()." files. Please check\n";
    }  
  }

  # write relevant xrefs to file
  $max_object_xref_id = $self->dump_core_xrefs(\%primary_xref_ids, $object_xref_id+1, $xref_id_offset, $object_xref_id_offset);


  # dump interpro table as well
  $self->dump_interpro();

  # dump direct xrefs
  $self->dump_direct_xrefs($xref_id_offset, $max_object_xref_id);

  # dump xrefs that don't appear in either the primary_xref or dependent_xref tables
  $self->dump_orphan_xrefs($xref_id_offset);


  # dump triage type data
  $self->dump_triage_data($xref_id_offset) if (!$notriage);


  # write comparison info. Can be removed after development
  ###writes to xref.txt.Do not want to do this if loading data afterwards
  ####  $self->dump_comparison();

}

sub get_stable_ids(){
  my ($self, $type, $string, $hashref) = @_;

  my $sql = "SELECT ".$type."_id ,stable_id ";
  $sql .=      "FROM ".$type."_stable_id ";
  $sql .=          "WHERE ".$type."_id IN (".$string.")";

  my $sth = $self->core->dbc->prepare($sql);
  $sth->execute();
  my ($trans, $stable);
  $sth->bind_columns(\$trans,\$stable);
  while($sth->fetch()){
    $hashref->{$trans} = $stable;
  }
  $sth->finish;
}

sub dump_triage_data() {
  my ($self, $xref_id_offset) = @_;

  print "Dumping triage data\n";
  my $translation="";
  my $translation_count=0;
  my $transcript="";
  my $transcript_count=0;
  my $batch_size=200;
  my %translation_2_stable=();
  my %transcript_2_stable=();
  foreach my $temp (values %failed_xref_mappings){
    my ($id, $type) = split(/\|/,$temp);
    if($type =~ /Translation/){
      $translation_count++;
      if($translation_count > $batch_size){
	my $ex = $self->get_stable_ids("translation",$translation.$id,\%translation_2_stable);
	$translation_count = 0;
	$translation = "";
      }
      else{
	$translation .= "$id,";
      }
    }
    elsif($type =~ /Transcript/){
      $transcript_count++;
      if($transcript_count > $batch_size){
	$self->get_stable_ids("transcript",$transcript.$id,\%transcript_2_stable);
	$transcript_count=0;
	$transcript="";
      }
      else{
	$transcript .= "$id,";
      }
    }
    else{
      die "Unknown type *".$type."*\n".$temp."\n";
    }
  }
  if($transcript_count){
    chop $transcript; # remove last , from list
    $self->get_stable_ids("transcript",$transcript,\%transcript_2_stable);
  }
  if($translation_count){
    chop $translation;
    $self->get_stable_ids("translation",$translation,\%translation_2_stable);
  }


  my $primary_sql= (<<PSQL);
    SELECT DISTINCT(s.source_id) 
      FROM source s, primary_xref px, xref x 
	WHERE x.xref_id = px.xref_id
	  AND s.source_id = x.source_id
PSQL

  my $psth = $self->xref->dbc->prepare($primary_sql) || die "prepare failed";
  $psth->execute() || die "execute failed";

  my @primary_sources =();

  my ($prim);
  $psth->bind_columns(\$prim);
  while($psth->fetch()){
    push @primary_sources, $prim;
  }

  open (XREF, ">>" . $self->core->dir() . "/xref.txt");
  open (XREF_MISSED, ">" . $self->core->dir() . "/triage.txt");

#  foreach my $source ("%RefSeq%","UniProt%","Anopheles_symbol"){
  foreach my $source (@primary_sources){
    my $sql = "select x.xref_id, x.accession, x.version, x.label, x.description, x.source_id, x.species_id from xref x, source s where x.source_id = s.source_id";
#    my $sql = "select x.xref_id, x.accession, x.version, x.label, x.description, x.source_id, x.species_id from xref x, source s where s.name like '".$source."' and x.source_id = s.source_id";
    my $sth = $self->xref->dbc->prepare($sql);
    $sth->execute();
    
    my ($xref_id, $accession, $version, $label, $description, $source_id, $species_id);
    $sth->bind_columns(\$xref_id, \$accession, \$version, \$label, \$description, \$source_id, \$species_id);
    while($sth->fetch()){
      if (!$xrefs_written{$xref_id}) {
#	$xrefs_written{$xref_id} = 1;
	my $external_db_id = $source_to_external_db{$source_id};
	if(!defined($updated_source{$external_db_id})){
	  $self->cleanup_sources_file($external_db_id);
	}
	print XREF ($xref_id+$xref_id_offset) . "\t" . $external_db_id . "\t" . $accession . "\t" . $label . "\t" . $version . "\t" . $description . "\n";

#dump out dependencies aswell
	
        $self->dump_all_dependencies($xref_id, $xref_id_offset);

	if(defined($failed_xref_mappings{$xref_id})){
	  my ($ensembl_id,$type,$q_perc,$t_perc,$q_cut,$t_cut) =  split(/\|/,$failed_xref_mappings{$xref_id});
	  print XREF_MISSED  ($xref_id+$xref_id_offset) . "\thighest match is $accession (".$q_perc."%) to ";
	  if($type  =~ /Translation/){
	    print XREF_MISSED $translation_2_stable{$ensembl_id};
	  }
	  elsif($type  =~ /Transcript/){
	    print XREF_MISSED $transcript_2_stable{$ensembl_id};
	  }
	  else{
	    die "type=*".$type."*\n".$failed_xref_mappings{$xref_id}."\n";
	  }
	  print XREF_MISSED " (".$t_perc."\%) which are below their respective cutoffs off $q_cut\% and $t_cut\%.\n";
	}
	else{
	  print XREF_MISSED  ($xref_id+$xref_id_offset) . "\t$accession no match to Ensembl\n";
	}
      }
    }
    $sth->finish;
  }
  close(XREF);
  close(XREF_MISSED);
}

# dump xrefs that don't appear in either the primary_xref or dependent_xref tables
# e.g. Interpro xrefs

sub dump_orphan_xrefs() {

  my ($self, $xref_id_offset) = @_;

  my $count;

  open (XREF, ">>" . $self->core->dir() . "/xref.txt");
  open (XREF_MISSED, ">>" . $self->core->dir() . "/xref_failed.txt");

  # need a double left-join
  my $sql = "SELECT x.xref_id, x.accession, x.version, x.label, x.description, x.source_id, x.species_id FROM xref x LEFT JOIN primary_xref px ON px.xref_id=x.xref_id LEFT JOIN dependent_xref dx ON dx.dependent_xref_id=x.xref_id WHERE px.xref_id IS NULL AND dx.dependent_xref_id IS NULL";

  my $sth = $self->xref->dbc->prepare($sql);
  $sth->execute();

  my ($xref_id, $accession, $version, $label, $description, $source_id, $species_id);
  $sth->bind_columns(\$xref_id, \$accession, \$version, \$label, \$description, \$source_id, \$species_id);

  while ($sth->fetch()) {

    my $external_db_id = $source_to_external_db{$source_id};
    if ($external_db_id) { # skip "unknown" sources
      if (!$xrefs_written{$xref_id}) {
	if(!defined($updated_source{$external_db_id})){
	  $self->cleanup_sources_file($external_db_id);
	}
	print XREF ($xref_id+$xref_id_offset) . "\t" . $external_db_id . "\t" . $accession . "\t" . $label . "\t" . $version . "\t" . $description . "\n";

#	$xrefs_written{$xref_id} = 1;
	$count++;
      }
    }

  }
  $sth->finish();

  close(XREF);
  close(XREF_MISSED);

  print "Wrote $count xrefs that are neither primary nor dependent\n";

}


# Dump direct xrefs. Need to do stable ID -> internal ID mapping.

sub dump_direct_xrefs {

  my ($self, $xref_id_offset, $max_object_xref_id) = @_;
  my $object_xref_id = $max_object_xref_id + 1;

  print "Writing direct xrefs\n";

  my $count = 0;

  open (XREF, ">>" . $self->core->dir() . "/xref.txt");
  open (OBJECT_XREF, ">>" . $self->core->dir() . "/object_xref.txt");
  open (GO_XREF, ">>" .$self->core->dir(). "/go_xref.txt");


  my $go_source_id;
  my $worm_source_id = undef;
  # Will need to look up translation stable ID from transcript stable ID, build hash table
  print "Building transcript stable ID -> translation stable ID lookup table\n";
  my %transcript_stable_id_to_translation_stable_id;
  my $trans_sth = $self->core->dbc->prepare("SELECT tss.stable_id as transcript, tls.stable_id AS translation FROM translation tl, translation_stable_id tls, transcript_stable_id tss WHERE tss.transcript_id=tl.transcript_id AND tl.translation_id=tls.translation_id");
  $trans_sth->execute();
  my ($transcript_stable_id, $translation_stable_id);
  $trans_sth->bind_columns(\$transcript_stable_id, \$translation_stable_id);
  while ($trans_sth->fetch()) {
    $transcript_stable_id_to_translation_stable_id{$transcript_stable_id} = $translation_stable_id;
  }
  $trans_sth->finish();

  # Will need lookup tables for gene/transcript/translation stable ID to internal ID
  my $stable_id_to_internal_id = $self->build_stable_id_to_internal_id_hash();


  # SQL / statement handle for getting all direct xrefs
  my $xref_sql = "SELECT dx.general_xref_id, dx.ensembl_stable_id, dx.type, dx.linkage_xref, x.accession, x.version, x.label, x.description, x.source_id, x.species_id FROM direct_xref dx, xref x WHERE dx.general_xref_id=x.xref_id";
  my $xref_sth = $self->xref->dbc->prepare($xref_sql);

  $xref_sth->execute();

  my ($xref_id, $ensembl_stable_id, $type, $linkage_xref, $accession, $version, $label, $description, $source_id, $species_id);
  $xref_sth->bind_columns(\$xref_id, \$ensembl_stable_id, \$type, \$linkage_xref,\ $accession, \$version, \$label, \$description, \$source_id, \$species_id);

  while ($xref_sth->fetch()) {

    my $external_db_id = $source_to_external_db{$source_id};
    if ($external_db_id) {

      # In the case of CCDS xrefs, direct_xref is to transcript but we want
      # the mapping in the core db to be to the *translation*
      if ($source_id == get_source_id_from_source_name($self->xref(), "CCDS")) {
	$type = 'translation';
	my $tmp_esid = $ensembl_stable_id;
	$ensembl_stable_id = $transcript_stable_id_to_translation_stable_id{$tmp_esid};
	warn "Can't find translation for transcript $tmp_esid" if (!$ensembl_stable_id);
	#print "CCDS: transcript $tmp_esid -> translation $ensembl_stable_id\n";
      }

      my $ensembl_internal_id = $stable_id_to_internal_id->{$type}->{$ensembl_stable_id};

      if ($ensembl_internal_id) {

	if (!$xrefs_written{$xref_id}) {
	  if(!defined($updated_source{$external_db_id})){
	    $self->cleanup_sources_file($external_db_id);
	  }
	  print XREF ($xref_id+$xref_id_offset) . "\t" . $external_db_id . "\t" . $accession . "\t" . $label . "\t" . $version . "\t" . $description . "\n";
	  $xrefs_written{$xref_id} = 1;
	}
	print OBJECT_XREF "$object_xref_id\t$ensembl_internal_id\t" . ucfirst($type) . "\t" . ($xref_id+$xref_id_offset) . "\n";
	$object_xref_id++;
	$count++;

      } else {
	if(!defined($worm_source_id)){
	  $worm_source_id = get_source_id_from_source_name($self->xref(), "wormpep_id");
	  $go_source_id = get_source_id_from_source_name($self->xref(), "GO" );
	}
	# deal with UTR transcripts in Elegans and potentially others
	# Need to link xrefs that are listed as linking to e.g. ZK829.4
	# to each of ZK829.4.1, ZK829.4.2, ZK829.4.3
	my $old_object_xref_id = $object_xref_id;
	if ($source_id == $worm_source_id || $source_id == $go_source_id) {

	  # search for matching stable IDs
	  my $pat = $ensembl_stable_id .  '\..+';
	  foreach my $stable_id (keys %{$stable_id_to_internal_id->{$type}}) {

	    if ($stable_id =~ /$pat/) {

	      if (!$xrefs_written{$xref_id}) {
		if(!defined($updated_source{$external_db_id})){
		  $self->cleanup_sources_file($external_db_id);
		}
		print XREF ($xref_id+$xref_id_offset) . "\t" . $external_db_id . "\t" . $accession . "\t" . $label . "\t" . $version . "\t" . $description . "\n";
		$xrefs_written{$xref_id} = 1;
	      }
	      $ensembl_internal_id = $stable_id_to_internal_id->{$type}->{$stable_id};
	      print OBJECT_XREF "$object_xref_id\t$ensembl_internal_id\t" . ucfirst($type) . "\t" . ($xref_id+$xref_id_offset) . "\n";
	      if( $source_id == $go_source_id){
		print GO_XREF $object_xref_id . "\t" . $linkage_xref . "\n";
	      }
	      $object_xref_id++;

	    }

	  } # foreach stable_id

	} # if source_id

	# if we haven't changed $object_xref_id, nothing was written
	print STDERR "Can't find $type corresponding to stable ID $ensembl_stable_id in ${type}_stable_id, not writing record for xref $accession\n" if ($object_xref_id == $old_object_xref_id);

      }

     

    }

  }

  close(OBJECT_XREF);
  close(XREF);
  close(GO_XREF);

  $xref_sth->finish();

  print "Wrote $count direct xrefs\n";

}


# Dump the interpro table from the xref database
sub dump_interpro {

  my $self = shift;

  open (INTERPRO, ">" .  $self->core->dir() . "/interpro.txt");

  my $sth = $self->xref->dbc->prepare("SELECT * FROM interpro");
  $sth->execute();

  my ($interpro, $pfam);
  $sth->bind_columns(\$interpro, \$pfam);
  while ($sth->fetch()) {
    print INTERPRO $interpro . "\t" . $pfam . "\n";
  }
  $sth->finish();

  close (INTERPRO);

}

sub build_stable_id_to_internal_id_hash {

  my ($self) = @_;

  my %stable_id_to_internal_id;

  foreach my $type ('gene', 'transcript', 'translation') { # Add exon here if required

    print "Caching stable ID -> internal ID links for ${type}s\n";

    my $core_sql = "SELECT ${type}_id, stable_id FROM ${type}_stable_id" ;
    my $sth = $self->core->dbc->prepare($core_sql);
    $sth->execute();
    my ($internal_id, $stable_id);
    $sth->bind_columns(\$internal_id, \$stable_id);

    while ($sth->fetch) {

      $stable_id_to_internal_id{$type}{$stable_id} = $internal_id;

    }

  }

  return \%stable_id_to_internal_id;

}

sub get_ensembl_object_type {

  my $filename = shift;
  my $type;

  $filename = basename($filename);

  if ($filename =~ /_dna_/i) {

    $type = "Transcript";

  } elsif ($filename =~ /_peptide_/i) {

    $type = "Translation";

  } else {

    print STDERR "Cannot deduce Ensembl object type from filename $filename\n";
  }

  return $type;

}

sub get_parts {

  my $filename = shift;

  $filename = basename($filename);

  my ($method,$type,$part) = $filename =~ /^(.*)_(dna|peptide)_(\d+)\.map/;
  if ($type eq "dna" ) {

    $type = "Transcript";

  } elsif ($type eq "peptide") {

    $type = "Translation";

  } else {

    print STDERR "Cannot deduce Ensembl object type from filename $filename\n";
  }

  return ($method, $type, $part);

}

sub get_method {

  my $filename = shift;

  $filename = basename($filename);

  my ($method) = $filename =~ /^(.*)_(dna|peptide)_\d+\.map/;

  return $method;

}

sub get_analysis_id {

  my ($self, $ensembl_type) = @_;

  my %typeToLogicName = ( 'transcript' => 'XrefExonerateDNA',
			  'translation' => 'XrefExonerateProtein' );

  my $logic_name = $typeToLogicName{lc($ensembl_type)};

  my $sth = $self->core->dbc->prepare("SELECT analysis_id FROM analysis WHERE logic_name='" . $logic_name ."'");
  $sth->execute();

  my $analysis_id;

  if (my @row = $sth->fetchrow_array()) {

    $analysis_id = $row[0];
#    print "Found exising analysis ID ($analysis_id) for $logic_name\n";

  } else {

    print "No analysis with logic_name $logic_name found, creating ...\n";
    $sth = $self->core->dbc->prepare("INSERT INTO analysis (logic_name, created) VALUES ('" . $logic_name. "', NOW())");
    # TODO - other fields in analysis table
    $sth->execute();
    $analysis_id = $sth->{'mysql_insertid'};
    print "Done (analysis ID=" . $analysis_id. ")\n";

  }
  $sth->finish();

  return $analysis_id;

}


sub dump_core_xrefs {

  my ($self, $xref_ids_hashref, $start_object_xref_id, $xref_id_offset, $object_xref_id_offset) = @_;

  my @xref_ids = keys %$xref_ids_hashref;
  my %xref_to_objects = %$xref_ids_hashref;

  my $dir = $self->core->dir();

  open (XREF, ">$dir/xref.txt");
  open (OBJECT_XREF, ">>$dir/object_xref.txt");
  open (EXTERNAL_SYNONYM, ">$dir/external_synonym.txt");
  open (GO_XREF, ">>$dir/go_xref.txt");

  # Cache synonyms for later use
  # Do one big query to get a list of all the synonyms; note each xref may have
  # more than one synonym so they are stored in a hash of lists
  my $syn_count;
  my %synonyms;
  my $syn_sth = $self->xref->dbc->prepare("SELECT xref_id, synonym FROM synonym");
  $syn_sth->execute();

  my ($sxref_id, $synonym);
  $syn_sth->bind_columns(\$sxref_id, \$synonym);
  while ($syn_sth->fetch()) {

    push @{$synonyms{$sxref_id}}, $synonym;

  }

  # keep a unique list of source IDs to build the external_db table later
  my %source_ids;

  my $object_xref_id = $start_object_xref_id;

  # build cache of source id -> external_db id; note %source_to_external_db is global
  %source_to_external_db = $self->map_source_to_external_db();

  # execute several queries with a max of 200 entries in each IN clause - more efficient
  my $batch_size = 200;

  # keep track of what xref_id & object_xref_ids have been written to prevent
  # duplicates; e.g. several dependent xrefs may be dependent on the same master xref.
  # Note %xrefs_written and %object_xrefs_written are global

  while(@xref_ids) {

    my @ids;
    if($#xref_ids > $batch_size) {
      @ids = splice(@xref_ids, 0, $batch_size);
    } else {
      @ids = splice(@xref_ids, 0);
    }

    my $id_str;
    if(@ids > 1)  {
      $id_str = "IN (" . join(',', @ids). ")";
    } else {
      $id_str = "= " . $ids[0];
    }


    my $sql = "SELECT * FROM xref WHERE xref_id $id_str";
    my $xref_sth = $self->xref->dbc->prepare($sql);
    $xref_sth->execute();

    my ($xref_id, $accession, $version, $label, $description, $source_id, $species_id, $master_xref_id, $linkage_annotation);
    $xref_sth->bind_columns(\$xref_id, \$accession, \$version, \$label, \$description, \$source_id, \$species_id);

    # note the xref_id we write to the file is NOT the one we've just read
    # from the internal xref database as the ID may already exist in the
    # core database so we add on $xref_id_offset
    while ($xref_sth->fetch()) {

      # make sure label is set to /something/ so that the website displays something
      $label = $accession if (!$label);

      if (!$xrefs_written{$xref_id}) {
	my $external_db_id = $source_to_external_db{$source_id};
	if ($external_db_id) { # skip "unknown" sources
	  if(!defined($updated_source{$external_db_id})){
	    $self->cleanup_sources_file($external_db_id);
	  }
	  print XREF ($xref_id+$xref_id_offset) . "\t" . $external_db_id . "\t" . $accession . "\t" . $label . "\t" . $version . "\t" . $description . "\n";
	  $xrefs_written{$xref_id} = 1;
	  $source_ids{$source_id} = $source_id;
	}
      }
    }

    # Now get the dependent xrefs for each of these xrefs and write them as well
    # Store the go_linkage_annotations as we go along (need for dumping go_xref)
    my $go_source_id = get_source_id_from_source_name($self->xref, "GO");

    $sql = "SELECT DISTINCT(x.xref_id), dx.master_xref_id, x.accession, x.label, x.description, x.source_id, x.version, dx.linkage_annotation FROM dependent_xref dx, xref x WHERE x.xref_id=dx.dependent_xref_id AND master_xref_id $id_str";

    my $dep_sth = $self->xref->dbc->prepare($sql);
    $dep_sth->execute();

    $dep_sth->bind_columns(\$xref_id, \$master_xref_id, \$accession, \$label, \$description, \$source_id, \$version, \$linkage_annotation);
    while ($dep_sth->fetch()) {


      my $external_db_id = $source_to_external_db{$source_id};
      next if (!$external_db_id);


      $label = $accession if (!$label);

      if (!$xrefs_written{$xref_id}) {
	if(!defined($updated_source{$external_db_id})){
	  $self->cleanup_sources_file($external_db_id);
	}
	
	print XREF ($xref_id+$xref_id_offset) . "\t" . $external_db_id . "\t" . $accession . "\t" . $label . "\t" . $version . "\t" . $description . "\tDEPENDENT\n";
	$xrefs_written{$xref_id} = 1;
	$source_ids{$source_id} = $source_id;
      }

      # create an object_xref linking this (dependent) xref with any objects it maps to
      # write to file and add to object_xref_mappings
      if (defined $xref_to_objects{$master_xref_id}) {
	my @ensembl_object_ids = keys( %{$xref_to_objects{$master_xref_id}} ); 
	foreach my $object_id_key (@ensembl_object_ids) {
	  my ($object_id, $type) = split /\|/, $object_id_key;
	  my $full_key = $type."|".$object_id."|".$xref_id;
	  if (!$object_xrefs_written{$full_key}) {
	    print OBJECT_XREF "$object_xref_id\t$object_id\t$type\t" . ($xref_id+$xref_id_offset) . "\tDEPENDENT\n";
	    # Add this mapping to the list - note NON-OFFSET xref_id is used
	    my $key = $type . "|" . $object_id;
	    push @{$object_xref_mappings{$key}}, $xref_id;
	    $object_xrefs_written{$full_key} = 1;

	    # Also store *parent's* query/target identity for dependent xrefs
	    print GO_XREF $object_xref_id . "\t" . $linkage_annotation . "\n"  if ($source_id == $go_source_id);
	    $object_xref_identities{$key}->{$xref_id}->{"target_identity"} = $object_xref_identities{$key}->{$master_xref_id}->{"target_identity"};
	    $object_xref_identities{$key}->{$xref_id}->{"query_identity"} = $object_xref_identities{$key}->{$master_xref_id}->{"query_identity"};

	    # write a go_xref with the appropriate linkage type
	    $object_xref_id++;

	  }
	}
      }
    }

    #print "source_ids: " . join(" ", keys(%source_ids)) . "\n";

  } # while @xref_ids

  # Dump any synonyms for xrefs we've written
  # Now write the synonyms we want to the file
  foreach my $xref_id (keys %synonyms) {
    foreach my $syn (@{$synonyms{$xref_id}}) {
      print EXTERNAL_SYNONYM ($xref_id+$xref_id_offset) . "\t" . $syn . "\n";
      $syn_count++;
    }
  }

  print "Wrote $syn_count synonyms\n";

  close(XREF);
  close(OBJECT_XREF);
  close(EXTERNAL_SYNONYM);
  close(GO_XREF);

  # calculate display_xref_ids for transcripts and genes
  my $transcript_display_xrefs = $self->build_transcript_display_xrefs($xref_id_offset);

  $self->build_genes_to_transcripts();

  $self->build_gene_display_xrefs($transcript_display_xrefs);

  # now build gene descriptions
  $self->build_gene_descriptions();

  return $object_xref_id;

}


# produce output for comparison with existing ensembl mappings
# format is (with header)
# xref_accession ensembl_type ensembl_id

sub dump_comparison {

  my $self = shift;

  my $dir = $self->core->dir();

  print "Dumping comparison data\n";

  open (COMPARISON, ">comparison/xref_mappings.txt");
  print COMPARISON "xref_accession" . "\t" . "ensembl_type" . "\t" . "ensembl_id\n";

  # get the xref accession for each xref as the xref_ids are ephemeral
  # first read all the xrefs that were dumped and get an xref_id->accession map
  my %xref_id_to_accesson;
  open (XREF, "$dir/xref.txt");
  while (<XREF>) {
    my ($xref_id,$external_db_id,$accession,$label,$version,$description) = split;
    $xref_id_to_accesson{$xref_id} = $accession;
  }
  close (XREF);

  open (OBJECT_XREF, "$dir/object_xref.txt");
  while (<OBJECT_XREF>) {
    my ($object_xref_id,$object_id,$type,$xref_id) = split;
    print COMPARISON $xref_id_to_accesson{$xref_id} . "\t" . $type . "\t" . $object_id . "\n";
  }

  close (OBJECT_XREF);
  close (COMPARISON);

}

sub build_transcript_display_xrefs {

  my ($self, $xref_id_offset) = @_;

  my $dir = $self->core->dir();

  # get a list of xref sources; format:
  # key: xref_id value: source_name
  # lots of these; if memory is a problem, just get the source ID (not the name)
  # and look it up elsewhere
  # note %xref_to_source is global
  print "Building xref->source mapping table\n";
  my $sql = "SELECT x.xref_id, s.name FROM source s, xref x WHERE x.source_id=s.source_id";
  my $sth = $self->xref->dbc->prepare($sql);
  $sth->execute();

  my ($xref_id, $source_name);
  $sth->bind_columns(\$xref_id, \$source_name);

  while ($sth->fetch()) {
    $xref_to_source{$xref_id} = $source_name;
  }

  print "Got " . scalar(keys %xref_to_source) . " xref-source mappings\n";

  # Cache the list of translation->transcript mappings & vice versa
  # Note variables are global
  print "Building translation to transcript mappings\n";
  $sth = $self->core->dbc->prepare("SELECT translation_id, transcript_id FROM translation");
  $sth->execute();

  my ($translation_id, $transcript_id);
  $sth->bind_columns(\$translation_id, \$transcript_id);

  while ($sth->fetch()) {
    $translation_to_transcript{$translation_id} = $transcript_id;
    $transcript_to_translation{$transcript_id} = $translation_id if ($translation_id);
  }

  print "Building transcript display_xrefs\n";
  my @priorities = $self->transcript_display_xref_sources();

  my $n = 0;

  # go through each object/xref mapping and store the best ones as we go along
  my %obj_to_best_xref;

  foreach my $key (keys %object_xref_mappings) {

    my ($type, $object_id) = split /\|/, $key;

    next if ($type !~ /(Transcript|Translation)/i);

    # if a transcript has more than one associated xref,
    # use the one with the highest priority, i.e. lower list position in @priorities
    my @xrefs = @{$object_xref_mappings{$key}};
    my ($best_xref, $best_xref_priority_idx);
    # store best query & target identities for each source for this object
    my %best_qi;
    my %best_ti;
    $best_xref_priority_idx = 99999;
    foreach my $xref (@xrefs) {

      my $source = $xref_to_source{$xref};
      if ($source) {
	my $i = find_in_list($source, @priorities);

	my $s = $source . "|" . $object_id;
	my $key = $type . "|" . $object_id;
	my $query_identity  = $object_xref_identities{$key}->{$xref}->{"query_identity"};
	my $target_identity = $object_xref_identities{$key}->{$xref}->{"target_identity"};

	# Check if this source has a better priority than the current best one
	# Note if 2 sources are the same priority, the mappings are compared on
	# query_identity then target_identity
	if ($i > -1 && $i <= $best_xref_priority_idx &&
	    (($query_identity > $best_qi{$s}) ||
	    ($query_identity == $best_qi{$s} && $target_identity > $best_ti{$s}))) {
	#if ($i > -1 && $i <= $best_xref_priority_idx && $query_identity > $best_qi{$s}) {
	  $best_xref = $xref;
	  $best_xref_priority_idx = $i;
	  $best_qi{$s} = $query_identity;
	  $best_ti{$s} = $target_identity;
	}
      } else {
	warn("Couldn't find a source for xref id $xref " . $xref_accessions{$xref_id} . "\n");
      }
    }
    # store object type, id, and best xref id and source priority
    if ($best_xref) {
      $obj_to_best_xref{$key} = $best_xref . "|" . $best_xref_priority_idx;
    }

  }


  # Now go through each of the calculated best xrefs and convert any that are
  # calculated against translations to be associated with their transcript,
  # if the priority of the translation xref is higher than that of the transcript
  # xref.
  # Needs to be done this way to avoid clobbering higher-priority transcripts.

  # hash keyed on transcript id, value is xref_id|source prioirity index
  my %transcript_display_xrefs;

  # Write a .sql file that can be executed, and a .txt file that can be processed
  open (TRANSCRIPT_DX, ">$dir/transcript_display_xref.sql");
  open (TRANSCRIPT_DX_TXT, ">$dir/transcript_display_xref.txt");

  foreach my $key (keys %obj_to_best_xref) {

    my ($type, $object_id) = split /\|/, $key;

    my ($best_xref, $best_xref_priority_idx) = split /\|/, $obj_to_best_xref{$key};

    # If transcript has a translation, use the best xref out of the transcript & translation

    my $transcript_id;
    my $translation_id;
    if ($type =~ /Transcript/i) {
      $transcript_id = $object_id;
      $translation_id = $transcript_to_translation{$transcript_id};
    }
    elsif ($type =~ /Translation/i) {
      $translation_id = $object_id;
      $transcript_id = $translation_to_transcript{$translation_id};
      $object_id = $transcript_id;
    }
    else{
      print "Cannot deal with type $type\n";
      next;
    }
    if ($translation_id) {
      my ($translation_xref, $translation_priority) = split /\|/, $obj_to_best_xref{"Translation|$translation_id"};
      my ($transcript_xref, $transcript_priority)   = split /\|/, $obj_to_best_xref{"Transcript|$transcript_id"};
      my $transcript_qi = $object_xref_identities{"Transcript|$object_id"}->{$transcript_xref}->{"query_identity"};
      my $translation_qi = $object_xref_identities{"Translation|$object_id"}->{$translation_xref}->{"query_identity"};

      #print "transcript xref $transcript_xref transcript pri $transcript_priority transcript qi: $transcript_qi translation xref $translation_xref translation pri $translation_priority translation qi $translation_qi\n";
      if(!$translation_xref){
	$best_xref = $transcript_xref;
	$best_xref_priority_idx = $transcript_priority;
      }
      if(!$transcript_xref){
	$best_xref = $translation_xref;
	$best_xref_priority_idx = $translation_priority;
      }
      elsif ($translation_priority < $transcript_priority) {
	$best_xref = $translation_xref;
	$best_xref_priority_idx = $translation_priority;
      } else {
	$best_xref = $transcript_xref;
	$best_xref_priority_idx = $transcript_priority;
      }

    }
    if ($best_xref) {

      # Write record with xref_id_offset
      print TRANSCRIPT_DX "UPDATE transcript SET display_xref_id=" . ($best_xref+$xref_id_offset) . " WHERE transcript_id=" . $object_id . ";\n";
      print TRANSCRIPT_DX_TXT ($best_xref+$xref_id_offset) . "\t" . $object_id . "\n";
      $n++;

      my $value = ($best_xref+$xref_id_offset) . "|" . $best_xref_priority_idx;
      $transcript_display_xrefs{$object_id} = $value;

    }

  }

  close(TRANSCRIPT_DX);
  close(TRANSCRIPT_DX_TXT);

  print "Wrote $n transcript display_xref entries to transcript_display_xref.sql\n";

  return \%transcript_display_xrefs;

}


# Assign display_xrefs to genes based on transcripts
# Gene gets the display xref of the highest priority of all of its transcripts
# If more than one transcript with the same priority, longer transcript is used

sub build_gene_display_xrefs {

  my ($self, $transcript_display_xrefs) = @_;

  my $dir = $self->core->dir();

  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-dbconn => $self->core->dbc);

  my $ta = $db->get_TranscriptAdaptor();

  print "Assigning display_xrefs to genes\n";

  open (GENE_DX, ">$dir/gene_display_xref.sql");
  open (GENE_DX_TXT, ">$dir/gene_display_xref.txt");
  my $hit = 0;
  my $miss = 0;
  my $trans_no_xref = 0;
  my $trans_xref = 0;
  foreach my $gene_id (keys %genes_to_transcripts) {

    my @transcripts = @{$genes_to_transcripts{$gene_id}};

    my $best_xref=undef;
    my $best_xref_priority_idx = 99999;
    my $best_transcript_length = -1;
    foreach my $transcript_id (@transcripts) {
      if (!$transcript_display_xrefs->{$transcript_id}) {
	$trans_no_xref++;
	next;
      } else {
	$trans_xref++;
      }
      my ($xref_id, $priority) = split (/\|/, $transcript_display_xrefs->{$transcript_id});

      # 2 separate if clauses to avoid having to fetch transcripts unnecessarily

      if (($priority < $best_xref_priority_idx)) {

	$best_xref_priority_idx = $priority;
	$best_xref = $xref_id;

      } elsif ($priority == $best_xref_priority_idx) {

	# compare transcript lengths and use longest
	my $transcript = $ta->fetch_by_dbID($transcript_id);
	my $transcript_length = $transcript->length();
	if ($transcript_length > $best_transcript_length) {
	  $best_transcript_length = $transcript_length;
	  $best_xref_priority_idx = $priority;
	  $best_xref = $xref_id;
	}
      }
    }

    if (defined($best_xref)) {
      # Write record
      print GENE_DX "UPDATE gene g, analysis a SET g.display_xref_id=" . $best_xref . " WHERE g.gene_id=" . $gene_id . " AND g.analysis_id=a.analysis_id AND a.logic_name != \"ncRNA\";\n";
      print GENE_DX_TXT $best_xref . "\t" . $gene_id ."\n";
      $hit++;
    } else {
      $miss++;
    }

  }

  close (GENE_DX);
  close (GENE_DX_TXT);
  print "Transcripts with no xrefs: $trans_no_xref with xrefs: $trans_xref\n";
  print "Wrote $hit gene display_xref entries to gene_display_xref.sql\n";
  print "Couldn't find display_xrefs for $miss genes\n" if ($miss > 0);
  print "Found display_xrefs for all genes\n" if ($miss == 0);

  return \%genes_to_transcripts;

}

# Display xref sources to be used for transcripts *in order of priority*
# Source names used must be identical to those in the source table.

sub transcript_display_xref_sources {

  return ('HUGO',
	  'MarkerSymbol',
#	  'wormbase_transcript',
	  'flybase_symbol',
	  'Anopheles_symbol',
	  'Genoscope_annotated_gene',
	  'Genoscope_predicted_transcript',
	  'Genoscope_predicted_gene',
	  'Uniprot/SWISSPROT',
	  'RefSeq_peptide',
	  'RefSeq_dna',
	  'Uniprot/SPTREMBL',
	  'RefSeq_peptide_predicted',
	  'RefSeq_dna_predicted',
	  'EntrezGene');

}

# Get transcripts associated with each gene

sub build_genes_to_transcripts {

  my ($self) = @_;

  print "Getting transcripts for all genes\n";

  my $sql = "SELECT gene_id, transcript_id FROM transcript";
  my $sth = $self->core->dbc->prepare($sql);
  $sth->execute();

  my ($gene_id, $transcript_id);
  $sth->bind_columns(\$gene_id, \$transcript_id);

  # Note %genes_to_transcripts is global
  while ($sth->fetch()) {
    push @{$genes_to_transcripts{$gene_id}}, $transcript_id;
  }

  print "Got " . scalar keys(%genes_to_transcripts) . " genes\n";

}

# Find the index of an item in a list(ref), or -1 if it's not in the list.
# Only look for exact matches (case insensitive)

sub find_in_list {

  my ($item, @list) = @_;

  for (my $i = 0; $i < scalar(@list); $i++) {
    if (lc($list[$i]) eq lc($item)) {
      return $i;
    }
  }

  return -1;

}

# Take a string and a list of regular expressions
# Find the index of the highest matching regular expression
# Return the index, or -1 if not found.

sub find_match {

 my ($str, @list) = @_;

 my $str2 = $str;
 my $highest_index = -1;

  for (my $i = 0; $i < scalar(@list); $i++) {
    my $re = $list[$i];
    if ($str2 =~ /$re/i) {
      $highest_index = $i;
    }
  }

  return $highest_index;

}

# Build a map of source id (in xref database) to external_db (in core database)

sub map_source_to_external_db {

  my $self = shift;

  my %source_to_external_db;

  # get all sources
  my $sth = $self->xref->dbc->prepare("SELECT source_id, name FROM source");
  $sth->execute();
  my ($source_id, $source_name);
  $sth->bind_columns(\$source_id, \$source_name);

  while($sth->fetchrow_array()) {

    # find appropriate external_db_id for each one
    my $sql = "SELECT external_db_id FROM external_db WHERE db_name=?";
    my $core_sth = $self->core->dbc->prepare($sql);
    $core_sth->execute($source_name);

    my @row = $core_sth->fetchrow_array();

    if (@row) {

      $source_to_external_db{$source_id} = $row[0];
      #print "Source name $source_name id $source_id corresponds to core external_db_id " . $row[0] . "\n";

    } else {

      print STDERR "Can't find external_db entry for source name $source_name; xrefs for this source will not be written. Consider adding $source_name to external_db\n"

    }

  } # while source

  return %source_to_external_db;
}


sub cleanup_sources_file{
  my ($self,$id) = @_;

  $updated_source{$id} =1;

  my $dir = $self->core->dir();
  open (DEL, ">>$dir/cleanup.sql") || die "Could not open $dir/cleanup.sql\n";


  print DEL "DELETE external_synonym ";
  print DEL     "FROM external_synonym, xref ";
  print DEL       "WHERE external_synonym.xref_id = xref.xref_id ";
  print DEL         "AND xref.external_db_id = $id\n";


  print DEL "DELETE identity_xref ";
  print DEL     "FROM identity_xref, object_xref, xref ";
  print DEL       "WHERE identity_xref.object_xref_id = object_xref.object_xref_id ";
  print DEL         "AND object_xref.xref_id = xref.xref_id ";
  print DEL         "AND xref.external_db_id = $id \n";


  print DEL "DELETE object_xref ";
  print DEL     "FROM object_xref, xref ";
  print DEL       "WHERE object_xref.xref_id = xref.xref_id ";
  print DEL         "AND xref.external_db_id = $id\n";


  print DEL "DELETE FROM xref WHERE xref.external_db_id = $id \n";

  close DEL;

}

# Upload .txt files and execute .sql files.

sub do_upload {

  my ($self) = @_;

  my $ensembl = $self->core;
  my $core_db = $ensembl->dbc;
  # xref.txt etc

  print "Deleting existing data\n";
  my $file = $ensembl->dir() . "/cleanup.sql";
  open(CLEAN,"<$file") || die "could not open $file for reading \n";
  while(<CLEAN>){
    chomp;
    my $sth = $core_db->prepare($_);
    $sth->execute() or die "Couldn't execute statement: " . $sth->errstr;
  }     
  close CLEAN;
  
  
  foreach my $table ("go_xref", "interpro") {
    
    my $sth = $core_db->prepare("DELETE FROM $table");
    print "Deleting existing data in $table\n";
    $sth->execute();
    
  }
  
  # gene & transcript display_xrefs
  my $sth = $core_db->prepare(<<GADES);
  UPDATE gene g, analysis a 
    SET g.display_xref_id=NULL 
      WHERE g.analysis_id=a.analysis_id 
	AND a.logic_name != "ncRNA"
GADES
  print "Setting all existing display_xref_id in gene to null\n";
  $sth->execute();
  

  my $sth = $core_db->prepare(<<TRAN);
  UPDATE transcript t, gene g, analysis a 
    SET t.display_xref_id=NULL
      WHERE g.analysis_id = a.analysis_id
	AND a.logic_name != "ncRNA"
	  AND g.gene_id =t.gene_id
TRAN
  print "Setting all existing display_xref_id in transcript to null\n";
  $sth->execute();

  # gene descriptions
  my $sth = $core_db->prepare(<<GENE);
  UPDATE gene g, analysis a 
    SET g.description=NULL 
      WHERE g.analysis_id=a.analysis_id 
	AND a.logic_name != "ncRNA";
GENE
  print "Setting all existing descriptions in gene table to null\n";
  $sth->execute();

  print "Uploading new data\n";
  foreach my $table ("xref", "object_xref", "identity_xref", "external_synonym", "go_xref", "interpro") {

    my $file = $ensembl->dir() . "/" . $table . ".txt";

    # don't seem to be able to use prepared statements here
    my $sth = $core_db->prepare("LOAD DATA INFILE \'$file\' IGNORE INTO TABLE $table");
    print "Uploading data in $file to $table\n";
    $sth->execute();

  }

  # gene_display_xref.sql etc
  foreach my $table ("gene", "transcript") {

    my $file = $ensembl->dir() . "/" . $table . "_display_xref.sql";

    print "Setting $table display_xrefs from $file\n";
    my $str = "mysql -u " .$core_db->username() ." -p" . $core_db->password() . " -h " . $core_db->host() ." -P " . $core_db->port() . " " .$core_db->dbname() . " < $file";
    system $str;

  }

  # gene descriptions
  my $file = $ensembl->dir() . "/gene_description.sql";
  print "Setting gene descriptions from $file\n";
  my $str = "mysql -u " .$core_db->username() ." -p" . $core_db->password() . " -h " . $core_db->host() ." -P " . $core_db->port() . " " .$core_db->dbname() . " < $file";
  system $str;

  # update meta table with timestamp saying when xrefs were last updated
  $file =  $ensembl->dir() . "/meta_timestamp.sql";
  open (FILE, ">$file");
  print FILE "DELETE FROM meta WHERE meta_key='xref.timestamp';\n";
  print FILE "INSERT INTO meta (meta_key,meta_value) VALUES ('xref.timestamp', NOW())\n";
  close(FILE);

  my $str = "mysql -u " .$core_db->username() ." -p" . $core_db->password() . " -h " . $core_db->host() ." -P " . $core_db->port() . " " .$core_db->dbname() . " < $file";
  system $str;

  # set gene and transcript statuses to KNOWN or NOVEL based on external_db status
  print "Setting gene and transcript status from external_db KNOWN/NOVEL\n";
  $file =  $ensembl->dir() . "/gene_transcript_status.sql";
  open (FILE, ">$file");
  print FILE "UPDATE transcript SET status=\'NOVEL\';\n";
  print FILE "UPDATE gene SET status=\'NOVEL\';\n";
  print FILE "UPDATE gene g, xref x, external_db e SET g.status = \'KNOWN\' WHERE g.display_xref_id = x.xref_id AND x.external_db_id = e.external_db_id AND e.status=\'KNOWN\';\n";
  print FILE "UPDATE transcript t, xref x, external_db e SET t.status = \'KNOWN\' WHERE t.display_xref_id = x.xref_id AND x.external_db_id = e.external_db_id AND e.status=\'KNOWN\';\n";
  close(FILE);

  my $str = "mysql -u " .$core_db->username() ." -p" . $core_db->password() . " -h " . $core_db->host() ." -P " . $core_db->port() . " " .$core_db->dbname() . " < $file";
  system $str;

}

# Assign gene descriptions
# Algorithm:
# foreach gene
#   get all transcripts & translations
#   get all associated xrefs
#   filter by regexp, discard blank ones
#   if, after filtering, the description is not blank, use the *original* description
#   order by source & keyword
#   assign description of best xref to gene
# }
#
# One gene may have several associated peptides; the one to use is decided as follows.
# In decreasing order of precedence:
#
# - Consortium xref, e.g. ZFIN for zebrafish
#
# - UniProt/SWISSPROT
#     If there are several, the one with the best %query_id then %target_id is used
#
# - RefSeq
#    If there are several, the one with the best %query_id then %target_id is used
#
# - UniProt/SPTREMBL
#    If there are several, precedence is established on the basis of the occurrence of 
#    regular expression patterns in the description.

sub build_gene_descriptions {

  my ($self) = @_;

  # TODO - don't call this from, but after, gene_display_xref

  # Get all xref descriptions, filtered by regexp.
  # Discard any that are blank (i.e. regexp has removed everything)

  print "Getting & filtering xref descriptions\n";
 # Note %xref_descriptions & %xref_accessions are global

  my $sth = $self->xref->dbc->prepare("SELECT xref_id, accession, description FROM xref");
  $sth->execute();
  my ($xref_id, $accession, $description);
  $sth->bind_columns(\$xref_id, \$accession, \$description);

  my $removed = 0;
  my @regexps = $self->gene_description_filter_regexps();
  while ($sth->fetch()) {
    if ($description) {
      my $filtered_description = $self->filter_by_regexp($description, \@regexps);
      if ($filtered_description ne "") {
	$xref_descriptions{$xref_id} = $description;
	$xref_accessions{$xref_id} = $accession;
      } else {
	$removed++;
      }
    }
  }

  print "Regexp filtering (" . scalar(@regexps) . " regexps) removed $removed descriptions, left with " . scalar(keys %xref_descriptions) . "\n";

  my $dir = $self->core->dir();
  open(GENE_DESCRIPTIONS,">$dir/gene_description.sql") || die "Could not open $dir/gene_description.sql";

  # Foreach gene, get any xrefs associated with its transcripts or translations

  # just in case species specific routines missed the population of the
  # needed hashes. Check an populate if miising.
  if(!scalar(keys %translation_to_transcript)){
    print "Building translation to transcript mappings\n";
    $sth = $self->core->dbc->prepare("SELECT translation_id, transcript_id FROM translation");
    $sth->execute();
    
    my ($translation_id, $transcript_id);
    $sth->bind_columns(\$translation_id, \$transcript_id);
    
    while ($sth->fetch()) {
      $translation_to_transcript{$translation_id} = $transcript_id;
      $transcript_to_translation{$transcript_id} = $translation_id if ($translation_id);
    }
  }
  if(!scalar(keys %xref_to_source)){
    print "Building xref->source mapping table\n";
    my $sql = "SELECT x.xref_id, s.name FROM source s, xref x WHERE x.source_id=s.source_id";
    my $sth = $self->xref->dbc->prepare($sql);
    $sth->execute();

    my ($xref_id, $source_name);
    $sth->bind_columns(\$xref_id, \$source_name);

    while ($sth->fetch()) {
      $xref_to_source{$xref_id} = $source_name;
    }

    print "Got " . scalar(keys %xref_to_source) . " xref-source mappings\n";

  }


  print "Assigning gene descriptions\n";


  foreach my $gene_id (keys %genes_to_transcripts) {

    my @gene_xrefs;

    my %local_xref_to_object;

    my @transcripts = @{$genes_to_transcripts{$gene_id}};
    foreach my $transcript (@transcripts) {

      my @xref_ids;

      my $key = "Transcript|$transcript";
      if ($object_xref_mappings{$key}) {

	@xref_ids = @{$object_xref_mappings{$key}};
	push @gene_xrefs, @xref_ids;
	foreach my $xref (@xref_ids) {
	  $local_xref_to_object{$xref} = $key;
	}

      }

      my $translation = $transcript_to_translation{$transcript};
      $key = "Translation|$translation";
      if ($object_xref_mappings{$key}) {

	@xref_ids = @{$object_xref_mappings{$key}};
	push @gene_xrefs, @xref_ids ;
	foreach my $xref (@xref_ids) {
	  $local_xref_to_object{$xref} = $key;
	}
      }

    }

    # Now sort through these and find the "best" description and write it

    if (@gene_xrefs) {

#      @gene_xrefs = @{$self->strip(\@gene_xrefs)};  # remove blank entrys

      @gene_xrefs = sort {compare_xref_descriptions($self->consortium(), $gene_id, \%local_xref_to_object)} @gene_xrefs;

#      my $best_xref = $gene_xrefs[-1];
      my $best_xref = $self->get_best(\@gene_xrefs,$gene_id);
      my $source = $xref_to_source{$best_xref};

      # only store the description if its source is one of the allowed ones
      if (find_in_list($source, $self->gene_description_sources()) > -1) {

	my $description = $xref_descriptions{$best_xref};
	my $acc = $xref_accessions{$best_xref};

	$description =~ s/\"//ig; # remove " as they will cause problems in .sql files

	my $desc = $description . " [Source:$source;Acc:$acc]";

	# prevent overwriting ncRNA gene descriptions as these are calculated by an external method
	print GENE_DESCRIPTIONS "UPDATE gene g, analysis a SET g.description=\"$desc\" WHERE a.analysis_id=g.analysis_id AND a.logic_name != \"ncRNA\" AND g.gene_id=$gene_id;\n" if ($description);

      }

    }

  } # foreach gene

  close(GENE_DESCRIPTIONS);

}
# Check if any .err files exist that have non-zero size;
# this indicates that something has gone wrong with the exonerate run

sub check_err {

  my ($self, $dir) = @_;

  foreach my $err (glob("$dir/*.err")) {

    print "\n\n*** Warning: $err has non-zero size; may indicate problems with exonerate run\n\n\n" if (-s $err);

  }
}

sub get_best {
  my ($self,$refxref,$gene_id) = @_;
  return $$refxref[-1];
}

sub strip{
  my ($self,$refarray) = @_;
  my @ret;

  foreach my $arr (@$refarray){
    if(defined($xref_accessions{$arr})){
      push @ret, $arr;
    }
  }
  return \@ret;
}

# remove a list of patterns from a string
sub filter_by_regexp {

  my ($self, $str, $regexps) = @_;

  foreach my $regexp (@$regexps) {
    $str =~ s/$regexp//ig;
  }

  return $str;

}

# Regexp used for filter out useless text from gene descriptions
# Method can be overridden in species-specific modules
sub gene_description_filter_regexps {

  return ();

}


# The "consortium" source for this species, should be the same as in
# source table

sub consortium {

  return "xxx"; # Default to something that won't be matched as a source

}



# move translation hits on to their transcripts if the same external_db_id
# has hits to both translation and transcript


sub cleanup_database {
  my ($self) = @_;
  my $ensembl_dbc = $self->core->dbc;

  # 1 make sure external_db do not lie on both translations and transcripts

  my $sql_check  = (<<ESQL);
  SELECT x.external_db_id, ox.ensembl_object_type, COUNT(*), e.db_name
    FROM xref x, object_xref ox, external_db e 
     WHERE x.xref_id = ox.xref_id AND e.external_db_id = x.external_db_id
       GROUP BY x.external_db_id, ox.ensembl_object_type
ESQL

  my $sth = $ensembl_dbc->prepare($sql_check);
  $sth->execute();
  my $previous_id = -1;
  my $previous_type ="";
  while(my @row = $sth->fetchrow_array()){
    my $external_db_id = $row[0];
    if($external_db_id == $previous_id){
      $self->fix_mart_prob($row[3],$external_db_id,$row[1],$previous_type);  
    }
    $previous_id = $external_db_id;
    $previous_type = $row[1];
  }

  # now recheck just incase :-)
  $sth->execute();
  $previous_id = -1;
  $previous_type ="";
  my $error =0;
  while(my @row = $sth->fetchrow_array()){
    my $external_db_id = $row[0];
    if($external_db_id == $previous_id){
      print "Problem: Still have multiple associations with ".$row[3]."\n";
      $error++;
    }
    $previous_id = $external_db_id;
    $previous_type = $row[1];
  }
  $sth->finish();
  if(!$error){
    print "External databases only associate to one ensembl type (PASS)\n";
  }
  else{
    print "External databases only associate to one ensembl type (FAIL)\n";
  }
}


sub fix_mart_prob{
  my ($self,$db_name,$db_id,$type1,$type2) = @_;
  my $ensembl_dbc = $self->core->dbc;

  print "$db_name is associated with both $type1 and $type2 object types\n";

  my $to;
  if($type1 eq "Gene" or $type2 eq "Gene"){
    $to = "Gene";
  }
  else{
    $to = "Transcript";
  }
    
  print "Therefore moving all associations to the ".$to."s\n";


  $ensembl_dbc->do("CREATE TABLE object_xref2 like object_xref");

  $ensembl_dbc->do("ALTER TABLE object_xref DROP INDEX ensembl_object_type");

# Move translations onto the transcripts
  my $sql =(<<EOF);
  UPDATE object_xref, translation, xref
     SET object_xref.ensembl_object_type = 'Transcript',
         object_xref.ensembl_id = translation.transcript_id 
     WHERE object_xref.ensembl_object_type = 'Translation' AND
           object_xref.ensembl_id = translation.translation_id AND
           xref.xref_id = object_xref.xref_id AND
           xref.external_db_id = $db_id;
EOF
  $ensembl_dbc->do($sql);
  
  if($to eq "Gene"){ #move transcripts to the gene
    my $sql =(<<GENE);
  UPDATE object_xref, transcript, xref
     SET object_xref.ensembl_object_type = 'Gene',
         object_xref.ensembl_id = transcript.gene_id 
     WHERE object_xref.ensembl_object_type = 'Transcript' AND
           object_xref.ensembl_id = transcript.transcript_id AND
           xref.xref_id = object_xref.xref_id AND
           xref.external_db_id = $db_id;
GENE
    $ensembl_dbc->do($sql);
    
  }
  
  #  print $sql."\n";
  
  #  $ensembl_dbc->do($sql);
  
  $ensembl_dbc->do("INSERT IGNORE INTO object_xref2 SELECT * FROM object_xref");
  
  $ensembl_dbc->do("DROP TABLE object_xref");

  $ensembl_dbc->do("ALTER TABLE object_xref2 RENAME object_xref");


}


# Sort a list of xrefs by the priority of their sources
# Assumed this function is called by Perl sort, passed with parameter
# See comment for build_gene_descriptions for how precedence is decided.

sub compare_xref_descriptions {

  my ($consortium, $gene_id, $xref_to_object) = @_;

  my @sources = gene_description_sources();
  push @sources, $consortium;

  my @words = qw(unknown hypothetical putative novel probable [0-9]{3} kDa fragment cdna protein);

  my $src_a = $xref_to_source{$a};
  my $src_b = $xref_to_source{$b};
  my $pos_a = find_in_list($src_a, @sources);
  my $pos_b = find_in_list($src_b, @sources);

  # If same source, need to do more work
  if ($pos_a == $pos_b) {

   if ($src_a eq "Uniprot/SWISSPROT" || $src_a =~ /RefSeq/) {

     # Compare on query identities, then target identities if queries are the same
     my $key_a = $xref_to_object->{$a}; # e.g. "Translation|1234"
     my $key_b = $xref_to_object->{$b};
     my ($type_a, $object_a) = split(/\|/, $key_a);
     my ($type_b, $object_b) = split(/\|/, $key_b);

     return 0 if ($type_a != $type_b); # only compare like with like

     my $query_identity_a = $object_xref_identities{$key_a}->{$a}->{"query_identity"};
     my $query_identity_b = $object_xref_identities{$key_b}->{$b}->{"query_identity"};

     #print "gene 78163 " . $xref_accessions{$a} . " key a $key_a qia $query_identity_a " . $xref_accessions{$b} . " key b $key_b qib $query_identity_b \n" if ($gene_id==78163);
     return ($query_identity_a <=> $query_identity_b) if ($query_identity_a != $query_identity_b);

     my $target_identity_a = $object_xref_identities{$key_a}->{$a}->{"target_identity"};
     my $target_identity_b = $object_xref_identities{$key_b}->{$b}->{"target_identity"};

     return ($target_identity_a <=> $target_identity_b);

   } elsif ($src_a eq "Uniprot/SPTREMBL") {

     # Compare on words
     my $wrd_idx_a = find_match($xref_descriptions{$a}, @words);
     my $wrd_idx_b = find_match($xref_descriptions{$b}, @words);
     return $wrd_idx_a <=> $wrd_idx_b;

   } else {

     return 0;

   }
    return 0;

  } else {

    return $pos_a <=> $pos_b;

  }
}

# list of sources to be used when building gene descriptions
# sorted into increasing order of priority

sub gene_description_sources {

  return ("RefSeq_dna_predicted",
	  "RefSeq_peptide_predicted",
	  "Uniprot/SPTREMBL",
	  "RefSeq_dna",
	  "RefSeq_peptide",
	  "Uniprot/SWISSPROT");

}

# load external_db (if it's empty) from ../external_db/external_dbs.txt

sub upload_external_db {
  my ($self) = @_;

  my $core_db = $self->core->dbc;
  $core_db->connect();
  my $row = @{$core_db->db_handle->selectall_arrayref("SELECT COUNT(*) FROM external_db")}[0];
  my $count = @{$row}[0];

  if ($count == 0) {
    my $edb = cwd() . "/../external_db/external_dbs.txt";
    print "external_db table is empty, uploading from $edb\n";
    my $edb_sth = $core_db->prepare("LOAD DATA INFILE \'$edb\' INTO TABLE external_db");
    $edb_sth->execute();
  } else {
    print "external_db table already has $count rows, will not change it\n";
   }

}

# Cache xref labels for debugging purposes

sub cache_xref_labels {

  my ($self,$xref_dbc) = @_;

  print "Caching xref labels\n";
  my $sth = $xref_dbc->prepare("SELECT xref_id, label FROM xref");
  $sth->execute();
  while(my @row = $sth->fetchrow_array()){
    $xref_labels{$row[0]} = $row[1];
  }

}
sub get_xref_descriptions{
  return \%xref_descriptions;
}

sub get_xref_accessions{
 return \%xref_accessions;
}

sub xref_id_offset{
 my  $self  = shift;

  $self->{'xref_id_offset'} = shift if( @_ );
  if( exists $self->{'xref_id_offset'} ) {
    return $self->{'xref_id_offset'};
  }
 return undef;
}

sub add_missing_pairs{
  my ($self) = @_;
  my $xref_id_offset = $self->xref_id_offset();
  #
  # add the pairs
  #
  # get current max object_xref_id
  my $row = @{$self->core->dbc->db_handle->selectall_arrayref("SELECT MAX(object_xref_id) FROM object_xref")}[0];
  my $max_object_xref_id = @{$row}[0];
  if (!defined $max_object_xref_id) {
    die ("No existing object_xref_ids, something very wrong\n");
  }
#  print "xref offset => $xref_id_offset\n";
#  print "max object xref => $max_object_xref_id \n";
  my $xref_sql = (<<EOS);
  SELECT x1.xref_id, x2.xref_id 
    FROM pairs p, xref x1, xref x2
      WHERE p.accession1 = x1.accession
	AND p.accession2 = x2.accession
	 AND p.source_id = x1.source_id
EOS
  my $xref_sth = $self->xref->dbc->prepare($xref_sql);
  $xref_sth->execute(); 
  my ($xref_id1,$xref_id2);
  $xref_sth->bind_columns(\$xref_id1, \$xref_id2);

  my %good2missed=();
  my $okay =0;
  my $both = 0;
  my $poss = 0;
  while ($xref_sth->fetch()) {
    if(!defined($xrefs_written{$xref_id1}) or !defined($xrefs_written{$xref_id2})){
      if (!defined($xrefs_written{$xref_id1}) and !defined($xrefs_written{$xref_id2})){
	$okay++;
      }
      elsif(!defined($xrefs_written{$xref_id2})){
	$poss++;
	$good2missed{$xref_id1+$xref_id_offset} = $xref_id2+$xref_id_offset;
      }
      else{
	$poss++;
	$good2missed{$xref_id2+$xref_id_offset} = $xref_id1+$xref_id_offset;
      }
    }
    else{
      $both++;
    }
  }
#  print "but apparently $okay have no matches at all and $both have two\n";
#  print "potential filler ins = $poss\n";
#  print "good2missed=>".scalar(%good2missed)."\n";
  open(OBJECT_XREF2, ">".$self->core->dir()."/pairs_object_xref.txt") || die "Could not open pairs_object_xref.txt";

  my $i=0;
  my $index;
  my $added = 0;
  my $sql = "SELECT xref_id, ensembl_id, ensembl_object_type FROM object_xref WHERE xref_id IN (";
  my ($goodxref, $ens_int_id,$type);
  my @list =();
  foreach my $key (keys %good2missed){
    if($i > 200){
      my $sth_ob = $self->core->dbc->prepare($sql.(join(',',@list)).")") || die @_;
      $sth_ob->execute();
      $sth_ob->bind_columns(\$goodxref,\$ens_int_id,\$type);
      while($sth_ob->fetch()){
	$max_object_xref_id++;
	$added++;
	print OBJECT_XREF2 "$max_object_xref_id\t$ens_int_id\t$type\t" .$good2missed{$goodxref} . "\tDEPENDENT\n";	
      }
      $sth_ob->finish();
      @list =();
      $i=0;
    }
    else{
      push @list, $key;
      $i++;
    }
  }
  if($i){
    my $sth_ob = $self->core->dbc->prepare($sql.(join(',',@list)).")") || die @_;
    $sth_ob->execute();
    $sth_ob->bind_columns(\$goodxref,\$ens_int_id,\$type);
    while($sth_ob->fetch()){
      $max_object_xref_id++;
      $added++;
      print OBJECT_XREF2 "$max_object_xref_id\t$ens_int_id\t$type\t" .$good2missed{$goodxref} . "\tDEPENDENT\n";	
    }
    $sth_ob->finish();
  }

  close OBJECT_XREF2;

  #
  # Now load the data into the database.
  #
 
  my $file = $self->core->dir()."/pairs_object_xref.txt";
  
  # don't seem to be able to use prepared statements here
  my $sth = $self->core->dbc->prepare("LOAD DATA INFILE \'$file\' IGNORE INTO TABLE object_xref");
  print "Uploading data in $file to object_xref\n";
  $sth->execute();
  
  print "$added new object xrefs added based on the Pairs\n";

}

sub dump_all_dependencies{
  my ($self, $master_id, $xref_id_offset) = @_;

  # Now get the dependent xrefs for this xref and write them
  
  my $sql = "SELECT DISTINCT(x.xref_id), dx.master_xref_id, x.accession, x.label, x.description, x.source_id, x.version, dx.linkage_annotation FROM dependent_xref dx, xref x WHERE x.xref_id=dx.dependent_xref_id AND master_xref_id = $master_id";
  
  my $dep_sth = $self->xref->dbc->prepare($sql);
  $dep_sth->execute();

   my ($xref_id, $accession, $version, $label, $description, $source_id, $species_id, $master_xref_id, $linkage_annotation);
  
  $dep_sth->bind_columns(\$xref_id, \$master_xref_id, \$accession, \$label, \$description, \$source_id, \$version, \$linkage_annotation);
  while ($dep_sth->fetch()) {
    
    
    my $external_db_id = $source_to_external_db{$source_id};
     next if (!$external_db_id);
    
    
    $label = $accession if (!$label);
    
    if (!$xrefs_written{$xref_id}) {
      if(!defined($updated_source{$external_db_id})){
	$self->cleanup_sources_file($external_db_id);
      }
      
      print XREF ($xref_id+$xref_id_offset) . "\t" . $external_db_id . "\t" . $accession . "\t" . $label . "\t" . $version . "\t" . $description . "\tDEPENDENT\n";
      $xrefs_written{$xref_id} = 1;
    }
  }
}
    


1;

