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

@ISA = qw{ XrefMapper::db };


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

  my ($self) = @_;

  $self->dump_xref();
  $self->dump_ensembl();

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
    push @dna, $self->ensembl_dna_file();
    push @list, \@dna;
    my @pep=();
    push @pep, $method;
    push @pep, $self->xref->dir."/xref_".$i."_peptide.fasta";
    push @pep, $self->ensembl_protein_file();
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
  my ($xref,$species) = @_;

  my $sql = "select species_id from species where name = '".$species."'";
  my $dbi = $xref->dbi();
  my $sth = $dbi->prepare($sql);
  $sth->execute();
  my @row = $sth->fetchrow_array();
  my $species_id;
  if (defined @row) {
    $species_id = $row[0];
  } else {
    print STDERR "Couldn't get ID for species ".$species."\n";
    print STDERR "It must be one of :-\n";
    $sql = "select name from species";
    $sth = $dbi->prepare($sql);
    $sth->execute();
    while(my @row = $sth->fetchrow_array()){
      print STDERR $row[0]."\n";
    }
    die("Please try again :-)\n");
  }
  $sth->finish();
  $dbi->disconnect();
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

  #  return [["ExonerateGappedBest1", ["homo_sapiens","Uniprot/SWISSPROT"]]];
  return [["ExonerateGappedBest1", ["homo_sapiens","*"]]];

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
  my ($xref, $source) = @_;
  my $source_id;
  
  my $sql = "select source_id from source where name = '".$source."'";
  my $dbi = $xref->dbi();
  my $sth = $dbi->prepare($sql);
  $sth->execute();
  my @row = $sth->fetchrow_array();
  if (defined $row[0] and $row[0] ne '') {
    $source_id = $row[0];
#    print $source."\t*".$row[0]."*\n";
  } else {
    print STDERR "Couldn't get ID for source ".$source."\n";
    print STDERR "It must be one of :-\n";
    $sql = "select name from source";
    $sth = $dbi->prepare($sql);
    $sth->execute();
    while(my @row = $sth->fetchrow_array()){
      print STDERR $row[0]."\n";
    }
    die("Please try again :-)\n");
  }  
  $sth->finish();
  $dbi->disconnect();
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
      $method[$i] = shift @$list;
      if(!-e $xref->dir()."/xref_".$i."_dna.fasta"){ 
	$skip = 0;
      }
      if(!-e $xref->dir()."/xref_".$i."_peptide.fasta"){ 
	$skip = 0;
      }
      $i++;
    }
    if($skip){
      $self->method(\@method);
      return;
    }
  }

  $i=0;
  foreach my $list (@lists){
#    print "method->".@$list[0]."\n";
    $method[$i] = shift @$list;
    my $j = 0;
    my @source_id=();
    my @species_id=();
    foreach my $element (@$list){
      while(my $species = shift(@$element)){
	#	print $j.")\t".$species."\n";
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
#	print $j."\t".$source. "\t".$source_id[$j] ."\n";
#	print $j."\t".$species."\t".$species_id[$j]."\n";
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

    my $sth = $xref->dbi()->prepare($sql);
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
  my ($self) = @_;

  $self->fetch_and_dump_seq();

}


=head2 fetch_and_dump_seq

  Description: Dumps the ensembl data to a file in fasta format.
  Returntype : none
  Exceptions : wil die if the are errors in db connection or file creation.
  Caller     : dump_ensembl

=cut

sub fetch_and_dump_seq{
  my ($self) = @_;

  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-species => $self->species(),
					      -dbname  => $self->dbname(),
					      -host    => $self->host(),
					      -port    => $self->port(),
					      -pass    => $self->password(),
					      -user    => $self->user(),
					      -group   => 'core');

  #
  # store ensembl dna file name and open it
  #
  
  # if no directory set then dump in the current directory.
  if(!defined($self->dir())){
    $self->dir(".");
  }
  $self->ensembl_dna_file($self->dir."/".$self->species."_dna.fasta");
  open(DNA,">".$self->ensembl_dna_file()) 
    || die("Could not open dna file for writing: ".$self->ensembl_dna_file."\n");

  #
  # store ensembl protein file name and open it
  #
  $self->ensembl_protein_file($self->dir."/".$self->species."_protein.fasta");

  if(defined($self->dumpcheck()) and -e $self->ensembl_protein_file() and -e $self->ensembl_dna_file()){
    return;
  }

  open(PEP,">".$self->ensembl_protein_file()) 
    || die("Could not open protein file for writing: ".$self->ensembl_protein_file."\n");

  my $gene_adap = $db->get_GeneAdaptor();
  my @gene_ids = @{$gene_adap->list_dbIDs()};
  my $max = undef;
  if(defined($self->maxdump())){
    $max = $self->maxdump();
  }
  my $i =0;
  foreach my $gene_id (@gene_ids){
    my $gene = $gene_adap->fetch_by_dbID($gene_id);
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



#=head2 xref_protein_file
# 
#  Arg [1]    : (optional) string $arg
#               the fasta file name for the protein xref
#  Example    : $file name = $xref->xref_protein_file();
#  Description: Getter / Setter for the protien xref fasta file 
#  Returntype : string
#  Exceptions : none
#
#=cut
#
#
#sub xref_protein_file{
#  my ($self, $arg) = @_;
#
#  (defined $arg) &&
#    ($self->{_xref_prot_file} = $arg );
#  return $self->{_xref_prot_file};
#}
#
#=head2 xref_dna_file
#
#  Arg [1]    : (optional) string $arg
#               the fasta file name for the dna xref
#  Example    : $file name = $xref->xref_dna_file();
#  Description: Getter / Setter for the dna xref fasta file 
#  Returntype : string
#  Exceptions : none
#
#=cut
#
#sub xref_dna_file{
#  my ($self, $arg) = @_;
#
#  (defined $arg) &&
#    ($self->{_xref_dna_file} = $arg );
#  return $self->{_xref_dna_file};
#}

=head2 ensembl_protein_file
 
  Arg [1]    : (optional) string $arg
               the fasta file name for the ensembl proteins 
  Example    : $file_name = $self->ensembl_protein_file();
  Description: Getter / Setter for the protien ensembl fasta file 
  Returntype : string
  Exceptions : none

=cut

sub ensembl_protein_file{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_ens_prot_file} = $arg );
  return $self->{_ens_prot_file};
}

=head2 ensembl_dna_file
 
  Arg [1]    : (optional) string $arg
               the fasta file name for the ensembl dna 
  Example    : $file_name = $self->ensembl_dna_file();
  Description: Getter / Setter for the protien ensembl fasta file 
  Returntype : string
  Exceptions : none

=cut

sub ensembl_dna_file{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_ens_dna_file} = $arg );
  return $self->{_ens_dna_file};
}

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


sub xref{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_xref} = $arg );
  return $self->{_xref};
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
  if (!defined($self->use_existing_mappings)) {
    my $dir = $self->dir();
    unlink (<$dir/*.map $dir/*.out $dir/*.err>);
  }

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
	my $job_name = $obj->run($queryfile, $targetfile, $self->dir());
	push @job_names, $job_name;
	sleep 1; # make sure unique names really are unique
      }
    }

  } # foreach method

  if (!defined($self->use_existing_mappings)) {
    # submit depend job to wait for all mapping jobs
    submit_depend_job($self->dir, @job_names);
  }

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

  # one -wended clause for each main job
  foreach my $job (@job_names) {
    push @depend_bsub, "-wended($job)";
  }

  # rest of command
  push @depend_bsub, ('-q', 'small', '-o', "$root_dir/depend.out", '-e', "$root_dir/depend.err", '/bin/true');

  #print "depend bsub:\n" . join (" ", @depend_bsub) . "\n";

  my ($depend_wtr, $depend_rtr, $depend_etr, $depend_pid);
  $depend_pid = open3($depend_wtr, $depend_rtr, $depend_etr, @depend_bsub);
  my $depend_jobid;
  while (<$depend_rtr>) {
    if (/Job <([0-9]+)> is/) {
      $depend_jobid = $1;
      print "LSF job ID for depend job: $depend_jobid \n" ;
    }
  }
  if (!defined($depend_jobid)) {
    print STDERR "Error: could not get depend job ID\n";
  }



}

=head2 parse_mappings

  Arg[1]     : The target file used in the exonerate run. Used to work out the Ensembl object type.
  Arg[2]     :
  Example    : none
  Description: Parse exonerate output files and build files for loading into target db tables.
  Returntype : List of strings
  Exceptions : none
  Caller     : general

=cut

sub parse_mappings {

  my ($self, $xref) = @_;

  my $dir = $self->dir();

  # get current max object_xref_id
  # TODO use selectall_arrayref
  my $row = @{$self->dbi()->selectall_arrayref("SELECT MAX(object_xref_id) FROM object_xref")}[0];
  my $max_object_xref_id = @{$row}[0];
  if (!defined $max_object_xref_id) {
    print "Can't get highest existing object_xref_id, using 1\n";
  } else {
    print "Maximum existing object_xref_id = $max_object_xref_id\n";
    $max_object_xref_id = 1;
  }

  $row = @{$self->dbi->selectall_arrayref("SELECT MAX(xref_id) FROM xref")}[0];
  my $max_xref_id = @$row[0];
  if (!defined $max_xref_id) {
    print "Can't get highest existing xref_id, using 0\n)";
  } else {
    print "Maximum existing xref_id = $max_xref_id\n";
    $max_object_xref_id = 1;
  }
  my $xref_id_offset = $max_xref_id + 1;

  #my $ox_sth = $dbi->prepare("INSERT INTO object_xref(ensembl_id, ensembl_object_type, xref_id) VALUES(?,?,?)");

  #my $ix_sth = $dbi->prepare("INSERT INTO identity_xref VALUES(?,?,?,?,?,?,?,?,?,?,?)");

  # files to write table data to
  open (OBJECT_XREF,   ">$dir/object_xref.txt");
  open (IDENTITY_XREF, ">$dir/identity_xref.txt");

  my $total_lines = 0;
  my $last_lines = 0;
  my $total_files = 0;

  my $object_xref_id = $max_object_xref_id + 1;

  # keep a (unique) list of xref IDs that need to be written out to file as well
  # this is a hash of hashes, keyed on xref id that relates xrefs to e! objects (may be 1-many)
  my %primary_xref_ids = ();

  # also keep track of types of ensembl objects
  my %ensembl_object_types;

  # and a list of mappings of ensembl objects to xrefs
  # (primary now, dependent added in dump_core_xrefs)
  # this is required for display_xref generation later
  # format:
  #   key: ensembl object type:ensembl object id
  #   value: list of xref_id (with offset)
  my %object_xref_mappings;


  foreach my $file (glob("$dir/*.map")) {

    #print "Parsing results from " . basename($file) .  "\n";
    open(FILE, $file);
    $total_files++;

    # files are named Method_(dna|peptide)_N.map
    my $type = get_ensembl_object_type($file);

    my $method = get_method($file);

    # get or create the appropriate analysis ID
    # XXX restore when using writeable database
    #my $analysis_id = $self->get_analysis_id($type);
    my $analysis_id = 999;

    while (<FILE>) {

      $total_lines++;
      chomp();
      my ($label, $query_id, $target_id, $identity, $query_length, $target_length, $query_start, $query_end, $target_start, $target_end, $cigar_line, $score) = split(/:/, $_);
      $cigar_line =~ s/ //g;

      # calculate percentage identities
      my $query_identity = int (100 * $identity / $query_length);
      my $target_identity = int (100 * $identity / $target_length);

      # TODO make sure query & target are the right way around

      # only take mappings where there is a good match on one or both sequences
      next if ($query_identity  < $method_query_threshold{$method} &&
	       $target_identity < $method_target_threshold{$method});

      # note we add on $xref_id_offset to avoid clashes
      print OBJECT_XREF "$object_xref_id\t$target_id\t$type\t" . ($query_id+$xref_id_offset) . "\n";
      print IDENTITY_XREF join("\t", ($object_xref_id, $query_identity, $target_identity, $query_start+1, $query_end, $target_start+1, $target_end, $cigar_line, $score, "\\N", $analysis_id)) . "\n";

      # TODO - evalue?
      $object_xref_id++;

      $ensembl_object_types{$target_id} = $type;

      # store mapping for later - note NON-OFFSET xref_id is used
      my $key = $type . "|" . $target_id;
      my $xref_id = $query_id;
      push @{$object_xref_mappings{$key}}, $xref_id;

      # note the NON-OFFSET xref_id is stored here as the values are used in
      # a query against the original xref database
      $primary_xref_ids{$query_id}{$target_id} = $target_id;
	
      # Store in database
      # create entry in object_xref and get its object_xref_id
      #$ox_sth->execute($target_id, $type, $query_id) || warn "Error writing to object_xref table";
      #my $object_xref_id = $ox_sth->{'mysql_insertid'};
	
      # create entry in identity_xref
      #$ix_sth->execute($object_xref_id, $query_id, $target_id, $query_start, $query_end, $target_start, $target_end, $cigar_line, $score, undef, $analysis_id) || warn "Error writing to identity_xref table";

    }

    close(FILE);
    #print "After $file, lines read increased by " . ($total_lines-$last_lines) . "\n";
    $last_lines = $total_lines;
  }

  close(IDENTITY_XREF);
  close(OBJECT_XREF);

  print "Read $total_lines lines from $total_files exonerate output files\n";

  # write relevant xrefs to file
  print "passing object_xref_mappings to dump_core_xrefs with " . scalar (keys %object_xref_mappings) . "\n";
  $self->dump_core_xrefs(\%primary_xref_ids, $object_xref_id+1, $xref_id_offset, \%ensembl_object_types, \%object_xref_mappings);

  # write comparison info. Can be removed after development
  $self->dump_comparison();

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

sub get_method {

  my $filename = shift;

  $filename = basename($filename);

  my ($method) = $filename =~ /^(.*)_(dna|peptide)_\d+\.map/;

  return $method;

}

sub get_analysis_id {

  my ($self, $ensembl_type) = @_;

  my %typeToLogicName = ( 'dna' => 'XrefExonerateDNA',
			  'protein' => 'XrefExonerateProtein' );

  my $logic_name = $typeToLogicName{lc($ensembl_type)};

  my $sth = $self->dbi()->prepare("SELECT analysis_id FROM analysis WHERE logic_name='" . $logic_name ."'");
  $sth->execute();

  my $analysis_id;

  if (my @row = $sth->fetchrow_array()) {

    $analysis_id = $row[0];
    print "Found exising analysis ID ($analysis_id) for $logic_name\n";

  } else {

    print "No analysis with logic_name $logic_name found, creating ...\n";
    $sth = $self->dbi()->prepare("INSERT INTO analysis (logic_name, created) VALUES ('" . $logic_name. "', NOW())");
    # TODO - other fields in analysis table
    $sth->execute();
    $analysis_id = $sth->{'mysql_insertid'};
    print "Done (analysis ID=" . $analysis_id. ")\n";

  }

  return $analysis_id;

}


sub dump_core_xrefs {

  my ($self, $xref_ids_hashref, $start_object_xref_id, $xref_id_offset, $ensembl_object_types_hashref, $object_xref_mappings) = @_;

  my @xref_ids = keys %$xref_ids_hashref;
  my %xref_to_objects = %$xref_ids_hashref;
  my %ensembl_object_types = %$ensembl_object_types_hashref;

  my $dir = $self->dir();

  open (XREF, ">$dir/xref.txt");
  open (OBJECT_XREF, ">>$dir/object_xref.txt");
  open (EXTERNAL_SYNONYM, ">$dir/external_synonym.txt");

  my $xref_dbi = $self->xref()->dbi();
  my $core_dbi = $self->dbi();

  # keep a unique list of source IDs to build the external_db table later
  my %source_ids;

  my $object_xref_id = $start_object_xref_id;

  # build cache of source id -> external_db id
  my %source_to_external_db = $self->map_source_to_external_db();

  # execute several queries with a max of 200 entries in each IN clause - more efficient
  my $batch_size = 200;

  # keep track of what xref_id & object_xref_ids have been written to prevent 
  # duplicates; e.g. several dependent xrefs may be dependent on the same master xref.
  my %xrefs_written;
  my %object_xrefs_written;

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
    my $xref_sth = $xref_dbi->prepare($sql);
    $xref_sth->execute();

    my ($xref_id, $accession, $version, $label, $description, $source_id, $species_id, $master_xref_id);
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
	  print XREF ($xref_id+$xref_id_offset) . "\t" . $external_db_id . "\t" . $accession . "\t" . $label . "\t" . $version . "\t" . $description . "\n";
	  $xrefs_written{$xref_id} = 1;
	  $source_ids{$source_id} = $source_id;
	}
      }
    }

    # Now get the dependent xrefs for each of these xrefs and write them as well
    $sql = "SELECT DISTINCT(x.xref_id), dx.master_xref_id, x.accession, x.label, x.description, x.source_id, x.version FROM dependent_xref dx, xref x WHERE x.xref_id=dx.dependent_xref_id AND master_xref_id $id_str";

    my $dep_sth = $xref_dbi->prepare($sql);
    $dep_sth->execute();

    $dep_sth->bind_columns(\$xref_id, \$master_xref_id, \$accession, \$label, \$description, \$source_id, \$version);
    while ($dep_sth->fetch()) {

      my $external_db_id = $source_to_external_db{$source_id};
      next if (!$external_db_id);

      $label = $accession if (!$label);

      if (!$xrefs_written{$xref_id}) {
	print XREF ($xref_id+$xref_id_offset) . "\t" . $external_db_id . "\t" . $accession . "\t" . $label . "\t" . $version . "\t" . $description . "\tDEPENDENT\n";
	$xrefs_written{$xref_id} = 1;
	$source_ids{$source_id} = $source_id;
      }

      # create an object_xref linking this (dependent) xref with any objects it maps to
      # write to file and add to object_xref_mappings
      if (defined $xref_to_objects{$master_xref_id}) {
	my @ensembl_object_ids = keys( %{$xref_to_objects{$master_xref_id}} ); 
	#print "xref $accession has " . scalar(@ensembl_object_ids) . " associated ensembl objects\n";
	foreach my $object_id (@ensembl_object_ids) {
	  my $type = $ensembl_object_types{$object_id};
	  my $full_key = $type."|".$object_id."|".$xref_id;
	  if (!$object_xrefs_written{$full_key}) {
	    print OBJECT_XREF "$object_xref_id\t$object_id\t$type\t" . ($xref_id+$xref_id_offset) . "\tDEPENDENT\n";
	    $object_xref_id++;
	    # Add this mapping to the list - note NON-OFFSET xref_id is used
	    my $key = $type . "|" . $object_id;
	    push @{$object_xref_mappings->{$key}}, $xref_id;
	    $object_xrefs_written{$full_key} = 1;
	  }
	}
      }
    }

    # Now get the synonyms for each of these xrefs and write them to the external_synonym table
    $sql = "SELECT DISTINCT xref_id, synonym FROM synonym WHERE xref_id $id_str";

    my $syn_sth = $xref_dbi->prepare($sql);
    $syn_sth->execute();

    $syn_sth->bind_columns(\$xref_id, \$accession);
    while ($syn_sth->fetch()) {

      print EXTERNAL_SYNONYM ($xref_id+$xref_id_offset) . "\t" . $accession . "\n";

    }

    #print "source_ids: " . join(" ", keys(%source_ids)) . "\n";

  } # while @xref_ids

  close(XREF);
  close(OBJECT_XREF);
  close(EXTERNAL_SYNONYM);

  print "Before calling display_xref, object_xref_mappings size " . scalar (keys %{$object_xref_mappings}) . "\n";

  # calculate display_xref_ids for transcripts and genes
  my $transcript_display_xrefs = $self->build_transcript_display_xrefs($object_xref_mappings, $xref_id_offset);

  $self->build_gene_display_xrefs($transcript_display_xrefs);

}


# produce output for comparison with existing ensembl mappings
# format is (with header)
# xref_accession ensembl_type ensembl_id

sub dump_comparison {

  my $self = shift;

  my $dir = $self->dir();

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

  my ($self, $object_xref_mappings, $xref_id_offset) = @_;

  my $dir = $self->dir();

  # get a list of xref sources; format:
  # key: xref_id value: source_name
  # lots of these; if memory is a problem, just get the source ID (not the name)
  # and look it up elsewhere
  print "Building xref->source mapping table\n";
  my %xref_to_source;
  my $sql = "SELECT x.xref_id, s.name FROM source s, xref x WHERE x.source_id=s.source_id";
  my $sth = $self->xref->dbi()->prepare($sql);
  $sth->execute();

  my ($xref_id, $source_name);
  $sth->bind_columns(\$xref_id, \$source_name);

  while ($sth->fetch()) {
    $xref_to_source{$xref_id} = $source_name;
  }

  print "Got " . scalar(keys %xref_to_source) . " xref-source mappings\n";

  # Cache the list of translation->transcript mappings & vice versa
  print "Building translation to transcript mappings\n";
  my %translation_to_transcript;
  my %transcript_to_translation;
  my $sth = $self->dbi()->prepare("SELECT translation_id, transcript_id FROM translation");
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

  foreach my $key (keys %{$object_xref_mappings}) {

    my ($type, $object_id) = split /\|/, $key;

    next if ($type !~ /(Transcript|Translation)/i);

    # if a transcript has more than one associated xref,
    # use the one with the highest priority, i.e. lower list position in @priorities
    my @xrefs = @{$object_xref_mappings->{$key}};
    my ($best_xref, $best_xref_priority_idx);
    $best_xref_priority_idx = 99999;
    foreach my $xref (@xrefs) {

      my $source = $xref_to_source{$xref};
      if ($source) {
	my $i = find_in_list($source, @priorities);
	if ($i > -1 && $i < $best_xref_priority_idx) {
	  $best_xref = $xref;
	  $best_xref_priority_idx = $i;
	}
      } else {
	warn("Couldn't find a source for xref $xref \n");
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
    if ($type =~ /Transcript/i) {
      my $transcript_id = $object_id;
      my $translation_id = $transcript_to_translation{$transcript_id};

      if ($translation_id) {
	my ($translation_xref, $translation_priority) = split /\|/, $obj_to_best_xref{"Translation|$translation_id"};
	my ($transcript_xref, $transcript_priority)   = split /\|/, $obj_to_best_xref{"Transcript|$transcript_id"};

	if ($translation_priority < $transcript_priority) {
	  $best_xref = $translation_xref;
	  $best_xref_priority_idx = $translation_priority;
	} else {
	  $best_xref = $transcript_xref;
	  $best_xref_priority_idx = $transcript_priority;
	}

      }
    }

    if ($best_xref) {

      # Write record with xref_id_offset
      print TRANSCRIPT_DX "UPDATE transcript SET display_xref_id=" . ($best_xref+$xref_id_offset) . " WHERE transcript_id=" . $object_id . ";\n";
      print "wrote " . $best_xref . " (plus offset) for 95625\n" if ($object_id eq 95625);
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

  my $dir = $self->dir();

  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-species => $self->species(),
					      -dbname  => $self->dbname(),
					      -host    => $self->host(),
					      -port    => $self->port(),
					      -pass    => $self->password(),
					      -user    => $self->user(),
					      -group   => 'core');
  my $ta = $db->get_TranscriptAdaptor();

  print "Building gene display_xrefs\n";
  print "Getting transcripts for all genes\n";

  my $sql = "SELECT gene_id, transcript_id FROM transcript";
  my $sth = $self->dbi()->prepare($sql);
  $sth->execute();

  my ($gene_id, $transcript_id);
  $sth->bind_columns(\$gene_id, \$transcript_id);

  my %genes_to_transcripts;
  while ($sth->fetch()) {
    push @{$genes_to_transcripts{$gene_id}}, $transcript_id;
  }

  print "Got " . scalar keys(%genes_to_transcripts) . " genes\n";

  print "Assigning display_xrefs to genes\n";

  open (GENE_DX, ">$dir/gene_display_xref.sql");
  open (GENE_DX_TXT, ">$dir/gene_display_xref.txt");
  my $hit = 0;
  my $miss = 0;
  my $trans_no_xref = 0;
  my $trans_xref = 0;
  foreach my $gene_id (keys %genes_to_transcripts) {

    my @transcripts = @{$genes_to_transcripts{$gene_id}};

    my $best_xref;
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
      #print "gene $gene_id orig:" . $transcript_display_xrefs->{$transcript_id} . " xref id: " . $xref_id . " pri " . $priority . "\n";
      # 2 separate if clauses to avoid having to fetch transcripts unnecessarily

      if (($priority lt $best_xref_priority_idx)) {

	$best_xref_priority_idx = $priority;
	$best_xref = $xref_id;

      } elsif ($priority eq $best_xref_priority_idx) {

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

    if ($best_xref) {
      # Write record
      print GENE_DX "UPDATE gene SET display_xref_id=" . $best_xref . " WHERE gene_id=" . $gene_id . ";\n";
      print GENE_DX_TXT $best_xref . "\t" . $gene_id ."\n";
      $hit++;
    }

  }

  close (GENE_DX);
  close (GENE_DX_TXT);
  print "Transcripts with no xrefs: $trans_no_xref with xrefs: $trans_xref\n";
  print "Wrote $hit gene display_xref entries to gene_display_xref.sql\n";
  print "Couldn't find display_xrefs for $miss genes\n" if ($miss > 0);
  print "Found display_xrefs for all genes\n" if ($miss eq 0);

}

# Display xref sources to be used for transcripts *in order of priority*
# Source names used must be identical to those in the source table.

sub transcript_display_xref_sources {

  return ('HUGO',
	  'MarkerSymbol',
	  'wormbase_transcript',
	  'flybase_symbol',
	  'Anopheles_symbol',
	  'Genoscope_predicted_transcript',
	  'Genoscope_predicted_gene',
	  'Uniprot/SWISSPROT',
	  'RefSeq',
	  'Uniprot/SPTREMBL',
	  'LocusLink');

}


# Find the index of an item in a list(ref), or 999999 if it's not in the list.
# Only look for exact matches (case insensitive)

sub find_in_list {

  my ($item, @list) = @_;

  for (my $i = 0; $i < scalar(@list); $i++) {
    if (lc($list[$i]) eq lc($item)) {
      return $i;
    }
  }

  return 999999;

}

# Build a map of source id (in xref database) to external_db (in core database)

sub map_source_to_external_db {

  my $self = shift;

  my %source_to_external_db;

  # get all sources
  my $sth = $self->xref->dbi()->prepare("SELECT source_id, name FROM source");
  $sth->execute();
  my ($source_id, $source_name);
  $sth->bind_columns(\$source_id, \$source_name);

  while($sth->fetchrow_array()) {

    # find appropriate external_db_id for each one
    my $sql = "SELECT external_db_id FROM external_db WHERE db_name=?";
    my $core_sth = $self->dbi()->prepare($sql);
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

# Upload .txt files and execute .sql files.

sub do_upload {

  my ($self, $deleteexisting) = @_;

  # xref.txt etc

  # TODO warn if table not empty

  foreach my $table ("xref", "object_xref", "identity_xref", "external_synonym") {

    my $file = $self->dir() . "/" . $table . ".txt";
    my $sth;

    if ($deleteexisting) {

      $sth = $self->dbi()->prepare("DELETE FROM $table");
      print "Deleting existing data in $table\n";
      $sth->execute();

    }

    # don't seem to be able to use prepared statements here
    $sth = $self->dbi()->prepare("LOAD DATA INFILE \'$file\' INTO TABLE $table");
    print "Uploading data in $file to $table\n";
    $sth->execute();

  }

  # gene_display_xref.sql etc
  foreach my $table ("gene", "transcript") {

    my $file = $self->dir() . "/" . $table . "_display_xref.txt";
    my $sth;

    if ($deleteexisting) {

      $sth = $self->dbi()->prepare("UPDATE $table SET display_xref_id=NULL");
      print "Setting all existing display_xref_id in $table to null\n";
      $sth->execute();

    }

    print "Setting $table display_xrefs from $file\n";
    # TODO this better
    #my $str = "mysql -u " .$self->user() ." -p" . $self->password() . " -h " . $self->host() ." -P " . $self->port() . " " .$self->dbname() . " < $file";
    #system $str;

    $sth = $self->dbi()->prepare("UPDATE $table SET display_xref_id=? WHERE ${table}_id=?");
    open(DX_TXT, $file);
    while (<DX_TXT>) {
      my ($xref_id, $object_id) = split;
      $sth->execute($xref_id, $object_id);
    }

    close(DX_TXT);
  }

}
1;
