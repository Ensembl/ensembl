package XrefMapper::BasicMapper;

use strict;
use DBI;
use IPC::Open3;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Translation;
use  XrefMapper::db;

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
  my ($self, $xref) = @_;
  $self->dump_xref($xref);
  $self->dump_ensembl();
}



=head2 run_matching

  Arg[1]: xref object which holds info on method and files.

  Description: runs the mapping of the list of files with specied methods
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub run_matching{
  my ($self,$xref) = @_;

  my @list=();

  my $i = 0;
  foreach my $method (@{$self->method()}){
    my @dna=();
    push @dna, $method;
    push @dna, $xref->dir."/xref_".$i."_dna.fasta";
    push @dna, $self->ensembl_dna_file();
    push @list, \@dna;
    my @pep=();
    push @pep, $method;
    push @pep, $xref->dir."/xref_".$i."_prot.fasta";
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

#  return [["method1",["homo_sapiens","RefSeq"],["homo_sapiens","UniProtSwissProt"]],
#	  ["method2",[$self->species,"*"]],
#	  ["method3",["*","*"]]];

  return [["ExonerateBest1",["homo_sapiens","RefSeq"]]];

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
  my ($self,$xref) = @_;
  
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
  foreach my $list (@lists){
    print "method->".@$list[0]."\n";
    $method[$i] = shift @$list;
    my $j = 1;
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
	print $j."\t".$source. "\t".$source_id[$j] ."\n";
	print $j."\t".$species."\t".$species_id[$j]."\n";
	$j++;
      }
    }
    #method data fully defined now
    dump_subset($xref,\@species_id,\@source_id,$i);    
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
  my ($xref,$rspecies_id,$rsource_id,$index) = @_;
  
  open(XDNA,">".$xref->dir()."/xref_".$index."_dna.fasta") 
    || die "Could not open xref_".$index."_dna.fasta";

  my $sql = "select p.xref_id, p.sequence, x.species_id , x.source_id ";
  $sql   .= "  from primary_xref p, xref x ";
  $sql   .= "  where p.xref_id = x.xref_id and ";
  $sql   .= "      p.sequence_type ='dna' ";
  
  
  for (my $j =1; $j<scalar(@$rspecies_id); $j++){
    print $j."\t".$$rspecies_id[$j]."\t".$$rsource_id[$j]."\n";
  }
  #  return $xref->dir."/xref_".$i."_dna.fasta";
  
  my $sth = $xref->dbi()->prepare($sql);
  $sth->execute();
  my $i = 0;
  while(my @row = $sth->fetchrow_array()){
    my $pass = 0;
    for (my $j =1; $j<scalar(@$rspecies_id); $j++){
      if($$rspecies_id[$j] < 0 or $row[2] == $$rspecies_id[$j]){
	if($$rsource_id[$j] < 0 or  $row[3] == $$rsource_id[$j]){
	  $pass = 1;
	}
      }
    }
    if($pass){
      $i++;
      $row[1] =~ s/(.{60})/$1\n/g;
      print XDNA ">".$row[0]."\n".$row[1]."\n";
      if($i > 10){
	goto ENDDNA;
      }
    }
  }
 ENDDNA:
  close XDNA;


  open(XPRO,">".$xref->dir."/xref_".$index."_prot.fasta") 
    || die "Could not open xref_".$index."_prot.fasta";
  my $sql = "select p.xref_id, p.sequence, x.species_id , x.source_id ";
  $sql   .= "  from primary_xref p, xref x ";
  $sql   .= "  where p.xref_id = x.xref_id and ";
  $sql   .= "      p.sequence_type ='peptide' ";
  
  
  $sth = $xref->dbi()->prepare($sql);
  $sth->execute();
  $i = 0;
  while(my @row = $sth->fetchrow_array()){
    my $pass = 0;
    for (my $j =1; $j<scalar(@$rspecies_id); $j++){
      if($$rspecies_id[$j] < 0 or $row[2] == $$rspecies_id[$j]){
	if($$rsource_id[$j] < 0 or  $row[3] == $$rsource_id[$j]){
	  $pass = 1;
	}
      }
    }
    if($pass){
      $i++;
      $row[1] =~ s/(.{60})/$1\n/g;
      print XPRO ">".$row[0]."\n".$row[1]."\n";
      if($i > 10){
	goto ENDPRO;
      }
    }
  }
 ENDPRO:
  $sth->finish();
  close XPRO;

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
                           -password => $self->password(),
                           -user     => $self->user(),
                           -group    => 'core');

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
  open(PEP,">".$self->ensembl_protein_file()) 
    || die("Could not open dna file for writing: ".$self->ensembl_protein_file."\n");

  my $gene_adap = $db->get_GeneAdaptor();
  my @gene_ids = @{$gene_adap->list_dbIDs()};
  my $i =0;
  foreach my $gene_id (@gene_ids){
    $i++;
    my $gene = $gene_adap->fetch_by_dbID($gene_id);
#    print "gene ".$gene."\n";
    foreach my $transcript (@{$gene->get_all_Transcripts()}) {
      my $seq = $transcript->spliced_seq(); 
      $seq =~ s/(.{60})/$1\n/g;
      print DNA ">" . $transcript->dbID() . "\n" .$seq."\n";
      my $trans = $transcript->translation();
      my $translation = $transcript->translate();
#      print "tranlation ".$translation."\n";
      my $pep_seq = $translation->seq();
      $pep_seq =~ s/(.{60})/$1\n/g;
      print PEP ">".$trans->dbID()."\n".$pep_seq."\n";
    }
    if($i > 10){
      goto FIN;
    }
  }
FIN:
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

  # foreach method, submit the appropriate job & keep track of the job name
  my @job_names;

  foreach my $list (@$lists){

    my ($method, $queryfile ,$targetfile)  =  @$list;

    print "Method:$method Query:$queryfile Target:$targetfile\n";

    my $obj_name = "XrefMapper::Methods::$method";
    # check that the appropriate object exists
    eval "require $obj_name";
    if($@) {

      warn("Could not find object $obj_name corresponding to mapping method $method, skipping\n$@");

    } else {

      my $obj = $obj_name->new();
      my $job_name = $obj->run($queryfile, $targetfile, $self->dir());
      push @job_names, $job_name;
      print "Added LSF job $job_name to list\n";
      sleep 1; # make sure unique names really are unique

    }

  } # foreach method

  # submit depend job to wait for all mapping jobs
  submit_depend_job($self->dir, @job_names);


} # run_exonerate


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

=head2 store

  Arg[1]     : The target file used in the exonerate run. Used to work out the Ensembl object type.
  Arg[2]     :
  Example    : none
  Description: Parse exonerate output files and build files for loading into target db tables.
  Returntype : List of strings
  Exceptions : none
  Caller     : general

=cut

sub store {

  my ($self, $xref, $target_file_name) = @_;

  my $type = get_ensembl_object_type($target_file_name);

  # get or create the appropriate analysis ID
  my $analysis_id = $self->get_analysis_id($type);

  # get current max object_xref_id
  my $max_object_xref_id = 0;
  my $sth = $self->dbi()->prepare("SELECT MAX(object_xref_id) FROM object_xref");
  $sth->execute();
  my $max_object_xref_id = ($sth->fetchrow_array())[0];
  if (!defined $max_object_xref_id) {
    print "Can't get highest existing object_xref_id, using 1\n)";
  } else {
    print "Maximum existing object_xref_id = $max_object_xref_id\n";
  }

  #my $ox_sth = $dbi->prepare("INSERT INTO object_xref(ensembl_id, ensembl_object_type, xref_id) VALUES(?,?,?)");

  #my $ix_sth = $dbi->prepare("INSERT INTO identity_xref VALUES(?,?,?,?,?,?,?,?,?,?,?)");

  # files to write table data to
  open (OBJECT_XREF, ">object_xref.txt");
  open (IDENTITY_XREF, ">identity_xref.txt");

  my $total_lines = 0;
  my $total_files = 0;

  my $object_xref_id = $max_object_xref_id + 1;

  # keep a (unique) list of xref IDs that need to be written out to file as well
  my %primary_xref_ids;

  my $dir = $self->dir();
  foreach my $file (glob("$dir/*.map")) {

    print "Parsing results from $file \n";
    open(FILE, $file);
    $total_files++;

    while (<FILE>) {

      $total_lines++;
      chomp();
      my ($label, $query_id, $target_id, $query_start, $query_end, $target_start, $target_end, $cigar_line, $score) = split(/:/, $_);
      $cigar_line =~ s/ //;

      # TODO make sure query & target are the right way around

      print OBJECT_XREF "$object_xref_id\t$target_id\t$type\t$query_id\n";
      print IDENTITY_XREF "$object_xref_id\t$query_id\t$target_id\t$query_start\t$query_end\t$target_start\t$target_end\t$cigar_line\t$score\t\\N\t$analysis_id\n";
      # TODO - evalue?
      $object_xref_id++;

      $primary_xref_ids{$query_id} = $query_id;

      # Store in database
      # create entry in object_xref and get its object_xref_id
      #$ox_sth->execute($target_id, $type, $query_id) || warn "Error writing to object_xref table";
      #my $object_xref_id = $ox_sth->{'mysql_insertid'};

      # create entry in identity_xref
      #$ix_sth->execute($object_xref_id, $query_id, $target_id, $query_start, $query_end, $target_start, $target_end, $cigar_line, $score, undef, $analysis_id) || warn "Error writing to identity_xref table";

    }

    close(FILE);

  }

  close(IDENTITY_XREF);
  close(OBJECT_XREF);

  print "Read $total_lines lines from $total_files exonerate output files\n";

  # write relevant xrefs to file
  $self->dump_xrefs($xref, \%primary_xref_ids);

}


sub get_ensembl_object_type {

  my $filename = shift;
  my $type;

  if ($filename =~ /_dna\./i) {

    $type = "Transcript";

  } elsif ($filename =~ /_protein\./i) {

    $type = "Translation";

  } else {

    print STDERR "Cannot deduce Ensembl object type from filename $filename\n";
  }

print "###$filename   $type\n";
  return $type;

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


sub dump_xrefs {

  my ($self, $xref, $xref_ids_hashref) = @_;
  my @xref_ids = keys %$xref_ids_hashref;

  open (XREF, ">xref.txt");

  # TODO - get this from config
  my $xref_dbi = $xref->dbi();

  my $core_dbi = $self->dbi();

  # get current highest internal ID from xref
  my $max_xref_id = 0;
  my $core_sth = $core_dbi->prepare("SELECT MAX(xref_id) FROM xref");
  $core_sth->execute();
  my $max_xref_id = ($core_sth->fetchrow_array())[0];
  if (!defined $max_xref_id) {
    print "Can't get highest existing xref_id, using 0\n)";
  } else {
    print "Maximum existing xref_id = $max_xref_id\n";
  }
  my $core_xref_id = $max_xref_id + 1;

  # keep a unique list of source IDs to build the external_db table later
  my %source_ids;

  # execute several queries with a max of 200 entries in each IN clause - more efficient
  my $batch_size = 200;

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

    my ($xref_id, $accession, $label, $description, $source_id, $species_id);
    $xref_sth->bind_columns(\$xref_id, \$accession, \$label, \$description, \$source_id, \$species_id);

    # note the xref_id we write to the file is NOT the one we've just read
    # from the internal xref database as the ID may already exist in the core database
    while (my @row = $xref_sth->fetchrow_array()) {
      print XREF "$core_xref_id\t$accession\t$label\t$description\n";
      $source_ids{$source_id} = $source_id;
      $core_xref_id++;

    }

    # Now get the dependent xrefs for each of these xrefs and write them as well
    $sql = "SELECT x.accession, x.label, x.description, x.source_id FROM dependent_xref dx, xref x WHERE x.xref_id=dx.master_xref_id AND master_xref_id $id_str";
    my $dep_sth = $xref_dbi->prepare($sql);
    $dep_sth->execute();

    $dep_sth->bind_columns(\$accession, \$label, \$description, \$source_id);
    while (my @row = $dep_sth->fetchrow_array()) {
      print XREF "$core_xref_id\t$accession\t$label\t$description\tDEPENDENT\n";
      $source_ids{$source_id} = $source_id;
      $core_xref_id++;
    }
    #print "source_ids: " . join(" ", keys(%source_ids)) . "\n";

  } # while @xref_ids

  close(XREF);

  # now write the exernal_db file - the %source_ids hash will contain the IDs of the
  # sources that need to be written as external_dbs
  open(EXTERNAL_DB, ">external_db.txt");

  # get current highest internal ID from external_db
  my $max_edb_id = 0;
  my $core_sth = $core_dbi->prepare("SELECT MAX(external_db_id) FROM external_db");
  $core_sth->execute();
  my $max_edb_id = ($core_sth->fetchrow_array())[0];
  if (!defined $max_edb_id) {
    print "Can't get highest existing external_db_id, using 0\n)";
  } else {
    print "Maximum existing external_db_id = $max_edb_id\n";
  }
  my $edb_id = $max_edb_id + 1;

  my @source_id_array = keys %source_ids;
  my $source_id_str;
  if(@source_id_array > 1)  {
    $source_id_str = "IN (" . join(',', @source_id_array). ")";
  } else {
    $source_id_str = "= " . $source_id_array[0];
  }

  my $source_sql = "SELECT name, release FROM source WHERE source_id $source_id_str";
  my $source_sth = $xref_dbi->prepare($source_sql);
  $source_sth->execute();

  my ($source_name, $release);
  $source_sth->bind_columns(\$source_name, \$release);

  while (my @row = $source_sth->fetchrow_array()) {
    print EXTERNAL_DB "$edb_id\t$source_name\t$release\tXREF\n";
    # TODO knownxref etc??
    $edb_id++;
  }

  close(EXTERNAL_DB);



}

1;
