package XrefMapper::BasicMapper;

use strict;
use DBI;
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
  
  glenn_call(\@list);
  print "NOT done yet:-)\n";
}


sub glenn_call{
  my ($lists) = @_;
  foreach my $list (@$lists){
    my ($meth, $q ,$t)  =  @$list;
    print "HELLO: ".$meth."\n\t".$q."\n\t".$t."\n\n";
  }
}

sub store{
  print "NOT done yet Either :-)\n";
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

  return [["method1",["homo_sapiens","RefSeq"],["homo_sapiens","UniProtSwissProt"]],
	  ["method2",[$self->species,"*"]],
	  ["method3",["*","*"]]];
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

1;
