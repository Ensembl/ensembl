package BasicMapper;
 
use strict;
use DBI;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
 
 
sub new {
  my($class, $species, $host, $port, $dbname, $user, $password ,$dir) = @_;
 
  my $self ={};
  bless $self,$class;
 
  $self->species($species);
  $self->host($host);
  $self->port($port);
  $self->dbname($dbname);
  $self->$user($user);
  $self->$password($password);
  $self->$dir($dir);
}


sub dump_seqs{
  my ($self, $slice) = @_;
  $self->dump_ensembl($slice);
  $self->dump_xref();
}


sub dump_xref{
  my ($self) = @_;


  my $sql = "select species_id from species where name = '".$self->species."'";
  my $sth = dbi()->prepare($sql);
  $sth->execute();
  my @row = $sth->fetchrow_array();
  my $species_id;
  if (defined @row) {
    $species_id = $row[0];
  } else {
    print STDERR "Couldn't get ID for species ".$self->species."\n";
    print STDERR "It must be one of :-\n";
    $sql = "select name from species";
    $sth = dbi()->prepare($sql);
    $sth->execute();
    while(my @row = $sth->fetchrow_array()){
      print STDERR $row[0]."\n";
    }
    die("Please try again\n");
  }
                                                                                                                            

  open(XDNA,">".$self->dir."/xref_dna.fasta") || die "Could not open xref_dna.fasta";
  my $sql = "select p.xref_id, p.sequence from primary_xref p, xref x ";
  $sql   .= "where p.xref_id = x.xref_id and ";
  $sql   .= "      p.sequence_type ='dna' and ";
  $sql   .= "      x.species_id = ".$species_id." ";
  $sth = dbi()->prepare($sql);
  $sth->execute();
  while(my @row = $sth->fetchrow_array()){
    $row[1] =~ s/(.{60})/$1\n/g;
    print XDNA ">".$row[0]."\n".$row[1]."\n";
  }
  close XDNA;

  open(XPEP,">".$self->dir."/xref_prot.fasta") || die "Could not open xref_prot.fasta";
  $sql    = "select p.xref_id, p.sequence from primary_xref p, xref x ";
  $sql   .= "where p.xref_id = x.xref_id and ";
  $sql   .= "      p.sequence_type ='peptide' and ";
  $sql   .= "      x.species_id = ".$species_id." ";
  $sth = dbi()->prepare($sql);
  $sth->execute();
  while(my @row = $sth->fetchrow_array()){
    $row[1] =~ s/(.{60})/$1\n/g;
    print XPEP ">".$row[0]."\n".$row[1]."\n";
  }
  close XPEP;

}


sub dump_ensembl{
  my ($self) = @_;

  #create filename
  $self->ensembl_protein_file($self->dir."/".$self->species."_protein.fasta");
  $self->ensembl_dna_file($self->dir."/".$self->species."_dna.fasta");


  $self->fetch_and_dump_seq();
  
}


sub fetch_and_dump_seq{
  my ($self, $type, $adaptortype) = @_;
 
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-species => $self->species(),
                           -dbname  => $self->dbname(),
                           -host    => $self->host(),
                           -port    => $self->port(),
                           -password => $self->password(),
                           -username => $self->username(),
                           -group    => 'core');
 
  open(FILE,">".$self->ensembl_dna_file()) 
    || die("Could not open dna file for writing: ".$self->ensembl_dna_file."\n");

  $gene_adap = $reg->get_adaptor($self->species(),'core','Gene');
  my @genes = @{$gene_adap->list_dbIDs()};
  foreach my $gene (@genes){
    foreach my $transcript (@{$gene->get_all_Transcripts()}) {
      my $seq = $transcript->spliced_seq(); 
      $seq =~ s/(.{60})/$1\n/g;
      print FILE ">" . $transcript->dbID() . "\n" .$seq."\n";
    }
  }
  close FILE;

  #now do the translations.

}
                                      


# @transcript_ids = @{$transcript_adaptor->list_dbIDs()};
#  $transcript = $transcript_adaptor->fetch_by_dbID($trans_id);



#  $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(...);
#  $slice_adaptor = $db->get_SliceAdaptor();
#                                                                                
#  $transcript_adaptor = $db->get_TranscriptAdaptor();
#                                                                                
#  $transcript = $transcript_adaptor->fetch_by_dbID(1234);
#                                                                                
#  $transcript = $transcript_adaptor->fetch_by_stable_id('ENST00000201961');
#                                                                                
#  $slice = $slice_adaptor->fetch_by_region('chromosome', '3', 1, 1000000);
#  @transcripts = @{$transcript_adaptor->fetch_all_by_Slice($slice)};
#                                                                                
#  ($transcript) = @{$transcript_adaptor->fetch_all_by_external_name('BRCA2')};
                                                                                

sub get_ensembl_type{

  my %type; # 

  my $type{'Translation'} = "peptide";
  $type{'Transcript'} = "dna";

  return /%type;
}

sub species {
  my ($self, $arg) = @_;
 
  (defined $arg) &&
    ($self->{_species} = $arg );
  return $self->{_species};
}
 
sub host {
  my ($self, $arg) = @_;
 
  (defined $arg) &&
    ($self->{_host} = $arg );
  return $self->{_host};
}
 

sub port {
  my ($self, $arg) = @_;
 
  (defined $arg) &&
    ($self->{_port} = $arg );
  return $self->{_port};
}

sub dbname {
  my ($self, $arg) = @_;
 
  (defined $arg) &&
    ($self->{_dbname} = $arg );
  return $self->{_dbname};
}
 
sub user {
  my ($self, $arg) = @_;
 
  (defined $arg) &&
    ($self->{_user} = $arg );
  return $self->{_user};
}
 
sub password {
  my ($self, $arg) = @_;
 
  (defined $arg) &&
    ($self->{_password} = $arg );
  return $self->{_password};
}

sub dir {
  my ($self, $arg) = @_;
 
  (defined $arg) &&
    ($self->{_dir} = $arg );
  return $self->{_dir};
}
 
