package XrefMapper::BasicMapper;
 
use strict;
use DBI;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Translation;

my $xref_host = "ecs1g";
my $xref_port = 3306;
my $xref_database = "ianl_test_xref";
my $xref_user = "ensadmin";
my $xref_password = "ensembl";
 
 
sub new {
  my($class, $species, $host, $port, $dbname, $user, $password ,$dir) = @_;
 
  my $self ={};
  bless $self,$class;
 
  $self->species($species);
  $self->host($host);
  $self->port($port);
  $self->dbname($dbname);
  $self->user($user);
  $self->password($password);
  $self->dir($dir);


  return $self;
}




sub dump_seqs{
  my ($self, $slice) = @_;
  $self->dump_ensembl($slice);
  $self->dump_xref();
}


sub dump_xref{
  my ($self) = @_;


  my $sql = "select species_id from species where name = '".$self->species."'";
  my $dbi = $self->dbi();
  my $sth = $dbi->prepare($sql);
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
  my $i = 0;
  while(my @row = $sth->fetchrow_array()){
    $i++;
    $row[1] =~ s/(.{60})/$1\n/g;
    print XDNA ">".$row[0]."\n".$row[1]."\n";
    if($i > 10){
      goto ENDDNA;
    }
  }
ENDDNA:
  close XDNA;

  open(XPEP,">".$self->dir."/xref_prot.fasta") || die "Could not open xref_prot.fasta";
  $sql    = "select p.xref_id, p.sequence from primary_xref p, xref x ";
  $sql   .= "where p.xref_id = x.xref_id and ";
  $sql   .= "      p.sequence_type ='peptide' and ";
  $sql   .= "      x.species_id = ".$species_id." ";
  $sth = dbi()->prepare($sql);
  $sth->execute();
  $i = 0;
  while(my @row = $sth->fetchrow_array()){
    $i++;
    $row[1] =~ s/(.{60})/$1\n/g;
    print XPEP ">".$row[0]."\n".$row[1]."\n";
    if($i > 10){
      goto ENDXPEP;
    }
  }
ENDXPEP:
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
                           -user     => $self->user(),
                           -group    => 'core');
 
  open(DNA,">".$self->ensembl_dna_file()) 
    || die("Could not open dna file for writing: ".$self->ensembl_dna_file."\n");

  open(PEP,">".$self->ensembl_protein_file()) 
    || die("Could not open dna file for writing: ".$self->ensembl_protein_file."\n");

  my $gene_adap = $db->get_GeneAdaptor();
  my @gene_ids = @{$gene_adap->list_dbIDs()};
  my $i =0;
  foreach my $gene_id (@gene_ids){
    $i++;
    my $gene = $gene_adap->fetch_by_dbID($gene_id);
    print "gene ".$gene."\n";
    foreach my $transcript (@{$gene->get_all_Transcripts()}) {
      my $seq = $transcript->spliced_seq(); 
      $seq =~ s/(.{60})/$1\n/g;
      print DNA ">" . $transcript->dbID() . "\n" .$seq."\n";
      my $trans = $transcript->translation();
      my $translation = $transcript->translate();
      print "tranlation ".$translation."\n";
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
                                      


sub ensembl_protein_file{
  my ($self, $arg) = @_;
 
  (defined $arg) &&
    ($self->{_ens_prot_file} = $arg );
  return $self->{_ens_prot_file};
}

sub ensembl_dna_file{
  my ($self, $arg) = @_;
 
  (defined $arg) &&
    ($self->{_ens_dna_file} = $arg );
  return $self->{_ens_dna_file};
}


sub get_ensembl_type{
  my %type;

  $type{'Translation'} = "peptide";
  $type{'Transcript'} = "dna";

  return \%type;
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


sub dbi {
 
  my $self = shift;
  
  my $dbi = DBI->connect("dbi:mysql:host=$xref_host;port=$xref_port;database=$xref_database",
                        "$xref_user",
                        "$xref_password",
 			 {'RaiseError' => 1}) || die "Can't connect to database";
  
 
  return $dbi;
}
 
 
1;
