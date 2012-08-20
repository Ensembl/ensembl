package XrefMapper::db;

use vars '@ISA';
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use XrefParser::Database qw//;
use Cwd;

sub new{
  my($class, @args) = @_;

  my $self ={};
  bless $self,$class;

   $self->dbc(new Bio::EnsEMBL::DBSQL::DBConnection(@args));

  return $self;
} 

sub species{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_species} = $arg );
  return $self->{_species};

}

sub dbc{
  my $self  = shift;

  if(@_){
    my $arg = shift;
    if(defined($arg)){
      if(!$arg->isa('Bio::EnsEMBL::DBSQL::DBConnection')){
        throw("$arg is no a DBConnection\n");
      }
    }
    $self->{_dbc} = $arg;
  }
  return $self->{_dbc};
}

=head2 dba

Produce a DBAdaptor from this db object. Useful for when you need to work
with a DB. Group is optional but defaults to core.

=cut

sub dba {
  my ($self, $group) = @_;
  my $dbc = $self->dbc();
  my $species = $self->species();
  $group ||= 'core';
  my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -DBCONN => $dbc, -SPECIES => $species, -GROUP => $group
  );
}

=head2 dir
                                                                                
  Arg [1]    : (optional) string $arg
               The new value of the dir used 
  Example    : $dir = $db->dir()
  Description: Getter/Setter for the directory used in the creation of fasta file
  Returntype : string
  Exceptions : none
  Caller     : new
                                                                                
=cut
                                                                                
sub dir {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_dir} = process_dir($arg) );
  return $self->{_dir};

}

=head2 protein_file
 
  Arg [1]    : (optional) string $arg
               the fasta file name for the ensembl proteins 
  Example    : $file_name = $self->ensembl_protein_file();
  Description: Getter / Setter for the protien ensembl fasta file 
  Returntype : string
  Exceptions : none

=cut

sub protein_file{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_ens_prot_file} = $arg );
  return $self->{_ens_prot_file};
}

=head2 dna_file
 
  Arg [1]    : (optional) string $arg
               the fasta file name for the ensembl dna 
  Example    : $file_name = $self->ensembl_dna_file();
  Description: Getter / Setter for the protien ensembl fasta file 
  Returntype : string
  Exceptions : none

=cut

sub dna_file{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_ens_dna_file} = $arg );
  return $self->{_ens_dna_file};
}

sub process_dir {
  my ($dir) = @_;

  if($dir =~ "^\/" ) { # if it start with / then its not from pwd
    if(! -d $dir){
      die "directory does not exist $dir\n";
    }
  }
  elsif($dir eq "."){
    $dir = cwd();
  }
  elsif($dir =~ "^\.\/"){
    my $tmp = $dir;
    $dir = cwd() . "/" . substr( $tmp, 2 );
    if(! -d $dir){
      die "directory does not exist $dir\n";
    }
  }
  else{
    die "directory does not exist $dir\n";
  }
  return $dir;
}

########## START OF XrefParser::Database compatiblity mode ##########
sub host {
  my ($self) = @_;
  return $self->dbc()->host();
}

sub port {
  my ($self) = @_;
  return $self->dbc()->port();
}

sub user {
  my ($self) = @_;
  return $self->dbc()->username();
}

sub pass {
  my ($self) = @_;
  return $self->dbc()->password();
}

sub dbname {
  my ($self) = @_;
  return $self->dbc()->dbname();
}

sub verbose {
  my ($self, $verbose) = @_;
  $self->{'verbose'} = $verbose if defined $verbose;
  return $self->{'verbose'};
}

sub dbi {
  my ($self) = @_;
  return XrefParser::Database::dbi($self);
}

sub create {
  my ($self, $sql_dir, $force, $drop_db) = @_;
  return XrefParser::Database::create($self, $sql_dir, $force, $drop_db);
}

########## END OF XrefParser::Database compatiblity mode ##########


1;
