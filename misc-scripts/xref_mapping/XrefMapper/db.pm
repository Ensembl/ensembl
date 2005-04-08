package XrefMapper::db;

use vars '@ISA';
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::ConfigRegistry;


@ISA = qw{Bio::EnsEMBL::DBSQL::DBAdaptor};


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

sub dumpcheck {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_dumpcheck} = $arg );
  return $self->{_dumpcheck};
}

sub maxdump {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_maxdump} = $arg );
  return $self->{_maxdump};
}

sub use_existing_mappings {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_use_existing_mappings} = $arg );
  return $self->{_use_existing_mappings};
}

sub process_dir {
  my ($dir) = @_;

  if($dir =~ "^\/" ) { # if it start with / then its not from pwd
    if(! -e $dir){
      die "directory does not exist $dir\n";
    }
  }
  elsif($dir eq "."){
    $dir = $ENV{PWD};
  }
  elsif($dir =~ "^\.\/"){
    my $tmp = $dir;
    $dir = $ENV{PWD}."/".substr($tmp,2);
    if(! -e $dir){
      die "directory does not exist $dir\n";
    }
  }
  else{
    die "directory does not exist $dir\n";
  }
  return $dir;
}


1;
