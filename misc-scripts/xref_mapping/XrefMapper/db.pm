package XrefMapper::db;


=head2 new

  Example    : $db = new XrefMapper::db(); 
  Description: Creates new db object.
  Returntype : XrefMapper::db;
  Exceptions : none
  Caller     : general
 
=cut


sub new {
  my($class) = @_;

  my $self ={};
  bless $self,$class;

  return $self;
}


=head2 species
 
  Arg [1]    : (optional) string $arg
               The new value of the species used by this connection.
  Example    : $host = $db->species()
  Description: Getter/Setter for the species of to use for
               this connection.  There is currently no point in setting
               this value after the connection has already been established
               by the constructor.
  Returntype : string
  Exceptions : none
  Caller     : new
 
=cut
 
sub species {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_species} = $arg );
  return $self->{_species};
}

=head2 host
 
  Arg [1]    : (optional) string $arg
               The new value of the host used by this connection.
  Example    : $host = $db->host()
  Description: Getter/Setter for the domain name of the database host use by
               this connection.  There is currently no point in setting
               this value after the connection has already been established
               by the constructor.
  Returntype : string
  Exceptions : none
  Caller     : new
 
=cut
 
sub host {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_host} = $arg );
  return $self->{_host};
}

=head2 port
 
  Arg [1]    : (optional) int $arg
               the TCP or UDP port to use to connect to the database
  Example    : $port = $db->port();
  Description: Getter / Setter for the port this connection uses to communicate
               to the database daemon.  There currently is no point in
               setting this value after the connection has already been
               established by the constructor.
  Returntype : string
  Exceptions : none
  Caller     : new
 
=cut
 

sub port {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_port} = $arg );
  return $self->{_port};
}

=head2 dbname
 
  Arg [1]    : (optional) string $arg
               The new value of the database name used by this connection.
  Example    : $dbname = $db->dbname()
  Description: Getter/Setter for the name of the database used by this
               connection.  There is currently no point in setting this value
               after the connection has already been established by the
               constructor.
  Returntype : string
  Exceptions : none
  Caller     : new
 
=cut
 
sub dbname {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_dbname} = $arg );
  return $self->{_dbname};
}

=head2 user
                                                                                
  Arg [1]    : (optional) string $arg
               The new value of the username used by this connection.
  Example    : $user = $db->user()
  Description: Getter/Setter for the username used by this
               connection.  There is currently no point in setting this value
               after the connection has already been established by the
               constructor.
  Returntype : string
  Exceptions : none
  Caller     : new
                                                                                
=cut
                                                                                
sub user {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_user} = $arg );
  return $self->{_user};
}

=head2 password
 
  Arg [1]    : (optional) string $arg
               The new value of the password used by this connection.
  Example    : $pass = $db->password()
  Description: Getter/Setter for the password of to use for
               this connection.  There is currently no point in setting
               this value after the connection has already been established
               by the constructor.
  Returntype : string
  Exceptions : none
  Caller     : new
 
=cut
 

sub password {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_password} = $arg );
  return $self->{_password};
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
  print STDERR $dir."\n";
  return $dir;
}

=head2 dbi
                                                                                
  Example    : $dbi = $db->dbi()
  Description: creates a new dbi connection
  Returntype : DBI
  Exceptions : will die if database is not connected to properly.
  Caller     : general
                                                                                
=cut
                                                                                
sub dbi {
  
  my $self = shift;
  
  my $dbi = DBI->connect("dbi:mysql:host=".$self->host().";port=".$self->port().";database=".$self->dbname(),
			 $self->user,
                        $self->password,
 			 {'RaiseError' => 1}) || die "Can't connect to database";


  return $dbi;
}

1;
