package XrefMapper::db;


sub new {
  my($class) = @_;

  my $self ={};
  bless $self,$class;

  print "creating new *$class*\n";
  return $self;
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

  my $dbi = DBI->connect("dbi:mysql:host=".$self->host().";port=".$self->port().";database=".$self->dbname(),
                        $self->user,
                        $self->password,
 			 {'RaiseError' => 1}) || die "Can't connect to database";


  return $dbi;
}

1;
