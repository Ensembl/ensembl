=head1 NAME - Bio::EnsEMBL::DBSQL::DBConnection

=head1 SYNOPSIS

    $db = Bio::EnsEMBL::DBSQL::DBConnection->new(
        -user   => 'anonymous',
        -dbname => 'homo_sapiens_core_20_34c',
        -host   => 'ensembldb.ensembl.org',
        -driver => 'mysql',
        );


   You should use this as a base class for all objects (DBAdaptor) that 
   connect to a database.  SQL statements should be created/executed through
   this modules prepare() and do() methods.

   $sth = $db->prepare( "SELECT something FROM yourtable" );

=head1 DESCRIPTION

  This class is a wrapper around DBIs datbase handle.  It provides some
  additional functionality such as the ability to automatically disconnect
  when inactive and reconnect when needed, and a way to obtain ObjectAdaptors.
  This class is intended to be inherited from rather than used directly.
  See Bio::EnsEMBL::DBSQL::DBAdaptor.

=head1 CONTACT

  This module is part of the Ensembl project: www.ensembl.org

  Ensembl development mailing list: <ensembl-dev@ebi.ac.uk>

=head1 METHODS

=cut


package Bio::EnsEMBL::DBSQL::DBConnection;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Container;
use Bio::EnsEMBL::Root;
use DBI;

use Bio::EnsEMBL::DBSQL::StatementHandle;

use Bio::EnsEMBL::Utils::Exception qw(throw info warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw(Bio::EnsEMBL::Root); # for backwards compatibility


=head2 new

  Arg [DBNAME] : string
                 The name of the database to connect to.
  Arg [HOST] : (optional) string
               The domain name of the database host to connect to.  
               'localhost' by default. 
  Arg [USER] : string
               The name of the database user to connect with 
  Arg [PASS] : (optional) string
               The password to be used to connect to the database
  Arg [PORT] : int
               The port to use when connecting to the database
               3306 by default.
  Arg [DRIVER] : (optional) string
                 The type of database driver to use to connect to the DB
                 mysql by default.
  Arg [DBCONN] : (optional)
                 Open another handle to the same database as another connection
                 If this argument is specified, no other arguments should be
                 specified.
  Arg [DISCONNECT_WHEN_INACTIVE]: (optional) boolean
                 If set to true, the database connection will be disconnected
                 everytime there are no active statement handles. This is
                 useful when running a lot of jobs on a compute farm
                 which would otherwise keep open a lot of connections to the
                 database.  Database connections are automatically reopened
                 when required.
  Example    : $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new
                  (-user   => 'anonymous',
                   -dbname => 'homo_sapiens_core_20_34c',
                   -host   => 'ensembldb.ensembl.org',
                   -driver => 'mysql');

  Description: Constructor for a DatabaseConenction. Any adaptors that require
               database connectivity should inherit from this class.
  Returntype : Bio::EnsEMBL::DBSQL::DBConnection
  Exceptions : thrown if USER or DBNAME are not specified, or if the database
               cannot be connected to.
  Caller     : Bio::EnsEMBL::DBSQL::DBAdaptor

=cut

sub new {
  my $class = shift;

  my ($db,$host,$driver,$user,$password,$port, $inactive_disconnect, $dbconn) =
    rearrange([qw(DBNAME HOST DRIVER USER PASS PORT
                  DISCONNECT_WHEN_INACTIVE DBCONN)], @_);

  my $self = {};
  bless $self, $class;

  if($dbconn) {
    if($db || $host || $driver || $password || $port || $inactive_disconnect) {
      throw("Cannot specify other arguments when -DBCONN argument used.");
    }

    $self->dbname($dbconn->dbname());
    $self->username($dbconn->username());
    $self->host($dbconn->host());
    $self->password($dbconn->password());
    $self->port($dbconn->port());
    $self->driver($dbconn->driver());

    $self->connect();

    if($dbconn->disconnect_when_inactive()) {
      $self->disconnect_when_inactive(1) 
    }
  } else {

    $db   || throw("-DBNAME argument is required.");
    $user || throw("-USER argument is required.");

    $driver ||= 'mysql';
    $host   ||= 'mysql';
    $port   ||= 3306;

    my $dsn = "DBI:$driver:database=$db;host=$host;port=$port";

    $self->username( $user );
    $self->host( $host );
    $self->dbname( $db );
    $self->password( $password );
    $self->port($port);
    $self->driver($driver);

    $self->connect();

    if($inactive_disconnect) {
      $self->disconnect_when_inactive($inactive_disconnect);
    }

    # be very sneaky and actually return a container object which is outside
    # of the circular reference loops and will perform cleanup when all 
    # references to the container are gone.
  }

  return new Bio::EnsEMBL::Container($self);
}


=head2 connect

  Example    : $dbcon->connect()
  Description: Connects to the database using the connection attribute 
               information.
  Returntype : none
  Exceptions : none
  Caller     : new, db_handle

=cut

sub connect {
  my $self = shift;

  my $dsn = "DBI:" . $self->driver() .
            ":database=". $self->dbname() .
            ";host=" . $self->host() .
            ";port=" . $self->port();

  my $dbh;
  eval{
    $dbh = DBI->connect($dsn,
                        $self->username(),
                        $self->password(),
                        {RaiseError => 1});
  };

  $dbh || throw("Could not connect to database " . $self->dbname() .
                " as user " . $self->username() .
                " using [$dsn] as a locator:\n" . $DBI::errstr);

  $self->db_handle($dbh);

  return;
}


=head2 driver

  Arg [1]    : (optional) string $arg
               the name of the driver to use to connect to the database
  Example    : $driver = $db_connection->driver()
  Description: Getter / Setter for the driver this connection uses.
               Right now there is no point to setting this value after a
               connection has already been established in the constructor.
  Returntype : string
  Exceptions : none
  Caller     : new

=cut

sub driver {
  my($self, $arg ) = @_;

  (defined $arg) &&
    ($self->{_driver} = $arg );
  return $self->{_driver};
}


=head2 port

  Arg [1]    : (optional) int $arg
               the TCP or UDP port to use to connect to the database
  Example    : $port = $db_connection->port();
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
  Example    : $dbname = $db_connection->dbname()
  Description: Getter/Setter for the name of the database used by this 
               connection.  There is currently no point in setting this value
               after the connection has already been established by the 
               constructor.
  Returntype : string
  Exceptions : none
  Caller     : new

=cut

sub dbname {
  my ($self, $arg ) = @_;
  ( defined $arg ) &&
    ( $self->{_dbname} = $arg );
  $self->{_dbname};
}


=head2 username

  Arg [1]    : (optional) string $arg
               The new value of the username used by this connection. 
  Example    : $username = $db_connection->username()
  Description: Getter/Setter for the username used by this 
               connection.  There is currently no point in setting this value
               after the connection has already been established by the 
               constructor.
  Returntype : string
  Exceptions : none
  Caller     : new

=cut

sub username {
  my ($self, $arg ) = @_;
  ( defined $arg ) &&
    ( $self->{_username} = $arg );
  $self->{_username};
}


=head2 host

  Arg [1]    : (optional) string $arg
               The new value of the host used by this connection. 
  Example    : $host = $db_connection->host()
  Description: Getter/Setter for the domain name of the database host use by 
               this connection.  There is currently no point in setting 
               this value after the connection has already been established 
               by the constructor.
  Returntype : string
  Exceptions : none
  Caller     : new

=cut

sub host {
  my ($self, $arg ) = @_;
  ( defined $arg ) &&
    ( $self->{_host} = $arg );
  $self->{_host};
}


=head2 password

  Arg [1]    : (optional) string $arg
               The new value of the password used by this connection. 
  Example    : $host = $db_connection->password()
  Description: Getter/Setter for the password of to use for 
               this connection.  There is currently no point in setting 
               this value after the connection has already been established 
               by the constructor.
  Returntype : string
  Exceptions : none
  Caller     : new

=cut

sub password {
  my ($self, $arg ) = @_;
  ( defined $arg ) &&
    ( $self->{_password} = $arg );
  $self->{_password};
}


=head2 disconnect_when_inactive

  Arg [1]    : (optional) boolean $newval
  Example    : $db->disconnect_when_inactive(1);
  Description: Getter/Setter for the disconnect_when_inactive flag.  If set
               to true this DBConnection will continually disconnect itself
               when there are no active statement handles and reconnect as
               necessary.  Useful for farm environments when there can be
               many (often inactive) open connections to a database at once.
  Returntype : boolean
  Exceptions : none
  Caller     : Pipeline

=cut

sub disconnect_when_inactive {
  my $self = shift;

  if(@_) {
    my $val = shift;
    $self->{'disconnect_when_inactive'} = $val;
    if($val) {
      $self->disconnect_if_idle();
    } elsif(!$self->db_handle->ping()) {
      # reconnect if the connect went away
      $self->connect();
    }
  }

  return $self->{'disconnect_when_inactive'};
}



=head2 locator

  Arg [1]    : none
  Example    : $locator = $dbc->locator;
  Description: Constructs a locator string for this database connection
               that can, for example, be used by the DBLoader module
  Returntype : string
  Exceptions : none
  Caller     : general

=cut


sub locator {
  my $self = shift;

  my $ref = ref($self);

  return "$ref/host=".$self->host.";port=".$self->port.";dbname=".
    $self->dbname.";user=".$self->username.";pass=".$self->password;
}


=head2 db_handle

  Arg [1]    : DBI Database Handle $value
  Example    : $dbh = $db_connection->db_handle() 
  Description: Getter / Setter for the Database handle used by this
               database connection.
  Returntype : DBI Database Handle
  Exceptions : none
  Caller     : new, DESTROY

=cut

sub db_handle {
   my $self = shift;

   $self->{'db_handle'} = shift if(@_);

   return $self->{'db_handle'};
}



=head2 prepare

  Arg [1]    : string $string
               the SQL statement to prepare
  Example    : $sth = $db_connection->prepare("SELECT column FROM table");
  Description: Prepares a SQL statement using the internal DBI database handle
               and returns the DBI statement handle.
  Returntype : DBI statement handle
  Exceptions : thrown if the SQL statement is empty, or if the internal
               database handle is not present
  Caller     : Adaptor modules

=cut

sub prepare {
   my ($self,$string) = @_;

   if( ! $string ) {
     throw("Attempting to prepare an empty SQL query.");
   }
   if( ! defined($self->db_handle) ) {
     throw("Database object has lost its database handle.");
   }

   if($self->disconnect_when_inactive() && !$self->db_handle->ping()) {
     # reconnect if we have been disconnected
     $self->connect();
   }

   #info("SQL(".$self->dbname."):$string");

   my $sth = $self->db_handle->prepare($string);

   # return an overridden statement handle that provides us with
   # the means to disconnect inactive statement handles automatically
   bless $sth, "Bio::EnsEMBL::DBSQL::StatementHandle";

   $sth->dbc($self);

   return $sth;
}


=head2 do

  Arg [1]    : string $string
               the SQL statement to prepare
  Example    : $sth = $db_connection->do("SELECT column FROM table");
  Description: Executes a SQL statement using the internal DBI database handle.
  Returntype : Result of DBI dbh do() method
  Exceptions : thrown if the SQL statement is empty, or if the internal
               database handle is not present.
  Caller     : Adaptor modules

=cut

sub do {
   my ($self,$string) = @_;

   if( ! $string ) {
     throw("Attempting to do an empty SQL query.");
   }
   if( ! defined($self->db_handle) ) {
     throw("Database object has lost its database handle.");
   }

   if($self->disconnect_when_inactive() && !$self->db_handle->ping()) {
     # reconnect if we have been disconnected
     $self->connect();
   }

   #info("SQL(".$self->dbname."):$string");

   my $result = $self->db_handle->do($string);

   # disconnect if the disconnect when inactive flag is set and
   # there are no active statement handles

   if($self->disconnect_when_inactive()) {
     $self->disconnect_if_idle();
   }

   return $result;
}



=head2 disconnect_if_idle

  Arg [1]    : none
  Example    : $dbc->disconnect_if_idle();
  Description: Disconnects from the database if there are no currently active
               statement handles.  You should not call this method directly.
               It is called automatically by the DESTROY method of the
               Bio::EnsEMBL::DBSQL::SQL::StatementHandle if the
               disconect_when_inactive flag is set.
  Returntype : none
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::SQL::StatementHandle::DESTROY
               Bio::EnsEMBL::DBSQL::DBConnection::do

=cut

sub disconnect_if_idle {
  my $self = shift;

  if(!$self->db_handle()->{'ActiveKids'} &&
     !$self->db_handle()->{'InactiveDestroy'}) {
    $self->db_handle->disconnect();
  }
}


=head2 _get_adaptor

  Arg [1]    : string $module
               the fully qualified of the adaptor module to be retrieved
  Arg [2..n] : (optional) arbitrary list @args
               list of arguments to be passed to adaptors constructor
  Example    : $adaptor = $self->_get_adaptor("full::adaptor::name");
  Description: PROTECTED Used by subclasses to obtain adaptor objects
               for this database connection using the fully qualified
               module name of the adaptor. If the adaptor has not been 
               retrieved before it is created, otherwise it is retreived
               from the adaptor cache.
  Returntype : Adaptor Object of arbitrary type
  Exceptions : thrown if $module can not be instantiated
  Caller     : Bio::EnsEMBL::DBAdaptor

=cut

sub _get_adaptor {
  my( $self, $module, @args) = @_;

  my( $adaptor, $internal_name );

  #Create a private member variable name for the adaptor by replacing
  #:: with _

  $internal_name = $module;

  $internal_name =~ s/::/_/g;

  unless (defined $self->{'_adaptors'}->{$internal_name}) {
    eval "require $module";

    if($@) {
      warning("$module cannot be found.\nException $@\n");
      return undef;
    }

    $adaptor = "$module"->new($self, @args);

    $self->{'_adaptors'}->{$internal_name} = $adaptor;
  }

  return $self->{'_adaptors'}->{$internal_name};
}



=head2 add_db_adaptor

  Arg [1]    : string $name
               the name of the database to attach to this database
  Arg [2]    : Bio::EnsEMBL::DBSQL::DBConnection
               the db adaptor to attach to this database
  Example    : $db->add_db_adaptor('lite', $lite_db_adaptor);
  Description: Attaches another database instance to this database so 
               that it can be used in instances where it is required.
  Returntype : none
  Exceptions : none
  Caller     : EnsWeb

=cut

sub add_db_adaptor {
  my ($self, $name, $adaptor) = @_;

  unless($name && $adaptor && ref $adaptor) {
    throw('adaptor and name arguments are required');
  }

  #avoid circular references and memory leaks
  if($adaptor->isa('Bio::EnsEMBL::Container')) {
      $adaptor = $adaptor->_obj;
  }

  $self->{'_db_adaptors'}->{$name} = $adaptor;
}


=head2 remove_db_adaptor

  Arg [1]    : string $name
               the name of the database to detach from this database.
  Example    : $lite_db = $db->remove_db_adaptor('lite');
  Description: Detaches a database instance from this database and returns
               it.
  Returntype : none
  Exceptions : none
  Caller     : ?

=cut

sub remove_db_adaptor {
  my ($self, $name) = @_;

  my $adaptor = $self->{'_db_adaptors'}->{$name};
  delete $self->{'_db_adaptors'}->{$name};

  unless($adaptor) {
      return undef;
  }

  return $adaptor;
}


=head2 get_all_db_adaptors

  Arg [1]    : none
  Example    : @attached_dbs = values %{$db->get_all_db_adaptors()};
  Description: returns all of the attached databases as 
               a hash reference of key/value pairs where the keys are
               database names and the values are the attached databases  
  Returntype : hash reference with Bio::EnsEMBL::DBSQL::DBConnection values
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::ProxyAdaptor

=cut

sub get_all_db_adaptors {
  my ($self) = @_;

  unless(defined $self->{'_db_adaptors'}) {
    return {};
  }

  return $self->{'_db_adaptors'};
}



=head2 get_db_adaptor

  Arg [1]    : string $name
               the name of the attached database to retrieve
  Example    : $lite_db = $db->get_db_adaptor('lite');
  Description: returns an attached db adaptor of name $name or undef if
               no such attached database exists
  Returntype : Bio::EnsEMBL::DBSQL::DBConnection
  Exceptions : none
  Caller     : ?

=cut

sub get_db_adaptor {
  my ($self, $name) = @_;

  unless($self->{'_db_adaptors'}->{$name}) {
      return undef;
  }

  return $self->{'_db_adaptors'}->{$name};
}


#
# This method is automatically called by the Container DESTROY method to
# break circular references. Do not call it directly.
#
sub deleteObj {
  my $self = shift;

  #print STDERR "DBConnection::deleteObj : Breaking circular references:\n";

  if(exists($self->{'_adaptors'})) {
    foreach my $adaptor_name (keys %{$self->{'_adaptors'}}) {
      my $adaptor = $self->{'_adaptors'}->{$adaptor_name};

      #call each of the adaptor deleteObj methods
      if($adaptor && $adaptor->can('deleteObj')) {
        #print STDERR "\t\tdeleting adaptor\n";
        $adaptor->deleteObj();
      }

      #break dbadaptor -> object adaptor references
      delete $self->{'_adaptors'}->{$adaptor_name};
    }
  }

  #print STDERR "Cleaning up attached databases\n";

  #break dbadaptor -> dbadaptor references
  foreach my $db_name (keys %{$self->get_all_db_adaptors()}) {
    #print STDERR "\tbreaking reference to $db_name database\n";
    $self->remove_db_adaptor($db_name);
  }
}


=head2 DESTROY

  Arg [1]    : none
  Example    : none
  Description: Called automatically by garbage collector.  Should
               never be explicitly called.  The purpose of this destructor
               is to disconnect any active database connections.
  Returntype : none 
  Exceptions : none
  Caller     : Garbage Collector

=cut

sub DESTROY {
   my ($obj) = @_;

   #print STDERR "DESTROYING DBConnection\n";

   my $dbh = $obj->db_handle();

   if($dbh) {
     # don't disconnect if InactiveDestroy flag is set, this can really
     # screw up forked processes
     if(!$dbh->{'InactiveDestroy'}) {
       $dbh->disconnect();
     }
   }

   $obj->db_handle(undef);
}



1;
