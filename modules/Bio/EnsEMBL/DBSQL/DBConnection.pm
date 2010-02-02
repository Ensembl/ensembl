=head1 LICENSE

  Copyright (c) 1999-2010 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::DBConnection

=head1 SYNOPSIS

  $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(
    -user   => 'anonymous',
    -dbname => 'homo_sapiens_core_20_34c',
    -host   => 'ensembldb.ensembl.org',
    -driver => 'mysql',
  );

  # SQL statements should be created/executed through this modules
  # prepare() and do() methods.

  $sth = $dbc->prepare("SELECT something FROM yourtable");

  $sth->execute();

  # do something with rows returned ...

  $sth->finish();

=head1 DESCRIPTION

This class is a wrapper around DBIs datbase handle.  It provides some
additional functionality such as the ability to automatically disconnect
when inactive and reconnect when needed.

Generally this class will be used through one of the object adaptors or
the Bio::EnsEMBL::Registry and will not be instantiated directly.

=head1 METHODS

=cut


package Bio::EnsEMBL::DBSQL::DBConnection;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Root;
use DBI;

use Bio::EnsEMBL::DBSQL::StatementHandle;

use Bio::EnsEMBL::Utils::Exception qw(deprecate throw info warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw(Bio::EnsEMBL::Root); # for backwards compatibility

my $reg = "Bio::EnsEMBL::Registry";

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
  Arg [WAIT_TIMEOUT]: (optional) integer
                 Time in seconds for the wait timeout to happen. Time after which
                 the connection is deleted if not used. By default this is 28800 (8 hours)
                 on most systems. 
                 So set this to greater than this if your connection are getting deleted.
                 Only set this if you are having problems and know what you are doing.

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
  Caller     : Bio::EnsEMBL::::Utils::ConfigRegistry ( for newer code using the registry)
               Bio::EnsEMBL::DBSQL::DBAdaptor        ( for old style code)
  Status     : Stable

=cut

sub new {
  my $class = shift;

  my (
    $db,                  $host,     $driver,
    $user,                $password, $port,
    $inactive_disconnect, $dbconn,   $wait_timeout
    )
    = rearrange( [
      'DBNAME', 'HOST', 'DRIVER', 'USER', 'PASS', 'PORT',
      'DISCONNECT_WHEN_INACTIVE', 'DBCONN', 'WAIT_TIMEOUT'
    ],
    @_
    );

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

    if($dbconn->disconnect_when_inactive()) {
      $self->disconnect_when_inactive(1);
    }
  } else {
    $db   || throw("-DBNAME argument is required.");
    $user || throw("-USER argument is required.");

    $driver ||= 'mysql';
    $host   ||= 'mysql';
    
    if(!defined($port)){
      $port   = 3306;
      if($host eq "ensembldb.ensembl.org"){
	if( $db =~ /\w+_\w+_\w+_(\d+)/){
	  if($1 >= 48){
	    $port = 5306;
	  }
	}
      }
    }

    $wait_timeout   ||= 0;

    $self->username( $user );
    $self->host( $host );
    $self->dbname( $db );
    $self->password( $password );
    $self->port($port);
    $self->driver($driver);
    $self->timeout($wait_timeout);

    if($inactive_disconnect) {
      $self->disconnect_when_inactive($inactive_disconnect);
    }
  }

#  if(defined $dnadb) {
#    $self->dnadb($dnadb);
#  }
  return $self;
}


=head2 connect

  Example    : $dbcon->connect()
  Description: Connects to the database using the connection attribute 
               information.
  Returntype : none
  Exceptions : none
  Caller     : new, db_handle
  Status     : Stable

=cut

sub connect {
  my ($self) = @_;

  if ( $self->connected() ) { return }

  $self->connected(1);

  if ( defined( $self->db_handle() ) and $self->db_handle()->ping() ) {
    warning( "unconnected db_handle is still pingable, "
        . "reseting connected boolean\n" );
  }

  my ( $dsn, $dbh );

  if ( $self->driver() eq "Oracle" ) {

    $dsn = "DBI:Oracle:";

    eval {
      $dbh =
        DBI->connect( $dsn,
        sprintf( "%s@%s", $self->username(), $self->dbname() ),
        $self->password(), { 'RaiseError' => 1, 'PrintError' => 0 } );
    };

  } elsif ( $self->driver() eq "ODBC" ) {

    $dsn = sprintf( "DBI:ODBC:%s", $self->dbname() );

    eval {
      $dbh = DBI->connect(
        $dsn,
        $self->username(),
        $self->password(),
        {
          'LongTruncOk'     => 1,
          'LongReadLen'     => 2**16 - 8,
          'RaiseError'      => 1,
          'PrintError'      => 0,
          'odbc_cursortype' => 2
        } );
    };

  } elsif ( $self->driver() eq "Sybase" ) {

    $dsn =
      sprintf( "DBI:Sybase:server=%s;database=%s;tdsLevel=CS_TDS_495",
      $self->host(), $self->dbname() );

    eval {
      $dbh = DBI->connect(
        $dsn,
        $self->username(),
        $self->password(),
        {
          'LongTruncOk' => 1,
          'RaiseError'  => 1,
          'PrintError'  => 0
        } );
    };

  } else {

    $dsn = sprintf(
      "DBI:%s:database=%s;host=%s;port=%s",
      $self->driver(), $self->dbname(),
      $self->host(),   $self->port() );

    eval {
      $dbh =
        DBI->connect( $dsn, $self->username(), $self->password(),
        { 'RaiseError' => 1 } );
    };
  }

  if ( !$dbh || $@ || !$dbh->ping() ) {
    warn( "Could not connect to database "
        . $self->dbname()
        . " as user "
        . $self->username()
        . " using [$dsn] as a locator:\n"
        . $DBI::errstr );

    $self->connected(0);

    throw("Could not connect to database "
        . $self->dbname()
        . " as user "
        . $self->username()
        . " using [$dsn] as a locator:\n"
        . $DBI::errstr );
  }

  $self->db_handle($dbh);

  if ( $self->timeout() ) {
    $dbh->do( "SET SESSION wait_timeout=" . $self->timeout() );
  }

  #print("CONNECT\n");
} ## end sub connect


=head2 connected

  Example    : $dbcon->connected()
  Description: Boolean which tells if DBConnection is connected or not.
               State is set internally, and external processes should not alter state.
  Returntype : undef or 1
  Exceptions : none
  Caller     : db_handle, connect, disconnect_if_idle, user processes
  Status     : Stable

=cut

sub connected {
  my $self = shift;

  # Use the process id ($$) as part of the key for the connected flag.
  # This forces the opening of another connection in a forked subprocess.
  $self->{'connected'.$$} = shift if(@_);
  return $self->{'connected'.$$};
}

sub disconnect_count {
  my $self = shift;
  return $self->{'disconnect_count'} = shift if(@_);
  $self->{'disconnect_count'}=0 unless(defined($self->{'disconnect_count'}));
  return $self->{'disconnect_count'};
}

sub timeout{
  my($self, $arg ) = @_;

  (defined $arg) &&
    ($self->{_timeout} = $arg );

  return $self->{_timeout};

}

sub query_count {
  my $self = shift;
  return $self->{'_query_count'} = shift if(@_);
  $self->{'_query_count'}=0 unless(defined($self->{'_query_count'}));
  return $self->{'_query_count'};
}

=head2 equals

  

=cut

  
sub equals {
  my ( $self, $dbc ) = @_;

  if ( $dbc->host() eq $self->host()
    && $dbc->dbname()   eq $self->dbname()
    && $dbc->driver()   eq $self->driver()
    && $dbc->port()     eq $self->port()
    && $dbc->username() eq $self->username() )
  {
    return 1;
  }

  return 0;
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
  Status     : Stable

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
  Status     : Stable

=cut

sub port {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'_port'} = $value;
  }

  return $self->{'_port'};
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
  Status     : Stable

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
  Status     : Stable

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
  Status     : Stable

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
  Description: Getter/Setter for the password of to use for this
               connection.  There is currently no point in setting
               this value after the connection has already been
               established by the constructor.
  Returntype : string
  Exceptions : none
  Caller     : new
  Status     : Stable

=cut

sub password {
  my ( $self, $arg ) = @_;

  if ( defined($arg) ) {
    # Use an anonymous subroutine that will return the password when
    # invoked.  This will prevent the password from being accidentally
    # displayed when using e.g. Data::Dumper on a structure containing
    # one of these objects.

    $self->{_password} = sub { $arg };
  }

  return ( ref( $self->{_password} ) && &{ $self->{_password} } ) || '';
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
  Status     : Stable

=cut

sub disconnect_when_inactive {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'disconnect_when_inactive'} = $value;
    if ($value) {
      $self->disconnect_if_idle();
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
  Status     : Stable

=cut


sub locator {
  my ($self) = @_;

  return sprintf(
    "%s/host=%s;port=%s;dbname=%s;user=%s;pass=%s",
    ref($self),      $self->host(),     $self->port(),
    $self->dbname(), $self->username(), $self->password() );
}


=head2 db_handle

  Arg [1]    : DBI Database Handle $value
  Example    : $dbh = $db_connection->db_handle() 
  Description: Getter / Setter for the Database handle used by this
               database connection.
  Returntype : DBI Database Handle
  Exceptions : none
  Caller     : new, DESTROY
  Status     : Stable

=cut

sub db_handle {
   my $self = shift;

   # Use the process id ($$) as part of the key for the database handle
   # this makes this object fork safe.  fork() does not makes copies
   # of the open socket which creates problems when one of the forked
   # processes disconnects,
   return $self->{'db_handle'.$$} = shift if(@_);
   return $self->{'db_handle'.$$} if($self->connected);

   $self->connect();
   return $self->{'db_handle'.$$};
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
  Status     : Stable

=cut

sub prepare {
   my ($self,@args) = @_;

   if( ! $args[0] ) {
     throw("Attempting to prepare an empty SQL query.");
   }

#   print STDERR  "SQL(".$self->dbname."):" . join(' ', @args) . "\n";

   my $sth = $self->db_handle->prepare(@args);

   # return an overridden statement handle that provides us with
   # the means to disconnect inactive statement handles automatically
   bless $sth, "Bio::EnsEMBL::DBSQL::StatementHandle";
   $sth->dbc($self);
   $sth->sql($args[0]);

   $self->query_count($self->query_count()+1);
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
  Status     : Stable

=cut

sub do {
   my ($self,$string) = @_;

   if( ! $string ) {
     throw("Attempting to do an empty SQL query.");
   }

   #info("SQL(".$self->dbname."):$string");

   my $result = $self->db_handle->do($string);

   # disconnect if the disconnect when inactive flag is set and
   # there are no active statement handles

   if($self->disconnect_when_inactive()) {
     $self->disconnect_if_idle();
   }
  
   $self->query_count($self->query_count()+1);
 
   return $result;
}



=head2 disconnect_if_idle

  Arg [1]    : none
  Example    : $dbc->disconnect_if_idle();
  Description: Disconnects from the database if there are no currently active
               statement handles. 
               It is called automatically by the DESTROY method of the
               Bio::EnsEMBL::DBSQL::SQL::StatementHandle if the
               disconect_when_inactive flag is set.
               Users may call it whenever they want to disconnect. Connection will
               reestablish on next access to db_handle()
  Returntype : 0,1
               1=problem trying to disconnect while a statement handle was still active
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::SQL::StatementHandle::DESTROY
               Bio::EnsEMBL::DBSQL::DBConnection::do
  Status     : Stable

=cut

sub disconnect_if_idle {
  my $self = shift;

  return 0 if(!$self->connected());
  my $db_handle = $self->db_handle();
  return 0 unless(defined($db_handle));

  #printf("disconnect_if_idle : kids=%d activekids=%d\n",
  #       $db_handle->{Kids}, $db_handle->{ActiveKids});

  #If InactiveDestroy is set, don't disconnect.
  #To comply with DBI specification
  return 0 if($db_handle->{InactiveDestroy});

  #If any statement handles are still active, don't allow disconnection
  #In this case it is being called before a query has been fully processed
  #either by not reading all rows of data returned, or not calling ->finish
  #on the statement handle.  Don't disconnect, send warning
  if($db_handle->{ActiveKids} != 0) {
     warn("Problem disconnect : kids=",$db_handle->{Kids},
            " activekids=",$db_handle->{ActiveKids},"\n");
     return 1;
  }
  
  $db_handle->disconnect();
  $self->connected(undef);
  $self->disconnect_count($self->disconnect_count()+1);
  #print("DISCONNECT\n");
  $self->db_handle(undef);
  return 0;
}


=head2 add_limit_clause

  Arg [1]    : string $sql
  Arg [2]    : int $max_number
  Example    : my $new_sql = $dbc->add_limit_clause($sql,$max_number);
  Description: Giving an SQL statement, it adds a limit clause, dependent on the database 
               (in MySQL, should add a LIMIT at the end, MSSQL uses a TOP clause)									 
  Returntype : String containing the new valid SQL statement        
  Exceptions : none
  Caller     : general
  Status     : at risk

=cut


sub add_limit_clause{
    my $self = shift;
    my $sql = shift;
    my $max_number = shift;

    my $new_sql = '';
    if ($self->driver eq 'mysql'){
	$new_sql = $sql . ' LIMIT ' . $max_number;
    }
    elsif ($self->driver eq 'odbc'){
	#need to get anything after the SELECT statement
	$sql =~ /select(.*)/i;
	$new_sql = 'SELECT TOP ' . $max_number . $1;
    }
    else{
	warning("Not possible to convert $sql to an unknow database driver: ", $self->driver, " no limit applied");
	$new_sql = $sql;
    }
    return $new_sql;
}


=head2 from_date_to_seconds

  Arg [1]    : date $date
  Example    : my $string = $dbc->from_date_to_seconds($date);
  Description: Giving a string representing a column of type date 
                applies the database function to convert to the number of seconds from 01-01-1970
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : at risk

=cut

sub from_date_to_seconds{
    my $self=  shift;
    my $column = shift;

    my $string;
    if ($self->driver eq 'mysql'){
	$string = "UNIX_TIMESTAMP($column)";
    }
    elsif ($self->driver eq 'odbc'){
	$string = "DATEDIFF(second,'JAN 1 1970',$column)";
    }
    else{
	warning("Not possible to convert $column due to an unknown database driver: ", $self->driver);
	return '';
    }    
    return $string;
}


=head2 from_seconds_to_date

  Arg [1]    : int $seconds
  Example    : my $string = $dbc->from_seconds_to_date($seconds);
  Description: Giving an int representing number of seconds
                applies the database function to convert to a date 
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : at risk

=cut

sub from_seconds_to_date{
    my $self=  shift;
    my $seconds = shift;

    my $string;
    if ($self->driver eq 'mysql'){
	if ($seconds){
	    $string = "from_unixtime( ".$seconds.")";
	}
	else{
	    $string = "\"0000-00-00 00:00:00\"";
	}
    }
    elsif ($self->driver eq 'odbc'){
	if ($seconds){
	    $string = "DATEDIFF(date,'JAN 1 1970',$seconds)";
	}
	else{
	    $string = "\"0000-00-00 00:00:00\"";
	}
    }
    else{
	warning("Not possible to convert $seconds due to an unknown database driver: ", $self->driver);
	return '';

    }    
    return $string;
}

####
#deprecated functions
####

=head2

   group is no longer available in DBConnection and should be accessed if needed
   from an adaptor.

=cut

sub group {
  my ($self, $arg ) = @_;
  ( defined $arg ) &&
    ( $self->{_group} = $arg );
  deprecate "group should not be called from DBConnection but from an adaptor\n";
  return $self->{_group};
}

=head2 species

   species is no longer available in DBConnection and should be accessed if needed
   from an adaptor.

=cut

sub species {
  my ($self, $arg ) = @_;
  ( defined $arg ) &&
    ( $self->{_species} = $arg );
  deprecate "species should not be called from DBConnection but from an adaptor\n";
  return $self->{_species};
}

1;
