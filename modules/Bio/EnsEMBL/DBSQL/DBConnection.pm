=Head1 NAME - Bio::EnsEMBL::DBSQL::DBConnection

=head1 SYNOPSIS

    $db = Bio::EnsEMBL::DBSQL::DBConnection->new(
        -user   => 'root',
        -dbname => 'pog',
        -host   => 'caldy',
        -driver => 'mysql',
        );


   You should use this as a base class for all objects (DBAdaptor) that connect to
   database. 

   $sth = $db->prepare( "SELECT something FROM yourtable" );

   If you go through prepare you could log all your select statements.
    

=head1 DESCRIPTION

  This only wraps around the perl DBI->connect call, 
  so you dont have to remember how to do this.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::DBConnection;

use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::Root;
use DBI;


@ISA = qw(Bio::EnsEMBL::Root);

sub new {
  my $class = shift;

  my $self = {};
  bless $self, $class;

  my (
      $db,
      $host,
      $driver,
      $user,
      $password,
      $port,
     ) = $self->_rearrange([qw(
			       DBNAME
			       HOST
			       DRIVER
			       USER
			       PASS
			       PORT
			      )],@_);
    

  $db   || $self->throw("Database object must have a database name");
  $user || $self->throw("Database object must have a user");

  if( ! $driver ) {
    $driver = 'mysql';
  }
  if( ! $host ) {
    $host = 'localhost';
  }
  if ( ! $port ) {
    $port = 3306;
  }

  my $dsn = "DBI:$driver:database=$db;host=$host;port=$port";

  my $dbh;
  eval{
    $dbh = DBI->connect("$dsn","$user",$password, {RaiseError => 1});
  };
    
  $dbh || $self->throw("Could not connect to database $db user " .
		       "$user using [$dsn] as a locator\n" . $DBI::errstr);

  $self->db_handle($dbh);

  $self->username( $user );
  $self->host( $host );
  $self->dbname( $db );
  $self->password( $password);

  return $self;
}


sub driver {
  my($self, $arg ) = @_;

  (defined $arg) &&
    ($self->{_driver} = $arg );
  return $self->{_driver};
}

sub port {
  my ($self, $arg) = @_;

  (defined $arg) && 
    ($self->{_port} = $arg );
  return $self->{_port};
}

sub dbname {
  my ($self, $arg ) = @_;
  ( defined $arg ) &&
    ( $self->{_dbname} = $arg );
  $self->{_dbname};
}

sub username {
  my ($self, $arg ) = @_;
  ( defined $arg ) &&
    ( $self->{_username} = $arg );
  $self->{_username};
}

sub host {
  my ($self, $arg ) = @_;
  ( defined $arg ) &&
    ( $self->{_host} = $arg );
  $self->{_host};
}

sub password {
  my ($self, $arg ) = @_;
  ( defined $arg ) &&
    ( $self->{_password} = $arg );
  $self->{_password};
}


=head2 _get_adaptor

  Title   : _get_adaptor
  Usage   : $obj->get_adaptor("full::module::name", )
  Returns : An already existing, or a new instance of the specified DB adaptor 
  Args : the fully qualified name of the adaptor module to retrieve
         (optional) special additional args for the adaptors constructor

=cut

sub _get_adaptor {
  my( $self, $module, @args) = @_;

  my( $adaptor, $internal_name );
  
  #Create a private member variable name for the adaptor by replacing
  #:: with _
  
  $internal_name = $module;

  $internal_name =~ s/::/_/g;

  unless (defined $self->{'_adaptors'}{$internal_name}) {
    eval "require $module";
    
    if($@) {
      $self->warn("$module cannot be found.\nException $@\n");
      return undef;
    }
      
    $adaptor = "$module"->new($self, @args);

    unless($adaptor->isa('Bio::EnsEMBL::DBSQL::BaseAdaptor')) {
      $self->throw("$module is not a Bio::EnsEMBL::DBSQL::BaseAdaptor\n");
      return undef;
    }

    $self->{'_adaptors'}{$internal_name} = $adaptor;
  }

  return $self->{'_adaptors'}{$internal_name};
}


=head2 db_handle

 Title   : _db_handle
 Usage   : $obj->_db_handle($newval)
 Function: 
 Example : 
 Returns : value of db_handle
 Args    : newvalue (optional)

=cut

sub db_handle {
   my ($self,$value) = @_;

   if( defined $value) {
      $self->{'_db_handle'} = $value;
    }
    return $self->{'_db_handle'};
}


=head2 prepare

 Title   : prepare
 Usage   : $sth = $dbobj->prepare("select seq_start,seq_end from feature where analysis = \" \" ");
 Function: prepares a SQL statement on the DBI handle

           If the debug level is greater than 10, provides information into the
           DummyStatement object
 Example :
 Returns : A DBI statement handle object
 Args    : a SQL string

=cut

sub prepare {
   my ($self,$string) = @_;

   if( ! $string ) {
       $self->throw("Attempting to prepare an empty SQL query.");
   }
   if( !defined $self->{_db_handle} ) {
      $self->throw("Database object has lost its database handle.");
   }
      
   # should we try to verify the string?
   return $self->{_db_handle}->prepare($string);
} 


=head2 DESTROY

 Title   : DESTROY
 Usage   : Called automatically by garbage collector
 Function: Disconnects any active database connections
 Example : -
 Returns : -
 Args    : -

=cut

sub DESTROY {
   my ($obj) = @_;

   if( $obj->{'_db_handle'} ) {
       $obj->{'_db_handle'}->disconnect;
       $obj->{'_db_handle'} = undef;
   }
}


