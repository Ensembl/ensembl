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

  This only wraps around the perl DBI->connect call, so you dont have to remember
  how to do this.

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
			      )],@args);
    

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
  
  my( $dbh );
  eval{
    $dbh = DBI->connect("$dsn","$user",$password, {RaiseError => 1});
  };
  
  $dbh || $self->throw
    ( "Could not connect to database $db user $user using [$dsn] as a locator\n". 
     $DBI::errstr );
}



sub prepare {
}

sub host {
}

sub user {
}

sub dbname {
}

sub pass {
}

 
