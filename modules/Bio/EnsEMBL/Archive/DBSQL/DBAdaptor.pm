#
# BioPerl module for Bio::EnsEMBL::Archive::DBSQL::DBAdaptor
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Archive::DBSQL::DBAdaptor

Object representing an instance of the Archive DB

=head1 SYNOPSIS

    $db = Bio::EnsEMBL::Archive::DBSQL::DBAdaptor->new(
						  -user   => 'pig',
						  -dbname => 'pog',
						  -host   => 'pug',
						  -driver => 'mysql',
						  );

    $asad  = $db->get_VersionedSeqAdaptor();

    my @aseqs = $asad->fetch_by_object_id('ENSE00000456');

=head1 DESCRIPTION

This adaptor provides the connection to a database and a handle on other adaptors used in this database, such as ArchiveSeqAdaptor. 

The archive database holds a slice of data for older versions of proteins,
genes, and exons. The purpose of this database is to allow versioning in 
EnsEMBL, holding only the most recent of an entry in the main DBSQL database, and storing here only the relevant information of older versions.

=head1 CONTACT

email: elia@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Archive::DBSQL::DBAdaptor;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;
use DBI;
use Bio::EnsEMBL::DBSQL::SQL;
use Bio::EnsEMBL::Archive::DBSQL::VersionedSeqAdaptor;
use Bio::EnsEMBL::Archive::DBSQL::SeqAdaptor;

BEGIN {
    print STDERR "\n\n**** Loading Archive ******\n\n";
};

@ISA = qw(Bio::EnsEMBL::Root);

sub new {
  my($class, @args) = @_;
  
  my $self = {};
  bless $self,$class;
  
  my ($db,$host,$driver,$user,
      $password,$port,
      $readonly,$dnadb) = $self->_rearrange([qw(
					 DBNAME
					 HOST
					 DRIVER
					 USER
					 PASS
					 PORT
					 READONLY
					 DNADB
					 )],@args);
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
	
  my $dbh = DBI->connect("$dsn","$user",$password, {RaiseError => 1});
  
  $dbh || $self->throw("Could not connect to database $db user $user using [$dsn] as a locator");

  $self->_readonly($readonly);
  if($self->_readonly){
      $self->warn("Archive Database accessed in READONLY mode");
  }


  $self->_db_handle($dbh);
  $self->username( $user );
  $self->host( $host );
  $self->dbname( $db );
  $self->password( $password);
  $self->dnadb($dnadb);
  return $self; # success - we hope!
}

=head2 get_SeqAdaptor

 Title   : get_SeqAdaptor
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_SeqAdaptor {
    my ($self) = @_;
    
    return Bio::EnsEMBL::Archive::DBSQL::SeqAdaptor->new($self);
}

=head2 get_VersionedSeqAdaptor

 Title   : get_VersionedSeqAdaptor
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_VersionedSeqAdaptor {
    my ($self) = @_;
    
    return Bio::EnsEMBL::Archive::DBSQL::VersionedSeqAdaptor->new($self);
}

=head2

 Title   : get_new_stable_ids
 Usage   : @ids = $obj->get_new_stable_ids('exon',6);
 Function: 
 Example : 
 Returns : new id, ENS-style, as an array
 Args    : id type, one of 'exon','transcript','gene','translation'


=cut

my %stubhash = ( 'exon' => 'ENSE', 'transcript' => 'ENST','gene' => 'ENSG','translation' => 'ENSP' );

sub get_new_stable_ids {
    my ($self,$table,$number) = @_;

    if( !defined $number ) {
	$self->throw("Must call as table,number");
    }

    if( !defined $stubhash{$table} ) {
	$self->throw("Does not have $table as a valid stable id table");
    }
    
    my $stub= $stubhash{$table};
    $table .= "_stable";
    my @out;

    my $query="lock table $table write";

    $self->_execute($query);

    # wrap critical region in an eval so we can catch errors and release table

    eval {

	my $query = "select external_id as id from $table";
	
	my $sth   = $self->prepare($query);
	my $res   = $sth->execute;
	my $row   = $sth->fetchrow_hashref;
	my $id    = $row->{id};
	my $new=0;
	if (!defined($id) || $id eq "") {
	    $id = $stub . "00000000000";
	    $new=1;
	}
	
	if ($id =~ /\D+(\d+)$/) {
	    
	    my $newid  = $1;
	    my $i;
	    
	    foreach $i ( 1..$number ) {

		$newid++;
		
		
		if (length($newid) > 11) {
		    if ($newid =~ /^0/) {
			$newid =~ s/^0//;
		    } else {
			$self->throw("Can't truncate number string to generate new id [$newid]");
		    }
		}
		my $c = $stub . $newid;
		my $query;
		if (! $new) {
		    $query = "update $table set external_id = '$c', created=now()";
		}
		else {
		    $query = "insert into $table (internal_id,external_id,created) values (NULL,'$c',NOW())";
		}
		
		$self->_execute($query);
		
		push(@out,$c);
	    }
	    
	    
	} else {
	    $self->throw("[$id] does not look like an object id (e.g. ENST00000019784)");
	}
    };
    
    my $error = undef;

    if( $@ ) {
	$error = $@;
    }
    
    $query = "unlock tables";
    $self->_execute($query);

    if( defined $error ) {
	$self->throw("Problem in making IDs. Unlocked tables. \n\n Error $@");
    }
    
    return @out;
    
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
       $self->throw("Attempting to prepare an empty SQL query!");
   }
   if( !defined $self->_db_handle ) {
      $self->throw("Database object has lost its database handle! getting otta here!");
   }

   return $self->_db_handle->prepare($string);
}

=head2 _execute

 Title   : _execute
 Usage   :
 Function: Internal SQL prepare and execute function
           which does nothing if in readonly mode
 Example :
 Returns : 
 Args    :


=cut

sub execute{
   my ($self,$query) = @_;

   my $sth;
   if (($query =~ /insert|update/) && ($self->_readonly)){
       #print "READONLY: $query\n";
   }else{
       $sth = $self->prepare($query);
       $sth->execute;
   }
   return $sth;
}

=head2 _readonly

 Title   : _readonly
 Usage   : $obj->_readonly($newval)
 Function: 
 Example : 
 Returns : value of _readonly
 Args    : newvalue (optional)


=cut

sub _readonly{
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'_readonly'} = $value;
    }
    return $self->{'_readonly'};
}

=head2 _db_handle

 Title   : _db_handle
 Usage   : $obj->_db_handle($newval)
 Function: 
 Example : 
 Returns : value of _db_handle
 Args    : newvalue (optional)


=cut

sub _db_handle{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_db_handle'} = $value;
    }
    return $self->{'_db_handle'};

}

=head2 _lock_tables

 Title   : _lock_tables
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _lock_tables{
   my ($self,@tables) = @_;
   
   my $state;
   foreach my $table ( @tables ) {
       if( $self->{'_lock_table_hash'}->{$table} == 1 ) {
	   $self->warn("$table already locked. Relock request ignored");
       } else {
	   if( $state ) { $state .= ","; } 
	   $state .= "$table write";
	   $self->{'_lock_table_hash'}->{$table} = 1;
       }
   }

   my $sth = $self->prepare("lock tables $state");
   my $rv = $sth->execute();
   $self->throw("Failed to lock tables $state") unless $rv;

}

=head2 _unlock_tables

 Title   : _unlock_tables
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _unlock_tables{
   my ($self,@tables) = @_;

   my $sth = $self->prepare("unlock tables");
   my $rv  = $sth->execute();
   $self->throw("Failed to unlock tables") unless $rv;
   %{$self->{'_lock_table_hash'}} = ();
}


=head2 DESTROY

 Title   : DESTROY
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub DESTROY {
   my ($obj) = @_;

   #$obj->_unlock_tables();

   if( $obj->{'_db_handle'} ) {
       $obj->{'_db_handle'}->disconnect;
       $obj->{'_db_handle'} = undef;
   }
}

=head2 deleteObj

    Title   : deleteObj
    Usage   : $dbObj->deleteObj
    Function: Call when you are done with this object. Breaks links between objects. Necessary to clean up memory.
    Example : -
    Returns : -
    Args    : -


=cut

sub deleteObj {

  my  $self=shift;
  my $dummy;

  print STDERR "Destroying DB Obj!\n";       
  $self->DESTROY;
  
  foreach my $name ( keys %{$self} ) {
    eval {
      $dummy = $self->{$name}; 
      $self->{$name}  = undef;
      $dummy->deleteObj;
    };
  }
}

=head2 dnadb

 Title   : dnadb
 Usage   : my $dnadb = $db->dnadb;
 Function: returns the database adaptor where the dna lives
           Useful if you only want to keep one copy of the dna
           on disk but have other databases with genes and features in
 Returns : dna database adaptor
  Args    : Bio::EnsEMBL::Archive::DBSQL::DBAdaptor

=cut

sub dnadb {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    if (! $arg->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")) {
      $self->throw("[$arg] is not a Bio::EnsEMBL::DBSQL::DBAdaptor");
    }
    $self->{_dnadb} = $arg;
  }
  return $self->{_dnadb} || $self;
}


1;
