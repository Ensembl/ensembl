
#
# BioPerl module for Obj
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Obj - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package EnsEMBL::Obj;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;


@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;
  print "Got",join(',',@args),"\n";
  my ($db,$host,$driver,$user,$password,$debug) = 
      $self->_rearrange([qw(DBNAME
			    HOST
			    DRIVER
			    USER
			    PASS
			    DEBUG
			    )],@args);
  print "Got $db as db and $user as user\n";
  
  $db || $self->throw("Database object must have a database name");
  $user || $self->throw("Database object must have a user");
  
  if( $debug ) {
      $self->_debug($debug);
  } else {
      $self->_debug(0);
  }
  
  if( ! $driver ) {
      $driver = 'mysql';
  }
  if( ! $host ) {
      $host = 'localhost';
  }
  my $dsn = "DBI:$driver:database=$db;host=$host";
  
  if( $debug && $debug > 10 ) {
      $self->_db_handle("dummy dbh handle in debug mode $debug");
  } else {
      
      my $dbh = DBI->connect("$dsn","$user","$password",{RaiseError => 1});
      $dbh || $self->throw("Could not connect to database $db user $user using [$dsn] as a locator");
      
      if( $self->_debug > 3 ) {
	  $self->warn("Using connection $dbh");
      }
      
      $self->_db_handle($dbh);
  }
  
# set stuff in self from @args
 return $make; # success - we hope!
}

=head2 write_Bug

 Title   : write_Bug
 Usage   : $db->write_Bug ($Bug)
 Function: Writes a bug in the bugtrack database
 Example : $db->write_Bug ($bug)
 Returns : nothing
 Args    : Bug object, 

=cut

sub write_Bug{
   my ($self,$bug) = @_;
   
   if( ! $bug->isa("Bio::EnsEMBL::BUGTRACK::Bug") ) {
       $self->throw("$bug is not a Bio::EnsEMBL::BUGTRACK::Bug object!");
   }
   $bug->id || $self->throw("Attempting to write a bug without a bug id!");
   $bug->title || $self->throw("Attempting to write a bug without a bug title!");
   $bug->type || $self->throw("Attempting to write a bug without a bug type!");

   my $sth = $self->prepare("insert into bug (id,title,type) values ('".$bug->id()."','".$bug->title."','".$bug->type."')");
   $sth->execute();
   
   foreach my $note ($bug->each_Worknote()) {
       $note->author || $self->throw("Attempting to write a note without a note author!");
       $note->note || $self->throw("Attempting to write a note without any note text!");

       my $sth2 = $self->prepare("insert into worknote (bug,author,date,note) values (LAST_INSERT_ID(),'".$note->author."',now(),'".$note->note."')");
       $sth2->execute();
   }
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

sub prepare{
   my ($self,$string) = @_;

   if( ! $string ) {
       $self->throw("Attempting to prepare an empty SQL query!");
   }

   if( $self->_debug > 10 ) {
       print STDERR "Prepared statement $string\n";
       my $st = Bio::EnsEMBL::DBSQL::DummyStatement->new();
       $st->_fileh(\*STDERR);
       $st->_statement($string);
       return $st;
   }

   # should we try to verify the string?

   return $self->_db_handle->prepare($string);
}

=head2 _debug

 Title   : _debug
 Usage   : $obj->_debug($newval)
 Function: 
 Example : 
 Returns : value of _debug
 Args    : newvalue (optional)


=cut

sub _debug{
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'_debug'} = $value;
    }
    return $self->{'_debug'};
    
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

=head2 DESTROY

 Title   : DESTROY
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub DESTROY{
   my ($obj) = @_;

   if( $obj->{'_db_handle'} ) {
       $obj->{'_db_handle'}->disconnect;
       $obj->{'_db_handle'} = undef;
   }
}
