
#
# Ensembl module for Bio::EnsEMBL::DBSQL::CrossMatchDBAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::CrossMatchDBAdaptor - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::CrossMatchDBAdaptor;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::RootI
use DBI;
use Bio::Root::RootI;
use Bio::EnsEMBL::DBSQL::SymmetricContigFeatureContainer;

@ISA = qw(Bio::Root::RootI);

sub new {
    my($class,@args) = @_;
 
   
    my $self = {};
    bless $self,$class;

    my ($db,$host,$driver,$user,$password,$port) = 
	$self->_rearrange([qw(DBNAME
			      HOST
			      DRIVER
			      USER
			      PASS
			      PORT
			      )],@args);
    $db   || $self->throw("Database object must have a database name");
    $user || $self->throw("Database object must have a user");
    
    if( ! $driver ) {
	$driver = 'mysql';
    }
    
    if( ! $host ) {
	$host = 'localhost';
    }
    my $dsn = "DBI:$driver:database=$db;host=$host";
    my $dbh = DBI->connect("$dsn","$user",$password, {RaiseError => 1});

    $self->_db_handle($dbh);

    return $self;
}


=head2 new_dbobj

 Title   : new_dbobj
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub new_dbobj{
   my ($self,@args) = @_;

   my $t = $self->_new_dbobj;
   if( defined $t ) { return $t; }

   # yank it out ;)

   my $sth = $self->prepare("select newdatabase from dblocation");
   $sth->execute;
   my ($loc) = $sth->fetchrow_array;
   
   #print STDERR "New database locator: $loc\n";
   my $db = Bio::EnsEMBL::DBLoader->new($loc);

   $self->_new_dbobj($db);

   return $db;

}


=head2 old_dbobj

 Title   : old_dbobj
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub old_dbobj{
   my ($self,@args) = @_;


   my $t = $self->_old_dbobj;
   if( defined $t ) { return $t; }

   # yank it out ;)

   my $sth = $self->prepare("select olddatabase from dblocation");
   $sth->execute;
   my ($loc) = $sth->fetchrow_array;
   #print STDERR "Old database locator: $loc\n";
   my $db = Bio::EnsEMBL::DBLoader->new($loc);

   $self->_old_dbobj($db);

   return $db;

}

=head2 get_SymmetricContigFeatureContainer

 Title   : get_SymmetricContigFeatureContainer
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_SymmetricContigFeatureContainer{
   my ($self,@args) = @_;

   return Bio::EnsEMBL::DBSQL::SymmetricContigFeatureContainer->new($self);
}


=head2 _new_dbobj

 Title   : _new_dbobj
 Usage   : $obj->_new_dbobj($newval)
 Function: 
 Returns : value of _new_dbobj
 Args    : newvalue (optional)


=cut

sub _new_dbobj{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_new_dbobj'} = $value;
    }
    return $obj->{'_new_dbobj'};

}

=head2 _old_dbobj

 Title   : _old_dbobj
 Usage   : $obj->_old_dbobj($newval)
 Function: 
 Returns : value of _old_dbobj
 Args    : newvalue (optional)


=cut

sub _old_dbobj{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_old_dbobj'} = $value;
    }
    return $obj->{'_old_dbobj'};

}





=head2 prepare

 Title   : prepare
 Usage   : $sth = $dbobj->prepare("select seq_start,seq_end from feature where analysis = \" \" ");
 Function: prepares a SQL statement on the DBI handle

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

   # should we try to verify the string?

   return $self->_db_handle->prepare($string);
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
