
#
# Ensembl module for Bio::EnsEMBL::DBOLD::CrossMatchDBAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBOLD::CrossMatchDBAdaptor - DESCRIPTION of Object

=head1 SYNOPSIS

my $db=Bio::EnsEMBL::DBOLD::Obj->new(-dbname=>'july_dna',-host=>'ecs1c',-user=>'ensadmin');

my $cross=Bio::EnsEMBL::DBOLD::CrossMatchDBAdaptor->new(-dbname=>'crossmatch',-host=>'ecs1c',-user=>'ensadmin');

$db->_crossdb($cross);

=head1 DESCRIPTION

This Object is a database adapter for the crossmatch database, it loads the old and new databases, which are held in _new_db and _old_db. It also gets a SymmetricFeatureContainer object, where the methods for getting and returning crossmatches are.

The crossdb can then be added to a standard db in its _crossdb method.

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBOLD::CrossMatchDBAdaptor;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::EnsEMBL::Root
use DBI;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::DBOLD::SymmetricContigFeatureContainer;
use Bio::EnsEMBL::DBLoader;

@ISA = qw(Bio::EnsEMBL::Root);

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
 Usage   : $crossdb->new_dbobj
 Function: reads the dblocation table to load the new db
 Example : $crossdb->new_dbobj()
 Returns : Bio::EnsEMBL::DBOLD::Obj
 Args    : none


=cut

sub new_dbobj{
   my ($self) = @_;

   my $t = $self->_new_dbobj;
   if( defined $t ) { return $t; }

   # yank it out ;)

   my $sth = $self->prepare("select newdatabase from dblocation");
   $sth->execute;
   my ($loc) = $sth->fetchrow_array;
   
   #print STDERR "New database locator: $loc\n";
   my $db = Bio::EnsEMBL::DBLoader->new($loc);

   $self->_new_dbobj($db);
   $db->_crossdb($self);
   return $db;
}


=head2 old_dbobj

 Title   : old_dbobj
 Usage   : $crossdb->old_dbobj
 Function: reads the dblocation table to load the old db
 Example : $crossdb->old_dbobj()
 Returns : Bio::EnsEMBL::DBOLD::Obj
 Args    : none


=cut

sub old_dbobj{
   my ($self) = @_;

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
 Usage   : $crossdb->get_SymmetricContigFeatureContainer
 Function: Gets a Bio::EnsEMBL::DBOLD::SymmetricContigFeatureContainer
 Example : $crossdb->get_SymmetricContigFeatureContainer
 Returns : Bio::EnsEMBL::DBOLD::SymmetricContigFeatureContainer
 Args    : none


=cut

sub get_SymmetricContigFeatureContainer{
   my ($self) = @_;

   return Bio::EnsEMBL::DBOLD::SymmetricContigFeatureContainer->new($self);
}

=head2 get_clonelist

 Title   : get_clonelist
 Usage   : $crossdb->get_clonelist
 Function: Reads the clonelist table, to return list of clones with different
           versions between old and new
 Example : $crossdb->get_clonelist
 Returns : array of strings 
 Args    : none

=cut

sub get_clonelist{
   my ($self) = @_;

   my @clones;
   my $sth=$self->prepare("select clone from clonelist");
   $sth->execute;
   while (my $clone = $sth->fetchrow_array()) {
       push @clones,$clone;
   }
   return @clones;
}

=head2 _new_dbobj

 Title   : _new_dbobj
 Usage   : $obj->_new_dbobj($newval)
 Function: get/set for the new db adapter
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
 Function: get/set for the old db adapter
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
 Usage   : $sth = $dbobj->prepare($statement);
 Function: prepares a SQL statement on the DBI handle
 Example :$sth = $dbobj->prepare("select seq_start,seq_end from feature where analysis = \'example\' ");
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
 Function: get/set for the db handle
 Example : $obj->_db_handle($newval)
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

