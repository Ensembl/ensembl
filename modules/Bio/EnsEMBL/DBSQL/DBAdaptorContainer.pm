=Head1 NAME - Bio::EnsEMBL::DBSQL::DBAdaptorContainer

=head1 SYNOPSIS

$db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -user   => 'root',
        -dbname => 'pog',
        -host   => 'caldy',
        -driver => 'mysql'
        );

$db_adaptor_container = new Bio::EnsEMBL::DBSQL::DBAdaptorContainer($db);
    
=head1 DESCRIPTION

This object is a hack necessary to work around perls circular reference 
memory leak problems.  Its sole purpose is to channel calls to the 
dbadaptor that it holds on to, and to break all circular memory references
to the DBAdaptor upon this objects destruction. Sneaky DBAdaptors should
return a DBAdaptorContainer object which contains a DBAdaptor.  

=head1 CONTACT

Arne Stabenau - stabenau@ebi.ac.uk
Graham McVicker - mcvicker@ebi.ac.uk
Ewan Birney - birney@ebi.ac.uk 

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


use strict;

package Bio::EnsEMBL::DBSQL::DBAdaptorContainer;

use Bio::EnsEMBL::Container;

use vars ('@ISA', '$AUTOLOAD');

@ISA = qw(Bio::EnsEMBL::Container);


=head2 add_db_adaptor

  Arg [1]    : Bio::EnsEMBL::DBHolder $dbholder
  Example    : none
  Description: Overrides the DBAdaptor method to ensure that DBAdaptors and
               not DBAdaptorHolders are attached together
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub add_db_adaptor {
  my ($self, $name, $dbholder) = @_;
  
  # Check if they actually passed in a dbholder.  
  # We want to attach the contained DBAdaptor not the DBHolder!
  if($dbholder->isa('Bio::EnsEMBL::DBSQL::DBAdaptorContainer')) {
    return $self->_obj->add_db_adaptor($name, $dbholder->_obj);
  }

  #the passed in object was not a dbadpatorholder, 
  #assume it is a normal DBAdaptor
  $self->_dba->add_db_adaptor($name, $dbholder);
}


1;
