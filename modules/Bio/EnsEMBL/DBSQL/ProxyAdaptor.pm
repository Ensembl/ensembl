=Head1 NAME - Bio::EnsEMBL::DBSQL::ProxyAdaptor

=head1 SYNOPSIS

   $primary_adaptor = $db_adaptor->get_ObjectAdaptor;

   #the lite databases ObjectAdaptor is the primary adaptor for this proxy
   $proxy_adaptor = new ProxyAdaptor($db_adaptor, $primary_adaptor);

   #attach some other databases to the database_adaptor we are using
   $db_adaptor->add_db_adaptor('lite', $lite_db_adaptor);
   $db_adpator->add_db_adaptor('snp', $snp_db_adaptor);
   
   #the proxy will first try the primary adaptor for the fetch_by_Slice method
   #if it cannot be found then it will try the lite and snp db adaptors
   #before giving up
   my $object = $proxy_adaptor->fetch_by_Slice($slice);

=head1 DESCRIPTION

   This acts as Proxy for an adaptor class.  It fills in for a database 
   specific class and when a request for a particular method is called it 
   will make an decision as to which adaptor to use.  It will first try
   the primary adaptor which is passed in as a constructor argument and
   will then proceed to try other adaptors attached to the database adaptor
   before giving up and throwing an exception.

   If intelligent decisions are required for certain methods (e.g. a certain
   adaptor other than the primary adaptor should be used for certain requests)
   then this ProxyAdaptor should be extended, and the methods defined in the
   subclass.  

=head1 CONTACT

  Graham McVicker: mcvicker@ebi.ac.uk
  Arne Stabenau  : stabenau@ebi.ac.uk
  Ewan Birney    : birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::DBSQL::ProxyAdaptor;

use strict;
use vars qw($AUTOLOAD @ISA);
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 new

  Arg [1]    : Bio::EnsEMBL::DBConnection $db
               The database adaptor which created this adaptor.
  Arg [2]    : The primary adaptor which should be used first for
               the invocation of requests (unless methods are defined
               by a subclass). I.e. this is the default adaptor to try.
  Example    : $proxy = new Bio::EnsEMBL::DBSQL::ProxyAdaptor($db,$primary_ad);
  Description: Constructor for a ProxyAdaptor - creates a new adaptor. 
  Returntype : Bio::EnsEMBL::DBSQL::ProxyAdaptor
  Exceptions : thrown if the $primary_adaptor argument is not defined
  Caller     : Bio::EnsEMBL::DBSQL::DBAdaptor, subclasses

=cut

sub new {
  my ($class, $db, $primary_adaptor) = @_;

  #invoke superclass constructor
  my $self = $class->SUPER::new($db);

  unless($primary_adaptor) {
    throw("The primary_adaptor argument is required\n");
    return undef;
  }
  
  #determine the type of adaptor the proxy is filling in for
  $self->{'_proxy_type'} = ref($primary_adaptor);
  
  #strip out fully qualified package name
  $self->{'_proxy_type'} =~ s/.*:://;

  $self->{'_primary_adaptor'} = $primary_adaptor;

  return $self;
}



=head2 AUTOLOAD

  Arg [1]    : list of arbitrary values @args
               a list of arguments to pass to the request method
  Example    : -
  Description: AUTOLOAD method should not be called directly.  It is called
               implicitly when a method requested from this class cannot be
               found. This method first tries to execute the requested method
               in the primary adaptor.  If the method cannot be found then 
               it searches the other attached databases for equivalent adaptors
               and tries then one at a time.
  Returntype : arbitrary
  Exceptions : thrown if the requested method cannot be found on the primary
               adaptor or on any of the attached databases.
  Caller     : called implicitly by perl

=cut

sub AUTOLOAD {
  my ($self, @args) =  @_;

  #determine the method which was called
  my $method = $AUTOLOAD;
  
  #strip out fully qualified method name
  $method =~ s/.*:://;
  
  #
  # First try the primary adaptor
  #
  if($self->{'_primary_adaptor'}->can($method)) {
    #execute the request using the primary adaptor
    my $adaptor = $self->{'_primary_adaptor'};

    return $adaptor->$method(@args);
  } 

  #
  # The request could not be filled by the primary adaptor
  # try the same request using all of the attached databases
  #
  my @databases = values %{$self->db()->get_all_db_adaptors()};
  foreach my $database (@databases) {

    #Try to get the appropriate adaptor from the database
    my $get_adaptor = "get_" . $self->{'_proxy_type'};
    if($database->can($get_adaptor)) {

      #Try to invoke the request on the database's adaptor
      my $adaptor = eval "\$database->$get_adaptor";
      if($adaptor->can($method)) {
	return $adaptor->$method(@args);
      }
    }
  }

  #none of the attached adaptors could fulfill the request either
  throw("The requested method $method could not be found in the " 
        . $self->{'_proxy_type'} . " of the attached databases:" .
	       @databases);
  return undef;
}


=head2 deleteObj

  Args       : none
  Example    : none
  Description: breaks circular references and is recursivley called during 
               memory cleanup (hopefully)
  Returntype : none
  Exceptions : none
  Caller     : DBConnection->deleteObj()

=cut

sub deleteObj {
  my $self = shift;
  delete $self->{'_primary_adaptor'};
  $self->SUPER::deleteObj();
}
 



=head2 DESTROY

  Arg [1]    : none
  Example    : none
  Description: Called automatically by garbage collector. This method does not
               actually do anything, but if it was not defined, then the 
               AUTOLOAD method would be called during destruction by the 
               garbage collector instead.
  Returntype : none
  Exceptions : none
  Caller     : automatic

=cut

sub DESTROY {
  #do nothing
}


  
1;
