
=Head1 NAME - Bio::EnsEMBL::Lite::DBAdaptor

=head1 SYNOPSIS

    $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -user   => 'root',
        -dbname => 'pog',
        -host   => 'caldy',
        -driver => 'mysql',
        );


    $db->get_GeneAdaptor

=head1 DESCRIPTION

This DBAdaptor provides the link to the lite database. It allows the rapid creation of drawable objects, mainly of gene objects. The created Gene objects are not entirely complete, but they are connected to the DBSQL::GeneAdaptor to allow lazy loading of missing data. 

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Lite::DBAdaptor;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::DBConnection;

@ISA = qw(Bio::EnsEMBL::DBSQL::DBConnection);


#Override constructor inherited by Bio::EnsEMBL::DBSQL::DBConnection
sub new {
  my($class, @args) = @_;

  #call superclass constructor
  my $self = $class->SUPER::new(@args);
  
  return $self;
}

sub get_GeneAdaptor {
  my $self = shift;

    return 
      $self->_get_adaptor("Bio::EnsEMBL::Lite::GeneAdaptor");
}

sub get_LandmarkMarkerAdaptor {
  my $self = shift;

    return 
      $self->_get_adaptor("Bio::EnsEMBL::Lite::LandmarkMarkerAdaptor");
}

=head2 core_DBAdaptor

  Arg  [1]  : Bio::EnsEMBL::DBSQL::DBAdaptor $db
              The link to the core db
  Function  : Makes the link to coredb available. Objects can attach
              with this to the core db Adaptors for lazy loading
  Returntype: Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions: none
  Caller    : GeneAdaptor->fetch_by_Slice()

=cut

sub core_DBAdaptor {
  my($self, $arg ) = @_;

  (defined $arg) &&
    ($self->{_coreDBAdaptor} = $arg );
  return $self->{_coreDBAdaptor};
}



=head2 _get_adaptor

  Arg   1   : string $full_module_name
             "Bio::EnsEMBL::Lite::GeneAdaptor"
  Function  : Requires and creates Modules. Makes singletons of these,
              so it keeps a cache for each different module.
  Returntype: $full_module_name
  Caller    : intern

=cut

sub _get_adaptor {
  my( $self, $module) = @_;

  my( $adaptor, $internal_name );
  
  #Create a private member variable name for the adaptor by: 
  # (1)stripping out the the parent package names, 
  # (2)prepending an underscore, and (3) converting to all lower case 
  
  $module =~ /(::)?(\S+)$/;
  $internal_name = '_' . lc($2); 

  unless ($adaptor = $self->{'_int_adaptors'}{$internal_name}) {
    eval "require $module";
    
    if($@) {
      $self->warn("$module cannot be found.\nException $@\n");
      return undef;
    }
      
    $adaptor = "$module"->new($self);
    $self->{'_int_adaptors'}{$internal_name} = $adaptor;
  }

  return $adaptor;
}
