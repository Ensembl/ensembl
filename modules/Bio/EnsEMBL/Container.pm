
=head1 NAME - Bio::EnsEMBL::DBSQL::DBAdaptorHolder

=head1 SYNOPSIS

$container = new Bio::EnsEMBL::Container($obj);

=head1 DESCRIPTION

This object is a hack necessary to work around perls circular reference 
memory leak problems.  Its sole purpose is to channel calls to the 
object which is held onto by the container and to invoke the objects deleteObj
method to breaks all circular memory references at the correct time.

Eventually this problem may be solved through the use of WeakRef instead. 

=head1 CONTACT

Post questions to the EnsEMBL developer mailing list: <ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


use strict;

package Bio::EnsEMBL::Container;

use vars ('@ISA', '$AUTOLOAD');
use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Root);


=head2 new

  Arg [1]    : Bio::EnsEMBL::DBAdaptor $dba
               The dbadaptor to wrap this holder around.
  Example    : $dba_holder = new Bio::EnsEMBL::DBSQL::DBAdaptorHolder($dba);
  Description: Creates a new DBAdaptor holder object that forwards calls to 
               $dba and breaks circular references to $dba upon destruction. 
  Returntype : Bio::EnsEMBL::DBAdaptorHolder
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBAdaptor

=cut

sub new {
  my ($class, $object) = @_;

  my $self = $class->SUPER::new();

  unless($object) {
    $self->throw("object argument is required");
  }

  $self->_obj($object);

  return $self;
}


=head2 _obj

  Arg [1]    : Generic Object $obj (optional)
  Example    : $object = $self->_obj;
  Description: PRIVATE Getter/Setter for the object held by this container.
  Returntype : Generic Object
  Exceptions : none
  Caller     : internal

=cut

sub _obj {
  my ($self, $obj) = @_;

  if($obj) {
    $self->{_obj} = $obj;
  }

  return $self->{_obj};
}


=head2 isa

  Arg [1]    : string $module
  Example    : none
  Description: Overrides the base perl object isa so that this object is
               also considered to be the contained object. Very sneaky
               and a bit of a hack, but necessary to make this container
               completely transparent.  
  Returntype : boolean
  Exceptions : none
  Caller     : general

=cut

sub isa {
  my ($self, $module) = @_;
  
  if($module eq ref $self) {
    return 1;
  } 
  
  if($self->_obj()->isa($module)) {
    return 1;
  }

  return $self->SUPER::isa($module);
}

=head2 can

  Arg [1]    : string $method
  Example    : none
  Description: Like the isa method above, calls can on the contained
               _obj.  Without this can fails because it gets trapped
               by the AUTOLOAD method.
  Returntype : boolean
  Exceptions : none
  Caller     : general

=cut


sub can {
    my( $self, $method ) = @_;
    
    return $self->_obj->can($method);
}

=head2 AUTOLOAD

  Arg [1]    : @args list of arguments
  Example    : none
  Description: Automatically called to forward calls to the object in this 
               container
  Returntype : arbitrary
  Exceptions : none
  Caller     : perl

=cut


sub AUTOLOAD {
  my ($self, @args) = @_;

  my $method = $AUTOLOAD;
  $method =~ s/.*:://;
  
  # call the method on the contained object
  if ($self->_obj->can($method)) {
    # update the symbol table so AUTOLOAD is not called
    # the next time this method is called (faster this way)
    no strict 'refs';
    *{$AUTOLOAD} = 
      sub {
        my ($self, @args) = @_;  
        return $self->_obj->$method(@args);
      };
    
    return $self->_obj->$method(@args);
  } else {
    # Method does not exist
    $self->throw("method '$method' does not exist in '". ref($self->_obj) ."'");
  }
}


=head2 DESTROY

  Arg [1]    : none
  Example    : none
  Description: Automatically called when there are no more references to 
               this container.  This container calls deleteObj on the 
               contained object to break the circular
               references contained within the object which would otherwise
               prevent the garbage collection of the object and objects 
               referenced by the object.
  Returntype : none
  Exceptions : none
  Caller     : perl

=cut

sub DESTROY {
  my $self = shift;

 # print STDERR "Container::DESTROY : Breaking circular references:\n";

  my $obj = $self->_obj;

  if(!$obj) {
    warn("Bio::EnsEMBL::Container: potential memory leak, contained"
	 . " object is not defined during garbage collection.");
  } elsif($obj->can('deleteObj')) {
    $obj->deleteObj();
  }

  $self->{_obj} = undef;
}

1;
