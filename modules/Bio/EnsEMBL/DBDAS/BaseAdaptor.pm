
#
# BioPerl module for Bio::EnsEMBL::DBDAS::BaseAdaptor
#
# Cared for by Tony Cox <avc@sanger.ac.uk>
#
# Copyright Tony Cox
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBDAS::BaseAdaptor - Base Adaptor for DBSQL adaptors

=head1 SYNOPSIS

    # base adaptor provides
    
    # get of root db object
    $adaptor->db();

    # delete memory cycles
    $adaptor->deleteObj();

    # constructor, ok for inheritence
    $adaptor = Bio::EnsEMBL::DBDAS::ClassWhichInheritsFromBaseAdaptor->new($dbobj)

=head1 DESCRIPTION

This is a true base class for Adaptors in the Ensembl DBDAS
system. Original idea from Arne

=head1 CONTACT

Contact Tony on x4763

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal 
methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBDAS::BaseAdaptor;
use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Root);

sub new {
    my ($class,$dbobj) = @_;

    my $self = {};
    bless $self,$class;

    if( !defined $dbobj || !ref $dbobj ) {
	$self->throw("Don't have a db [$dbobj] for new adaptor");
    }

    $self->db($dbobj);

    return $self;
}

=head2 db

 Title   : db
 Usage   : $obj->db($newval)
 Function: 
 Returns : value of db
 Args    : newvalue (optional)


=cut

sub db{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'db'} = $value;
    }
    return $obj->{'db'};

}


=head2 deleteObj

 Title   : deleteObj
 Usage   : $obj->deleteObj
 Function: removes memory cycles. Probably triggered by Root deleteObj
 Returns : 
 Args    : none


=cut

sub deleteObj {
  my $self = shift;
  my @dummy = values %{$self};
  foreach my $key ( keys %$self ) {
    delete $self->{$key};
  }
  foreach my $obj ( @dummy ) {
    eval {
      $obj->deleteObj;
    }
  }
}


