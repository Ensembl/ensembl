
#
# BioPerl module for Bio::EnsEMBL::DBOLD::BaseAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBOLD::BaseAdaptor - Base Adaptor for DBOLD adaptors

=head1 SYNOPSIS

    # base adaptor provides
    
    # SQL prepare function
    $adaptor->prepare("sql statement");

    # get of root db object
    $adaptor->db();

    # delete memory cycles
    $adaptor->deleteObj();

    # constructor, ok for inheritence
    $adaptor = Bio::EnsEMBL::DBOLD::ClassWhichInheritsFromBaseAdaptor->new($dbobj)

=head1 DESCRIPTION

This is a true base class for Adaptors in the Ensembl DBOLD
system. Original idea from Arne

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBOLD::BaseAdaptor;
use vars qw(@ISA);
use strict;
use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI);

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

=head2 prepare

 Title   : prepare
 Usage   : $sth = $adaptor->prepare("select yadda from blabla")
 Function: provides a DBI statement handle from the adaptor. A convience
           function so you don't have to write $adaptor->db->prepare all the
           time
 Example :
 Returns : 
 Args    :


=cut

sub prepare{
   my ($self,$string) = @_;

   return $self->db->prepare($string);
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


