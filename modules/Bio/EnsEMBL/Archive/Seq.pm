#
# EnsEMBL module for Bio::EnsEMBL::Archive::Seq
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Archive::Seq - Simple archive seq object

=head1 SYNOPSIS

   my $seq = Bio::EnsEMBL::Archive::Seq->new(
					     -name => 'ENSE0006734',
					     -type => 'exon',
					     );

   my $versioned_seq = Bio::EnsEMBL::Archive::VersionedSeq->new(
				     -seq => $seq,
				     -version => 1,
				     -sequence => 'ATCGTAGAT',
				     -start_clone => 'AC000234',
                                     -seq_start => 23423,
				     -end_clone => 'AC123132',
				     -seq_end => 1243,
				     -release_number => 100,
				     );

=head1 DESCRIPTION

This object is a very bare representation of the seq table within the archive.
This object is only in used in conjunction with the versioned_seq object, since
an entry in the seq table can never exist without at least an entry in the versioned_seq table.

Multiple versioned objects will have the same underlying seq object.

=head1 CONTACT

Elia Stupka - elia@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Archive::Seq;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root
# (inherits for methods like throw and rearrange)

use Bio::EnsEMBL::Root;


@ISA = qw(Bio::EnsEMBL::Root);

sub new {
  my($class, @args) = @_;
  
  my $self = {};
  bless $self, $class;

   my ($dbid,$name,$type,$created,$adaptor) = $self->_rearrange([qw(
					  DBID
					  NAME
					  TYPE
					  CREATED
			                  ADAPTOR
					  )],@args);

  $name || $self->throw("An Archive Seq object must have a name");
  $type || $self->throw("An Archive seq object must have a type");
  $created || $self->throw("An Archive seq object must have a created date"); 
 
  if ($adaptor) {
      $dbid || $self->throw("Creating Archive Seq with adaptor but without a db_ID");
      $self->adaptor($adaptor);
      $self->db_ID($dbid);
  }
  $self->name($name);
  $self->type($type);
  $self->created($created);

  return $self;
}

=head2 name

 Title   : name
 Usage   : $obj->name($newval)
 Function: Getset for name value
 Returns : value of name
 Args    : newvalue (optional)


=cut

sub name{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'name'} = $value;
    }
    return $obj->{'name'};

}

=head2 type

 Title   : type
 Usage   : $obj->type($newval)
 Function: Getset for type value
 Returns : value of type
 Args    : newvalue (optional)


=cut

sub type{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'type'} = $value;
    }
    return $obj->{'type'};

}

=head2 created

 Title   : created
 Usage   : $obj->created($newval)
 Function: Getset for created value
 Returns : value of created
 Args    : newvalue (optional)


=cut

sub created{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'created'} = $value;
    }
    return $obj->{'created'};

}


=head2 adaptor

 Title   : adaptor
 Usage   : $obj->adaptor($newval)
 Function: Getset for adaptor value
 Returns : value of adaptor
 Args    : newvalue (optional)


=cut

sub adaptor{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'adaptor'} = $value;
    }
    return $obj->{'adaptor'};

}

=head2 db_ID

 Title   : db_ID
 Usage   : $obj->db_ID($newval)
 Function: Getset for db_ID value
 Returns : value of db_ID
 Args    : newvalue (optional)


=cut

sub db_ID{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'db_ID'} = $value;
    }
    return $obj->{'db_ID'};

}
