#
# EnsEMBL module for Bio::EnsEMBL::Storable
#

=head1 NAME

Bio::EnsEMBL::Storable

=head1 SYNOPSIS

  my $dbID = $storable_object->dbID();
  my $adaptor = $storable_object->adaptor();
  if($storable_object->is_stored($db_adaptor))) {
    ...
  }
=head1 DESCRIPTION

This is a simple object which contains a few coordinate system attributes:
name, internal identifier, version

=head1 AUTHOR - Graham McVicker

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=head1 APPENDIX

This is a storable base class.  All objects which are storable in the database
should inherit from this class.  It provides two getter/setters: dbID()
adaptor().  And a is_stored() method that can be used to determine if an
object is already stored in a database.

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Storable;

use Bio::EnsEMBL::Root;
use vars qw(@ISA);

#
# will eventually remove unneeded inheritance to Root
#
@ISA = qw(Bio::EnsEMBL::Root);


use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($adaptor, $dbID) = rearrange(['ADAPTOR', 'dbID'],@_);

  if($adaptor) {
    if(!ref($adaptor) || !$adaptor->isa('Bio::EnsEMBL::DBSQL::BaseAdaptor')) {
      throw('-ADAPTOR argument must be a Bio::EnsEMBL::DBSQL::BaseAdaptor');
    }
  }

  return bless({'dbID' => $dbID, 'adaptor' => $adaptor}, $class);
}


=head2 dbID

  Arg [1]    : int $dbID
  Example    : none
  Description: getter/setter for the database internal id
  Returntype : int
  Exceptions : none
  Caller     : general, set from adaptor on store

=cut

sub dbID {
  my $self = shift;
  $self->{'dbID'} = shift if(@_);
  return $self->{'dbID'};
}



=head2 adaptor

  Arg [1]    : Bio::EnsEMBL::DBSQL::BaseAdaptor $adaptor
  Example    : none
  Description: get/set for this objects Adaptor
  Returntype : Bio::EnsEMBL::DBSQL::ChromsomeAdaptor
  Exceptions : none
  Caller     : general, set from adaptor on store

=cut

sub adaptor {
  my $self = shift;

  if(@_) {
    my $ad = shift;
    if($ad && (!ref($ad) || !$ad->isa('Bio::EnsEMBL::DBSQL::BaseAdaptor'))) {
      throw('Adaptor argument must be a Bio::EnsEMBL::DBSQL::BaseAdaptor');
    }
    $self->{'adaptor'} = $ad;
  }

  return $self->{'adaptor'}
}


=head2 is_stored

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBConnection
  Example    : do_something if($object->is_stored($db));
  Description: Returns true if this object is stored in the provided database.
               This works under the assumption that if the adaptor and dbID are
               set and the database of the adaptor shares the port, dbname and
               hostname with the provided database, this object is stored in
               that database.
  Returntype : 1 or 0
  Exceptions : throw if dbID is set but adaptor is not
               throw if adaptor is set but dbID is not
               throw if incorrect argument is passed
  Caller     : store methods

=cut

sub is_stored {
  my $self = shift;
  my $db = shift;

  if(!$db || !ref($db) || !$db->isa('Bio::EnsEMBL::DBSQL::DBConnection')) {
    throw('db argument must be a Bio::EnsEMBL::DBSQL::DBConnection');
  }

  my $adaptor = $self->{'adaptor'};
  my $dbID = $self->{'dbID'};

  if($dbID && !$adaptor) {
    throw("Storable object has a dbID but not an adaptor.\n" .
          'Storable objects must have neither OR both.');
  }

  if($adaptor && !$dbID) {
    throw("Storable object has an adaptor but not a dbID.\n".
          "Storable objects must have neither OR both.");
  }

  return 0 if (!$adaptor && !dbID);

  my $cur_db = $adaptor->db();

  #
  # Databases are the same if they share the same port, host and username
  #
  if($db->port   == $cur_db->port &&
     $db->host   eq $cur_db->host &&
     $db->dbname eq $cur_db->dbname) {
    return 1;
  }

  return 0;
}



1;
