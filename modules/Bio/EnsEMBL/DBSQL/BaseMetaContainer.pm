#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::BaseMetaContainer
#
# Cared for by Arne Stabenau
#
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

  Bio::EnsEMBL::DBSQL::BaseMetaContainer - 
  Encapsulates all generic access to database meta information

=head1 SYNOPSIS

  my $meta_container = $db_adaptor->get_MetaContainer();

  my @mapping_info = @{$meta_container->list_value_by_key('assembly.mapping')};

=head1 DESCRIPTION

  An object that encapsulates access to db meta data

=head1 CONTACT

  Post questions to the EnsEMBL development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::BaseMetaContainer;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(deprecate);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

# new() is inherited from Bio::EnsEMBL::DBSQL::BaseAdaptor


=head2 list_value_by_key

  Arg [1]    : string $key
               the key to obtain values from the meta table with
  Example    : my @values = $meta_container->list_value_by_key($key);
  Description: gets a value for a key. Can be anything 
  Returntype : listref of strings 
  Exceptions : none
  Caller     : ?

=cut

sub list_value_by_key {
  my ($self,$key) = @_;
  my @result;

  $self->{'cache'} ||= {};
  if( exists $self->{'cache'}->{$key} ) {
    return $self->{'cache'}->{$key};
  }

  my $sth = $self->prepare( "SELECT meta_value 
                             FROM meta 
                             WHERE meta_key = ? ORDER BY meta_id" );
  $sth->execute( $key );
  while( my $arrRef = $sth->fetchrow_arrayref() ) {
    push( @result, $arrRef->[0] );
  }
  $self->{'cache'}->{$key} = \@result;

  return \@result;
}


=head2 store_key_value

  Arg [1]    : string $key
               a key under which $value should be stored
  Arg [2]    : string $value
               the value to store in the meta table
  Example    : $meta_container->store_key_value($key, $value);
  Description: stores a value in the meta container, accessable by a key
  Returntype : none
  Exceptions : none
  Caller     : ?

=cut

sub store_key_value {
  my ( $self, $key, $value ) = @_;

  my $sth = $self->prepare( "INSERT INTO meta( meta_key, meta_value) 
                             VALUES( ?, ? )" );

  my $res = $sth->execute( $key, $value );

  $self->{'cache'} ||= {};

  delete $self->{'cache'}->{$key};

  return;
}

=head2 update_key_value

  Arg [1]    : string $key
               a key under which $value should be updated
  Arg [2]    : string $value
               the value to update in the meta table
  Example    : $meta_container->update_key_value($key, $value);
  Description: update a value in the meta container, accessable by a key
  Returntype : none
  Exceptions : none
  Caller     : ?

=cut

sub update_key_value {
  my ( $self, $key, $value ) = @_;

  my $sth = $self->prepare( "UPDATE meta SET meta_value = ? WHERE meta_key = ?" );

  my $res = $sth->execute( $value, $key );
  return;
}


=head2 delete_key

  Arg [1]    : string $key
               The key which should be removed from the database.
  Example    : $meta_container->delete_key('sequence.compression');
  Description: Removes all rows from the meta table which have a meta_key
               equal to $key.
  Returntype : none
  Exceptions : none
  Caller     : dna_compress script, general

=cut

sub delete_key {
  my ($self, $key) = @_;

  my $sth = $self->prepare("DELETE FROM meta WHERE meta_key = ?");
  $sth->execute($key);
  $sth->finish();

  delete $self->{'cache'}->{$key};

  return;
}

1;
