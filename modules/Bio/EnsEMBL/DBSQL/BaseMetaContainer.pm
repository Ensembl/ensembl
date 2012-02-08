=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::BaseMetaContainer - Encapsulates all generic access
to database meta information

=head1 SYNOPSIS

  my $meta_container = $db_adaptor->get_MetaContainer();

  my @mapping_info =
    @{ $meta_container->list_value_by_key('assembly.mapping') };

=head1 DESCRIPTION

  An object that encapsulates access to db meta data

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::BaseMetaContainer;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

# new() is inherited from Bio::EnsEMBL::DBSQL::BaseAdaptor

=head2 get_schema_version

  Arg [1]    : none
  Example    : $schema_ver = $meta_container->get_schema_version();
  Description: Retrieves the schema version from the database meta table
  Returntype : int
  Exceptions : none
  Caller     : ?
  Status     : Medium risk

=cut

sub get_schema_version {
  my $self = shift;

  my $arrRef = $self->list_value_by_key('schema_version');

  if (@$arrRef) {
    my ($ver) = ( $arrRef->[0] =~ /^\s*(\d+)\s*$/ );
    if ( !defined($ver) ) {    # old style format
      return 0;
    }
    return $ver;
  } else {
    warning(
      sprintf(
        "Please insert meta_key 'schema_version' "
          . "in meta table on core database '%s'\n",
        $self->dbc()->dbname() ) );
  }

  return 0;
}


=head2 list_value_by_key

  Arg [1]    : string $key
               the key to obtain values from the meta table with
  Example    : my @values = @{ $meta_container->list_value_by_key($key) };
  Description: gets a value for a key. Can be anything 
  Returntype : listref of strings 
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub list_value_by_key {
  my ( $self, $key ) = @_;

  $self->{'cache'} ||= {};

  if ( exists $self->{'cache'}->{$key} ) {
    return $self->{'cache'}->{$key};
  }

  my $sth;

  if ( !$self->_species_specific_key($key) ) {
    $sth =
      $self->prepare( "SELECT meta_value "
        . "FROM meta "
        . "WHERE meta_key = ? "
        . "AND species_id IS NULL "
        . "ORDER BY meta_id" );
  } else {
    $sth =
      $self->prepare( "SELECT meta_value "
        . "FROM meta "
        . "WHERE meta_key = ? "
        . "AND species_id = ? "
        . "ORDER BY meta_id" );
    $sth->bind_param( 2, $self->species_id(), SQL_INTEGER );
  }

  $sth->bind_param( 1, $key, SQL_VARCHAR );
  $sth->execute();

  my @result;
  while ( my $arrRef = $sth->fetchrow_arrayref() ) {
    push( @result, $arrRef->[0] );
  }

  $sth->finish();
  $self->{'cache'}->{$key} = \@result;

  return \@result;
} ## end sub list_value_by_key

=head2 single_value_by_key

  Arg [1]    : string $key
               the key to obtain values from the meta table with
  Arg [2]    : boolean $warn
               If true will cause the code to warn the non-existence of a value
  Example    : my $value = $mc->single_value_by_key($key);
  Description: Gets a value for a key. Can be anything
  Returntype : Scalar
  Exceptions : Raised if more than 1 meta item is returned

=cut

sub single_value_by_key {
  my ($self, $key, $warn) = @_;
  my $results = $self->list_value_by_key($key);
  if(defined $results) {
    my $count = scalar(@{$results});
    if($count == 1) {
      my ($value) = @{$results};
      return $value;
    }
    elsif($count == 0) {
      if($warn) {
        my $group = $self->db()->group();
        my $msg = sprintf(qq{Please insert meta_key '%s' in meta table at %s db\n}, $key, $group);
        warning($msg);
      }
    }
    else {
      my $values = join(q{,}, @{$results});
      throw sprintf(q{Found the values [%s] for the key '%s'}, $values, $key);
    }
  }
  return;
} ## end sub single_value_by_key

=head2 store_key_value

  Arg [1]    : string $key
               a key under which $value should be stored
  Arg [2]    : string $value
               the value to store in the meta table
  Example    : $meta_container->store_key_value($key, $value);
  Description: stores a value in the meta container, accessable by a key
  Returntype : none
  Exceptions : Thrown if the key/value already exists.
  Caller     : ?
  Status     : Stable

=cut

sub store_key_value {
  my ( $self, $key, $value ) = @_;

  if ( $self->key_value_exists( $key, $value ) ) {
    warn(   "Key-value pair '$key'-'$value' "
          . "already exists in the meta table; "
          . "not storing duplicate" );
    return;
  }

  my $sth;

  if ( !$self->_species_specific_key($key) ) {
    $sth = $self->prepare(
          'INSERT INTO meta (meta_key, meta_value, species_id) '
        . 'VALUES(?, ?, \N)' );
  } else {
    $sth = $self->prepare(
          'INSERT INTO meta (meta_key, meta_value, species_id) '
        . 'VALUES (?, ?, ?)' );
    $sth->bind_param( 3, $self->species_id(), SQL_INTEGER );
  }

  $sth->bind_param( 1, $key,   SQL_VARCHAR );
  $sth->bind_param( 2, $value, SQL_VARCHAR );
  $sth->execute();

  $self->{'cache'} ||= {};

  delete $self->{'cache'}->{$key};
} ## end sub store_key_value

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
  Status     : Stable

=cut

sub update_key_value {
  my ( $self, $key, $value ) = @_;

  my $sth;

  if ( !$self->_species_specific_key($key) ) {
    $sth =
      $self->prepare( 'UPDATE meta SET meta_value = ? '
        . 'WHERE meta_key = ?'
        . 'AND species_id IS NULL' );
  } else {
    $sth =
      $self->prepare( 'UPDATE meta '
        . 'SET meta_value = ? '
        . 'WHERE meta_key = ? '
        . 'AND species_id = ?' );
    $sth->bind_param( 3, $self->species_id(), SQL_INTEGER );
  }

  $sth->bind_param( 1, $value, SQL_VARCHAR );
  $sth->bind_param( 2, $key,   SQL_VARCHAR );
  $sth->execute();

} ## end sub update_key_value


=head2 delete_key

  Arg [1]    : string $key
               The key which should be removed from the database.
  Example    : $meta_container->delete_key('sequence.compression');
  Description: Removes all rows from the meta table which have a meta_key
               equal to $key.
  Returntype : none
  Exceptions : none
  Caller     : dna_compress script, general
  Status     : Stable

=cut

sub delete_key {
  my ( $self, $key ) = @_;

  my $sth;

  if ( !$self->_species_specific_key($key) ) {
    $sth =
      $self->prepare( 'DELETE FROM meta '
        . 'WHERE meta_key = ?'
        . 'AND species_id IS NULL' );
  } else {
    $sth =
      $self->prepare( 'DELETE FROM meta '
        . 'WHERE meta_key = ? '
        . 'AND species_id = ?' );
    $sth->bind_param( 2, $self->species_id(), SQL_INTEGER );
  }

  $sth->bind_param( 1, $key, SQL_VARCHAR );
  $sth->execute();

  delete $self->{'cache'}->{$key};
}

=head2 delete_key_value

  Arg [1]    : string $key
               The key which should be removed from the database.
  Arg [2]    : string $value
               The value to be removed.
  Example    : $meta_container->delete_key('patch', 'patch_39_40_b.sql|xref_unique_constraint');
  Description: Removes all rows from the meta table which have a meta_key
               equal to $key, AND a meta_value equal to $value.
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub delete_key_value {
  my ( $self, $key, $value ) = @_;

  my $sth;

  if ( !$self->_species_specific_key($key) ) {
    $sth =
      $self->prepare( 'DELETE FROM meta '
        . 'WHERE meta_key = ? '
        . 'AND meta_value = ?'
        . 'AND species_id IS NULL' );
  } else {
    $sth =
      $self->prepare( 'DELETE FROM meta '
        . 'WHERE meta_key = ? '
        . 'AND meta_value = ? '
        . 'AND species_id = ?' );
    $sth->bind_param( 3, $self->species_id(), SQL_INTEGER );
  }

  $sth->bind_param( 1, $key,   SQL_VARCHAR );
  $sth->bind_param( 2, $value, SQL_VARCHAR );
  $sth->execute();

  delete $self->{'cache'}->{$key};
} ## end sub delete_key_value

=head2 key_value_exists

  Arg [1]    : string $key
               the key to check
  Arg [2]    : string $value
               the value to check
  Example    : if ($meta_container->key_value_exists($key, $value)) ...
  Description: Return true (1) if a particular key/value pair exists,
               false (0) otherwise.
  Returntype : boolean
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub key_value_exists {
  my ( $self, $key, $value ) = @_;

  my $sth;

  if ( !$self->_species_specific_key($key) ) {
    $sth =
      $self->prepare( 'SELECT meta_value '
        . 'FROM meta '
        . 'WHERE meta_key = ? '
        . 'AND meta_value = ?'
        . 'AND species_id IS NULL' );
  } else {
    $sth =
      $self->prepare( 'SELECT meta_value '
        . 'FROM meta '
        . 'WHERE meta_key = ? '
        . 'AND meta_value = ? '
        . 'AND species_id = ?' );
    $sth->bind_param( 3, $self->species_id(), SQL_INTEGER );
  }

  $sth->bind_param( 1, $key,   SQL_VARCHAR );
  $sth->bind_param( 2, $value, SQL_VARCHAR );
  $sth->execute();

  while ( my $arrRef = $sth->fetchrow_arrayref() ) {
    if ( $arrRef->[0] eq $value ) {
      $sth->finish();
      return 1;
    }
  }

  return 0;
} ## end sub key_value_exists

# This utility method determines whether the key is a species-specific
# meta key or not.  If the key is either 'patch' or 'schema_version',
# then it is not species-specific.

# FIXME variation team messed up in release 65 and added the ploidy
# entry without species_id - this will be corrected for release 66,
# for now, I've added it to the list of allowed non-species specific

sub _species_specific_key {
  my ( $self, $key ) = @_;

  return (    $key ne 'patch'
           && $key ne 'schema_version'
           && $key ne 'schema_type'
           && $key ne 'ploidy');
}

1;
