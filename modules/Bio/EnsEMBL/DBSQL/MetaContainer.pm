#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::MetaContainer
#
# Cared for by Arne Stabenau
#
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

  Bio::EnsEMBL::DBSQL::MetaContainer - 
  Encapsulates all access to database meta information

=head1 SYNOPSIS

  my $meta_container = $db_adaptor->get_MetaContainer();
  my $assembly_type = $meta_container->get_default_assembly();  

=head1 DESCRIPTION

  An object that encapsulates access to db meta data

=head1 CONTACT

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::DBSQL::MetaContainer;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::Species;


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



# add well known meta info get-functions below

=head2 get_Species

  Arg [1]    : none
  Example    : $species = $meta_container->get_Species();
  Description: Obtains the species from this databases meta table
  Returntype : Bio::Species
  Exceptions : none
  Caller     : ?

=cut

sub get_Species {
  my $self = shift;

  my $arrRef = $self->list_value_by_key( 'species.common_name' );
  my $common_name;
  if( @$arrRef ) {
    $common_name = $arrRef->[0];
  } else {
    return undef;
  }
  
  my $classification = $self->list_value_by_key( 'species.classification' );
  if( ! @$classification ) {
    return undef;
  }

  my $species = new Bio::Species;
  $species->common_name( $common_name );
  $species->classification( @$classification );

  return $species;
}


=head2 get_taxonomy_id

  Arg [1]    : none
  Example    : $tax_id = $meta_container->get_taxonomy_id();
  Description: Retrieves the taxonomy id from the database meta table
  Returntype : string
  Exceptions : none
  Caller     : ?

=cut

sub get_taxonomy_id {
  my $self = shift;

  my $arrRef = $self->list_value_by_key( 'species.taxonomy_id' );
  
  if( @$arrRef ) {
    return $arrRef->[0];
  } else {
    $self->warn("Please insert meta_key 'species.taxonomy_id' " .
		"in meta table at core db.\n");
  }
}



=head2 get_top_coord_system

  Arg [1]    : none
  Example    : $top_level = $meta_container->get_top_coord_system;
  Description: Retrieves the top-level coordinate system for this database.
               For a genome that is fully assembled this will be 'chromosome',
               but for a shattered genome this might be 'supercontig' or
               'clone'.
  Returntype : string
  Exceptions : none
  Caller     : web?, SliceAdaptor

=cut


sub get_top_coord_system {
  my $self = shift;

  my ($top) = @{$self->get_all_coord_systems};

  return $top;
}



=head2 is_valid_coord_system

  Arg [1]    : none
  Example    : do_something() if ($meta_container->is_valid_coord_system($cs));
  Description: Returns true if a provided coordinate system is valid.  I.e.
               that coordinate system is present in the meta table.
  Returntype : string
  Exceptions : none
  Caller     : web?, SliceAdaptor

=cut

sub is_valid_coord_system {
  my $self = shift;
  my $coord_system = lc(shift);

  if(!exists($self->{'valid_coord_systems'})) {
    my %vcs = ();
    foreach my $value (@{$self->get_all_coord_systems()}) {
      $vcs{lc($value)} = 1;
    }
    $self->{'valid_coord_systems'} = \%vcs;
  }

  return $self->{'valid_coord_systems'}->{$coord_system};
}


=head2 get_all_coord_systems

  Arg [1]    : none
  Example    : my @coord_systems = @{$meta_container->get_all_coord_systems};
  Description: Returns a listref of the valid coordinate systems for this
               database.
  Returntype : string
  Exceptions : none
  Caller     : web?

=cut

sub get_all_coord_systems {
  my $self = shift;
  return $self->list_value_by_key('coord_systems');
}





=head2 get_default_assembly

  Arg [1]    : none
  Example    : $assembly = $meta_container->get_default_assembly();
  Description: Retrieves the default assembly for this database from the 
               meta container
  Returntype : string
  Exceptions : none
  Caller     : ?

=cut

#
# TBD This method should be removed/deprecated
#
sub get_default_assembly {
  my $self = shift;

  my $arrRef = $self->list_value_by_key('assembly.default' );
 
  if( @$arrRef ) {
    return $arrRef->[0];
  } else {
    return undef;
  }
}

#
# TBD This method should be removed/deprecated
#
sub get_max_assembly_contig {
  my $self = shift;

  my $value_list = $self->list_value_by_key( "assembly.maxcontig" );
  if( @$value_list ) {
    return $value_list->[0];
  } else {
    return undef;
  }
}



1;

