=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::MetaContainer - Encapsulates all access to core
database meta information

=head1 SYNOPSIS

  my $meta_container =
    $registry->get_adaptor( 'Human', 'Core', 'MetaContainer' );

  my $species = $meta_container->get_Species();

  my @mapping_info =
    @{ $meta_container->list_value_by_key('assembly.mapping') };

=head1 DESCRIPTION

  An object that encapsulates specific access to core db meta data

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::MetaContainer;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseMetaContainer;
use Bio::EnsEMBL::Utils::Exception;

use Bio::Species;


@ISA = qw(Bio::EnsEMBL::DBSQL::BaseMetaContainer);

# add well known meta info get-functions below



=head2 get_Species

  Arg [1]    : none
  Example    : $species = $meta_container->get_Species();
  Description: Obtains the species from this databases meta table
  Returntype : Bio::Species
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub get_Species {
  my $self = shift;

  my $arrRef = $self->list_value_by_key( 'species.common_name' );
  my $common_name;
  if( @$arrRef ) {
    $common_name = $arrRef->[0];
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
  Status     : Stable

=cut

sub get_taxonomy_id {
  my $self = shift;

  my $arrRef = $self->list_value_by_key( 'species.taxonomy_id' );
  
  if( @$arrRef ) {
    return $arrRef->[0];
  } else {
    warning("Please insert meta_key 'species.taxonomy_id' " .
	    "in meta table at core db.\n");
  }
}



=head2 get_default_assembly

  Description: DEPRECATED. Use the version of the coordinate system you are
             interested in instead.

  Example:     #use this instead
               my ($highest_cs) = @{$db->get_CoordSystemAdaptor->fetch_all()};
               my $assembly = $highest_cs->version();

=cut

sub get_default_assembly {
  my $self = shift;

  deprecate("Use version of coordinate system you are interested in instead.\n".
            "Example:\n".
            '  ($cs) = @{$coord_system_adaptor->fetch_all()};'."\n" .
            '  $assembly = $cs->version();');

  my ($cs) = @{$self->db->get_CoordSystemAdaptor->fetch_all()};

  return $cs->version();
}


#
# TBD This method should be removed/deprecated
#
sub get_max_assembly_contig {
  my $self = shift;

  deprecate('This method should either be fixed or removed');

  my $value_list = $self->list_value_by_key( "assembly.maxcontig" );
  if( @$value_list ) {
    return $value_list->[0];
  } else {
    return undef;
  }
}

=head2 get_genebuild

  Arg [1]    : none
  Example    : $tax_id = $meta_container->get_genebuild();
  Description: Retrieves the genebuild from the database meta table
  Returntype : string
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub get_genebuild {
  my $self = shift;

  my $arrRef = $self->list_value_by_key( 'genebuild.start_date' );

  if( @$arrRef ) {
    return $arrRef->[0];
  } else {
    warning("Please insert meta_key 'genebuild.start_date' " .
            "in meta table at core db.\n");
  }
}


1;

