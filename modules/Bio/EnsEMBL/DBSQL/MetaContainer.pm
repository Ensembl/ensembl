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

Bio::EnsEMBL::DBSQL::MetaContainer - Encapsulates all access to core
database meta information

=head1 SYNOPSIS

  my $meta_container =
    $registry->get_adaptor( 'Human', 'Core', 'MetaContainer' );

  my @mapping_info =
    @{ $meta_container->list_value_by_key('assembly.mapping') };
  
  my $scientific_name = $meta_container->get_scientific_name();

=head1 DESCRIPTION

  An object that encapsulates specific access to core db meta data

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::MetaContainer;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw/deprecate/;
use Bio::Species;


use base qw/Bio::EnsEMBL::DBSQL::BaseMetaContainer/;

# add well known meta info get-functions below

=head2 get_production_name

  Args          : none
  Example       : $species = $meta_container->get_production_name();
  Description   : Obtains the name of the species in a form usable as, for
                  example, a table name, file name etc.
  Returntype    : string
  Exceptions    : none
  Status        : Stable

=cut

sub get_production_name {
  my ($self) = @_;
  return $self->single_value_by_key('species.production_name');
}

=head2 get_short_name

  Args          : none
  Example       : $species = $meta_container->get_short_name();
  Description   : Obtains the name of the species in a form usable as, for
                  example, a short label in a GUI.
  Returntype    : string
  Exceptions    : none
  Status        : Stable

=cut

sub get_short_name {
  my ($self) = @_;
  return $self->single_value_by_key('species.short_name');
}

=head2 get_common_name

  Args          : none
  Example       : $species = $meta_container->get_common_name();
  Description   : Obtains the common name of the species.
  Returntype    : string
  Exceptions    : none
  Status        : Stable

=cut

sub get_common_name {
  my ($self) = @_;
  return $self->single_value_by_key('species.common_name');
}

=head2 get_scientific_name

  Args          : none
  Example       : $species = $meta_container->get_scientific_name();
  Description   : Obtains the full scientific name of the species.
  Returntype    : string
  Exceptions    : none
  Status        : Stable

=cut
sub get_scientific_name {
  my ($self) = @_;
  return $self->single_value_by_key('species.scientific_name');
}

=head2 get_division

  Args          : none
  Example       : $div = $meta_container->get_division();
  Description   : Obtains the Ensembl Genomes division to which the species belongs.
  Returntype    : string
  Exceptions    : none
  Status        : Stable

=cut
sub get_division {
  my ($self) = @_;
  return $self->single_value_by_key('species.division');
}

=head2 get_Species

  Arg [1]    : none
  Example    : $species = $meta_container->get_Species();
  Description: Obtains the species from this databases meta table. Call is
               deprecated; please use other subroutines in this package
  Returntype : Bio::Species
  Exceptions : none
  Caller     : ?
  Status     : Deprecated

=cut

sub get_Species {
  my ($self) = @_;

  deprecate('Call is deprecated. Use $self->get_common_name() / $self->get_classification() / $self->get_scientific_name() instead');

  my $common_name = $self->get_common_name();
  my $classification =
    $self->list_value_by_key('species.classification');
  if ( !@$classification ) {
    return undef;
  }
  
  #Re-create the old classification data structure by adding the scientific
  #name back onto the classification but with species before genus e.g.
  # sapiens -> Homo -> Hominade
  my $scientific_name = $self->get_scientific_name();
  my ($genus, @sp) = split(/ /, $scientific_name);
  unshift(@{$classification}, join(q{ }, @sp), $genus);
  
  my $species = Bio::Species->new();
  $species->common_name($common_name);
  
  $species->classification($classification, 1); #always force it

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
  my ($self) = @_;
  return $self->single_value_by_key('species.taxonomy_id', 1);
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
  return $self->single_value_by_key('assembly.maxcontig');
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
  my ($self) = @_;
  return $self->single_value_by_key('genebuild.start_date', 1);
}

=head2 get_genebuild

  Example    : $classification = $meta_container->get_classification();
  Description: Retrieves the classification held in the backing database minus
               any species specific levels. This means that the first element
               in the array will be subfamily/family level ascending to
               superkingdom
  Returntype : ArrayRef[String]
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub get_classification {
  my ($self) = @_;
  my $classification = $self->list_value_by_key('species.classification');
  my $copy = [@{$classification}];
  splice(@{$copy}, 0, 1); # remove the Homo sapiens
  return $copy;
}


1;

