=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

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

=head2 get_display_name

  Args          : none
  Example       : $species = $meta_container->get_display_name();
  Description   : Obtains the name of the species in a form usable as, for
                  example, a short label in a GUI.
  Returntype    : string
  Exceptions    : none
  Status        : Stable

=cut

sub get_display_name {
  my ($self) = @_;
  return $self->single_value_by_key('species.display_name');
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

