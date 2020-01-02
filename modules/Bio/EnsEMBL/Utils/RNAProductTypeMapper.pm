
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Utils::RNAProductTypeMapper - Utility class for mapping
between RNA-product types used in the Ensembl database and respective
Perl API classes.

=head1 DESCRIPTION

The purpose of this class is to hide from the users the mappings
between classes representing various mature RNA products
(e.g. generic, miRNA, circRNA etc.) and their respective identifiers
in the Ensembl database. In principle there should be no need for
users to call the mapper directly; it is invoked internally whenever
such mappings are needed (e.g. in RNAProduct constructor or in
RNAProductAdaptor).

Note that the type_code<->class_name mappings are hardcoded here,
specifically in the constructor, instead of being fetched from the
database. This is so that it is not possible for someone to trigger
construction of arbitrary objects by having modified class names
stored in the database.

=head1 SYNOPSIS

  my $rpt_mapper = Bio::EnsEMBL::Utils::RNAProductTypeMapper->new();

  my $class_name = $rpt_mapper->type_code_to_class( 'generic' );
  my $type_code = $rpt_mapper->class_to_type_code( 'Bio::EnsEMBL::MicroRNA' );

=cut

package Bio::EnsEMBL::Utils::RNAProductTypeMapper;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw( throw warning );



my $type_mapper;


=head2 mapper

  Example    : my $mapper = Bio::EnsEMBL::Utils::RNAProductTypeMapper->mapper();
  Description: Retrieves an instance of RNAProductTypeMapper from its module,
               reusing an existing one should it exist.
  Returntype : Bio::EnsEMBL::Utils::RNAProductTypeMapper
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub mapper {
  if ( !defined($type_mapper) ) {
    $type_mapper = Bio::EnsEMBL::Utils::RNAProductTypeMapper->new();
  }
  return $type_mapper;
}


=head2 new

  Example    : my $mapper = Bio::EnsEMBL::Utils::RNAProductTypeMapper->new();
  Description: Constructor.  Creates a new RNAProductTypeMapper object.
               Not particularly useful for end users because all
               RNAProductTypeMapper objects are identical; use mapper()
               instead.
  Returntype : Bio::EnsEMBL::Utils::RNAProductTypeMapper
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub new {
  my ( $caller ) = @_;
  my $class = ref( $caller ) || $caller;

  # Declare this here rather than in the package scope so that the map
  # cannot be modified (we do not presently use Readonly in the Core
  # API code), accidentally or otherwise.
  my $class_attribute_cache_map = {
    'Bio::EnsEMBL::RNAProduct' => { },
    'Bio::EnsEMBL::MicroRNA'   => {
      'arm' => 'mirna_arm',
    },
  };
  my $type_to_class_map = {
    'generic' => 'Bio::EnsEMBL::RNAProduct',
    'miRNA'   => 'Bio::EnsEMBL::MicroRNA',
  };

  my $self = bless {
    'class_attribute_cache_map' => $class_attribute_cache_map,
    'type_to_class_map'         => $type_to_class_map,
    'class_to_type_map'         => undef,
  }, $class;

  return $self;
}


=head2 class_attribute_cache_map

  Arg [1]    : string $class_name - fully qualified RNA-product class name
  Example    : my $attr_cache_map
                 = $mapper->class_attribute_cache_map( 'Bio::EnsEMBL::MicroRNA' );
  Description: For the given name of a class representing a mature RNA
               product, returns the map indicating which local members variables
               should be synchronised with which Attributes.
  Returntype : hashref
  Exceptions : throw if the class does not represent known RNA-product type
  Caller     : internal
  Status     : Stable

=cut

sub class_attribute_cache_map {
  my ( $self, $class_name ) = @_;

  my %map = %{ $self->{'class_attribute_cache_map'} };
  if ( !exists $map{$class_name} ) {
    throw( "Unknown RNA-product class name " . $class_name );
  }

  return $map{$class_name};
}


=head2 class_to_type_code

  Arg [1]    : string $class_name - fully qualified RNA-product class name
  Example    : my $type_code
                 = $mapper->class_to_type_code( 'Bio::EnsEMBL::MicroRNA' );
  Description: For the given name of a class representing a mature RNA
               product, returns the type code used to represent it in the
               Ensembl database.
  Returntype : string
  Exceptions : throw if the class does not represent known RNA-product type
  Caller     : internal
  Status     : Stable

=cut

sub class_to_type_code {
  my ( $self, $class_name ) = @_;

  if ( !defined( $self->{'class_to_type_map'} ) ) {
    $self->_generate_reverse_map();
  }

  my %map = %{ $self->{'class_to_type_map'} };
  if ( !exists $map{$class_name} ) {
    throw( "Unknown RNA-product class name " . $class_name );
  }

  return $map{$class_name};
}


=head2 type_code_to_class

  Arg [1]    : string $type_code - type code of RNA product
  Example    : my $class_name = $mapper->class_to_type_code( 1 );
  Description: For the type code representing a mature RNA product in the
               Ensembl database, return its API class name
  Returntype : string
  Exceptions : throw if the code does not represent known RNA-product type
  Caller     : internal
  Status     : Stable

=cut

sub type_code_to_class {
  my ( $self, $type_code ) = @_;

  my %map = %{ $self->{'type_to_class_map'} };
  if ( !exists $map{$type_code} ) {
    throw( "Unknown RNA-product type ID " . $type_code );
  }

  return $map{$type_code};
}


=head2 _generate_reverse_map

  Description: PRIVATE generates class_name->type_code map from the
               type_code->class_name one.
  Returntype : none
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub _generate_reverse_map {
  my ($self) = @_;

  # Safe to use reverse because both keys and values are unique
  my %reversed_map = reverse %{$self->{'type_to_class_map'}};
  $self->{'class_to_type_map'} = \%reversed_map;

  return;
}


1;
