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

  Bio::EnsEMBL::DBSQL::BiotypeAdaptor - An adaptor which performs database
  interaction relating to the storage and retrieval of Biotypes

=head1 SYNOPSIS

  my $biotype = $db_adaptor->fetch_by_name_object_type('protein_coding', 'gene');

=head1 DESCRIPTION

    This adaptor provides a means to retrieve and store information related
    to Biotypes.  Primarily this involves the retrieval or storage of
    Bio::EnsEMBL::Biotype objects from a database.

    See Bio::EnsEMBL::Biotype for details of the Biotype class.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::BiotypeAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);
use Bio::EnsEMBL::Biotype;

use strict;
use warnings;

use parent qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 _tables

  Arg [1]    : none
  Description: PROTECTED implementation of superclass abstract method.
               Returns the names, aliases of the tables to use for queries.
  Returntype : list of arrays of strings
  Exceptions : none

=cut

sub _tables {
  my $self = shift;

  return (['biotype', 'b']);
}

=head2 _columns

  Arg [1]    : none
  Example    : none
  Description: PROTECTED implementation of superclass abstract method.
               Returns a list of columns to use for queries.
  Returntype : list of strings
  Exceptions : none

=cut

sub _columns {
  my $self = shift;

  return ('b.biotype_id', 'b.name', 'b.object_type', 'b.db_type', 'b.attrib_type_id', 'b.description', 'b.biotype_group', 'b.so_acc', 'b.so_term');
}

=head2 _objs_from_sth

  Arg [1]    : StatementHandle $sth
  Example    : none
  Description: PROTECTED implementation of abstract superclass method.
               responsible for the creation of ProteinFeatures
  Returntype : arrayref of Bio::EnsEMBL::Biotype objects
  Exceptions : none

=cut

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my ($dbID, $name, $object_type, $db_type, $attrib_type_id, $description, $biotype_group, $so_acc, $so_term);

  $sth->bind_columns(\$dbID, \$name, \$object_type, \$db_type, \$attrib_type_id, \$description, \$biotype_group, \$so_acc, \$so_term);

  my @biotypes;

  while($sth->fetch()) {
    push( @biotypes,
      my $feat = Bio::EnsEMBL::Biotype->new_fast( {
         'dbID'           => $dbID,
         'name'           => $name,
         'object_type'    => $object_type,
         'db_type'        => $db_type,
         'attrib_type_id' => $attrib_type_id,
         'description'    => $description,
         'biotype_group'  => $biotype_group,
         'so_acc'         => $so_acc,
         'so_term'        => $so_term,
      } )
    );
  }

  return \@biotypes;
}


=head2 fetch_by_name_object_type

  Arg [1]    : String $name
               The name of the biotype to retrieve
  Arg [2]    : String $object_type
               The object type of the biotype to retrieve (gene or transcript)
  Example    : $biotype = $biotype_adaptor->fetch_by_name_object_type('mRNA', 'gene');
  Description: Retrieves a biotype object from the database via its combined key (name, object_type).
               If the Biotype requested does not exist in the database, a new Biotype object is
               created with the provided name and object_type to be returned.
  Returntype : Bio::EnsEMBL::Biotype
  Exceptions : none

=cut

sub fetch_by_name_object_type {
  my ($self, $name, $object_type) = @_;

  my $biotype;

  # biotype table implemented on schema_version 93
  if ($self->schema_version > 92) {
    my $constraint = "b.name = ? AND b.object_type = ?";
    $self->bind_param_generic_fetch($name, SQL_VARCHAR);
    $self->bind_param_generic_fetch($object_type, SQL_VARCHAR);
    $biotype = shift @{$self->generic_fetch($constraint)};
  }

  # If request biotype does not exist in the table
  # create a new biotype object containing name and object_type only
  # this is required by genebuild in pipelines
  if (!defined $biotype) {
    $biotype = Bio::EnsEMBL::Biotype->new(
          -NAME          => $name,
          -OBJECT_TYPE   => $object_type,
      )
  }

  return $biotype;
}

=head2 fetch_all_by_object_type

  Arg [1]    : String $object_type
               The object_type of the biotypes to retrieve (gene or transcript).
  Example    : $biotypes = $biotype_adaptor->fetch_all_by_object_type('gene');
  Description: Retrieves an array reference of biotype objects from the database.
  Returntype : arrayref of Bio::EnsEMBL::Biotype objects or empty arrayref
  Warning    : If empty arrayref is to be returned
  Exceptions : none

=cut

sub fetch_all_by_object_type {
  my ($self, $object_type) = @_;

  my $constraint = "b.object_type = ?";
  $self->bind_param_generic_fetch($object_type, SQL_VARCHAR);
  my @biotypes = @{$self->generic_fetch($constraint)};

  if ( !@biotypes ) {
    warning("No objects retrieved. Check if object_type '$object_type' is correct.")
  }

  return \@biotypes;
}

=head2 fetch_all_by_name

  Arg [1]    : String $name
               The name of the biotype to retrieve
  Arg [2]    : (optional) String $object_type
               The object_type of the biotypes to retrieve (gene or transcript).
  Example    : $biotypes = $biotype_adaptor->fetch_all_by_name('lincRNA');
  Description: Retrieves an array reference of biotype objects from the database.
  Returntype : arrayref of Bio::EnsEMBL::Biotype objects or empty arrayref
  Warning    : If empty arrayref is to be returned
  Exceptions : none

=cut

sub fetch_all_by_name{
  my ($self, $name, $object_type) = @_;

  my $constraint = "b.name = ?";
  $self->bind_param_generic_fetch($name, SQL_VARCHAR);
  if (defined $object_type) {
    $constraint .= "AND b.object_type = ?";
    $self->bind_param_generic_fetch($object_type, SQL_VARCHAR);
  }
  my @biotypes = @{$self->generic_fetch($constraint)};

  if ( !@biotypes ) {
    warning("No objects retrieved. Check if name '$name' is correct.")
  }

  return \@biotypes;
}

=head2 fetch_all_by_group_object_db_type

  Arg [1]    : String $biotype_group
               The group of the biotypes to retrieve
  Arg [2]    : String $object_type
               The object type of the biotypes to retrieve (gene or transcript)
  Arg [3]    : (optional) String $db_type
               The db_type of the biotypes to retrieve. If not provided defaults to 'core'.
  Example    : $biotype = $biotype_adaptor->fetch_all_by_group_object_db_type('coding', 'gene');
  Description: Retrieves an array reference of biotype objects from the database of the provided
               biotype_group and object_type and core db_type.
  Returntype : arrayref of Bio::EnsEMBL::Biotype objects or empty arrayref
  Warning    : If empty arrayref is to be returned
  Exceptions : none

=cut

sub fetch_all_by_group_object_db_type {
  my ($self, $biotype_group, $object_type, $db_type) = @_;

  $db_type //= 'core';
  $db_type = '%' . $db_type . '%';

  my $constraint = "b.biotype_group = ? AND b.db_type LIKE ? ";
  $self->bind_param_generic_fetch($biotype_group, SQL_VARCHAR);
  $self->bind_param_generic_fetch($db_type, SQL_VARCHAR);
  if (defined $object_type) {
    $constraint .= " AND b.object_type = ?";
    $self->bind_param_generic_fetch($object_type, SQL_VARCHAR);
  }
  my @biotypes = @{$self->generic_fetch($constraint)};

  if ( !@biotypes ) {
    warning("No objects retrieved. Check if object_type '$object_type' is correct.")
  }

  return \@biotypes;
}

1;
