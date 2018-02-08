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

Bio::EnsEMBL::DBSQL::BiotypeAdaptor - Encapsulates all generic access
to database meta information

=head1 SYNOPSIS

  my $meta_container = $db_adaptor->get_MetaContainer();

  my @mapping_info =
    @{ $meta_container->list_value_by_key('assembly.mapping') };

=head1 DESCRIPTION

  An object that encapsulates access to db meta data

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::BiotypeAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);
use Bio::EnsEMBL::Biotype;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 _tables

  Arg [1]    : none
  Description: PROTECTED implementation of superclass abstract method.
               Returns the names, aliases of the tables to use for queries.
  Returntype : list of listrefs of strings
  Exceptions : none
  Caller     : internal
  Status     : Stable

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
  Caller     : internal
  Status     : Stable

=cut

sub _columns {
  my $self = shift;

  return ('b.biotype_id', 'b.name', 'b.object_type', 'b.db_type', 'b.attrib_type_id', 'b.description', 'b.biotype_group', 'b.so_acc');
}

=head2 _objs_from_sth

  Arg [1]    : StatementHandle $sth
  Example    : none
  Description: PROTECTED implementation of abstract superclass method.
               responsible for the creation of ProteinFeatures
  Returntype : listref of Bio::EnsEMBL::Biotype
  Exceptions : none
  Caller     : internal
  Status     : At Risk

=cut

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my ($dbID, $name, $object_type, $db_type, $attrib_type_id, $description, $biotype_group, $so_acc);

  $sth->bind_columns(\$dbID, \$name, \$object_type, \$db_type, \$attrib_type_id, \$description, \$biotype_group, \$so_acc);

  my @biotypes;

  while($sth->fetch()) {
    push( @biotypes,
      my $feat = Bio::EnsEMBL::Biotype->new(
          -BIOTYPE_ID    => $dbID,
          -NAME          => $name,
          -OBJECT_TYPE   => $object_type,
          -DB_TYPE       => $db_type,
          -ATTRIB_TYPE_ID => $attrib_type_id,
          -DESCRIPTION   => $description,
          -BIOTYPE_GROUP => $biotype_group,
          -SO_ACC        => $so_acc,
      )
    );
  }

  return \@biotypes;
}


=head2 fetch_by_name_object_type

  Arg [1]    : String $name
               The name of the biotype to retrieve
  Arg [2]    : String $object_type
               The object type of the biotype to retrieve (gene or transcript)
  Example    : $biotype = $biotype_adaptor->fetch_by_name_object_type('gene', 'protein_coding');
  Description: Retrieves a biotype object from the database via its combined key (name, object_type).
               If the Biotype requested does not exist in the database, a new Biotype object is
               created with the provided name and object_type to be returned.
  Returntype : Bio::EnsEMBL::Biotype
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_name_object_type {
  my ($self, $name, $object_type) = @_;

  my $constraint = "b.name = ? AND b.object_type = ?";
  $self->bind_param_generic_fetch($name, SQL_VARCHAR);
  $self->bind_param_generic_fetch($object_type, SQL_VARCHAR);
  my ($biotype) = @{$self->generic_fetch($constraint)};

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
               listref of $sources
               The object_type of the biotypes to retrieve (gene or transcript).
  Example    : $biotypes = $biotype_adaptor->fetch_all_by_object_type('gene');
  Description: Retrieves an array reference of biotype objects from the database.
  Returntype : listref of Bio::EnsEMBL::Biotype objects or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_object_type {
  my ($self, $object_type) = @_;

  my $constraint = "b.object_type = ?";
  $self->bind_param_generic_fetch($object_type, SQL_VARCHAR);
  my @biotypes = @{$self->generic_fetch($constraint)};

  return \@biotypes;
}

1;