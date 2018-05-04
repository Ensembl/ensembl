=head1 LICENSE

Copyright [2018] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::DBSQL::RNAProductAdaptor - Provides a means to fetch and store
RNAProduct objects from a database.

=head1 DESCRIPTION

This adaptor provides a means to retrieve and store
Bio::EnsEMBL::RNAProduct objects from/in a database.

RNAProduct objects only truly make sense in the context of their
transcripts so the recommended means to retrieve RNAProducts is
by retrieving the Transcript object first, and then fetching the
RNAProduct.

=head1 SYNOPSIS

  use Bio::EnsEMBL::Registry;

  Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
  );

  $rnaproduct_adaptor =
    Bio::EnsEMBL::Registry->get_adaptor( "human", "core",
    "rnaproduct" );

  ...

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::RNAProductAdaptor;

use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::RNAProduct;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Utils::Scalar qw( assert_ref );


use parent qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );



=head2 list_dbIDs

  Arg [1]    : none
  Example    : @rnaproduct_ids = @{$rnaproduct_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all rnaproducts in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub list_dbIDs {
   my ($self) = @_;

   return $self->_list_dbIDs("rnaproduct");
}


=head2 _list_dbIDs

  Arg[1]      : String $table
  Arg[2]      : String $column
  Example     : $rnaproduct_adaptor->_list_dbIDs('rnaproduct', 'rnaproduct_id');
  Description : Local reimplementation to ensure multi-species rnaproducts
                are limited to their species alone
  Returntype  : ArrayRef of specified IDs
  Caller      : Internal
  Status      : Unstable
=cut

sub _list_dbIDs {
  my ($self, $table, $column) = @_;
  my $ids;
  if($self->is_multispecies()) {
    # FIXME: test this, it might not have been fully converted from TranslationAdaptor yet
    $column ||= "${table}_id";
    my $sql = <<SQL;
select `rp`.`${column}`
from rnaproduct rp
join transcript t using (transcript_id)
join seq_region sr using (seq_region_id)
join coord_system cs using (coord_system_id)
where cs.species_id =?
SQL
    return $self->dbc()->sql_helper()->execute_simple(-SQL => $sql, -PARAMS => [$self->species_id()]);
  }
  else {
    $ids = $self->SUPER::_list_dbIDs($table, $column);
  }
  return $ids;
}


# _tables
#  Arg [1]    : none
#  Description: PROTECTED implementation of superclass abstract method.
#               Returns the names, aliases of the tables to use for queries.
#  Returntype : list of listrefs of strings
#  Exceptions : none
#  Caller     : internal
#  Status     : Stable

sub _tables {
  return (['rnaproduct', 'rp']);
}

1;
