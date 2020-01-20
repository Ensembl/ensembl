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
use warnings;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::MicroRNA;
use Bio::EnsEMBL::RNAProduct;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Utils::Scalar qw( assert_ref );


use parent qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );



=head2 fetch_all_by_Transcript

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript
  Example    : $rps = $rnaproduct_adaptor->fetch_by_Transcript($transcript);
  Description: Retrieves RNAProducts via their associated transcript.
               If no RNAProducts are found, an empty list is returned.
  Returntype : arrayref of Bio::EnsEMBL::RNAProducts
  Exceptions : throw on incorrect argument
  Caller     : Transcript
  Status     : Stable

=cut

sub fetch_all_by_Transcript {
  my ($self, $transcript) = @_;

  assert_ref($transcript, 'Bio::EnsEMBL::Transcript');

  return $self->_fetch_direct_query(['rp.transcript_id', $transcript->dbID(), SQL_INTEGER]);
}


=head2 fetch_all_by_external_name

  Arg [1]    : String $external_name
               An external identifier of the RNAProduct to be obtained
  Arg [2]    : (optional) String $external_db_name
               The name of the external database from which the
               identifier originates.
  Arg [3]    : Boolean override. Force SQL regex matching for users
               who really do want to find all 'NM%'
  Example    : my @rnaproducts =
                  @{ $rp_a->fetch_all_by_external_name('MIMAT0000416') };
               my @more_rnaproducts =
                  @{ $rp_a->fetch_all_by_external_name('hsa-miR-1-__') };
  Description: Retrieves all RNAProducts which are associated with
               an external identifier such as a GO term, miRBase
               identifer, etc. Usually there will only be a single
               RNAProduct returned in the list reference, but not
               always. If no RNAProducts with the external identifier
               are found, a reference to an empty list is returned.
               SQL wildcards % and _ are supported in the $external_name
               but their use is somewhat restricted for performance reasons.
               Users that really do want % and _ in the first three characters
               should use argument 3 to prevent optimisations
  Returntype : listref of Bio::EnsEMBL::RNAProduct
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_external_name {
  my ($self, $external_name, $external_db_name, $override) = @_;

  my $entry_adaptor = $self->db->get_DBEntryAdaptor();

  my @ids = $entry_adaptor->list_rnaproduct_ids_by_extids($external_name,
                                                          $external_db_name,
                                                          $override);

  my $transcript_adaptor = $self->db()->get_TranscriptAdaptor();

  my @reference;
  my @non_reference;
  foreach my $id (@ids) {
    my $transcript = $transcript_adaptor->fetch_by_rnaproduct_id($id);

    if ( defined($transcript) ) {
      my $rnaproduct = $self->fetch_by_dbID($id);
      if ( $transcript->slice()->is_reference() ) {
        push(@reference, $rnaproduct);
      }
      else {
        push(@non_reference, $rnaproduct);
      }
    }
  }

  return [@reference, @non_reference];
}


=head2 fetch_all_by_type

  Arg [1]    : string $type_code
  Example    : $rps = $rp_a->fetch_all_by_type('miRNA');
  Description: Retrieves RNAProducts via their type (e.g. miRNA, circRNA).
               If no matching RNAProducts are found, an empty list is
               returned.
  Returntype : arrayref of Bio::EnsEMBL::RNAProducts
  Exceptions : throws if type code is undefined
  Caller     : ?
  Status     : In Development

=cut

sub fetch_all_by_type {
  my ($self, $type_code) = @_;

  if ( !defined $type_code ) {
    throw("type code argument is required");
  }

  return ($self->_fetch_direct_query(['pt.code', $type_code, SQL_VARCHAR]));
}


=head2 fetch_by_dbID

  Arg [1]    : int $dbID
               The internal identifier of the RNAProduct to obtain
  Example    : $rnaproduct = $rnaproduct_adaptor->fetch_by_dbID(1234);
  Description: This fetches a RNAProduct object via its internal id.
               This is only debatably useful since RNAProducts do
               not make much sense outside of the context of their
               Transcript.  Consider using fetch_by_Transcript instead.
  Returntype : Bio::EnsEMBL::RNAProduct, or undef if the RNAProduct is not
               found.
  Caller     : ?
  Status     : Stable

=cut

sub fetch_by_dbID {
  my ($self, $dbID) = @_;

  if ( !defined $dbID ) {
    throw("dbID argument is required");
  }

  return ($self->_fetch_direct_query(['rp.rnaproduct_id', $dbID, SQL_INTEGER]))->[0];
}


=head2 fetch_by_stable_id

  Arg [1]    : string $stable_id
               The stable identifier of the RNAProduct to obtain
  Example    : $rnaproduct = $rnaproduct_adaptor->fetch_by_stable_id("ENSS00001");
  Description: This fetches a RNAProduct object via its stable id.
  Returntype : Bio::EnsEMBL::RNAProduct, or undef if the RNAProduct is not
               found.
  Caller     : ?
  Status     : Stable

=cut

sub fetch_by_stable_id {
  my ($self, $stable_id) = @_;

  if ( !defined $stable_id ) {
    throw("stable id argument is required");
  }

  return ($self->_fetch_direct_query(['rp.stable_id', $stable_id, SQL_VARCHAR]))->[0];
}


=head2 list_dbIDs

  Arg [1]    : none
  Example    : @rnaproduct_ids = @{$rnaproduct_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all RNAProducts in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub list_dbIDs {
   my ($self) = @_;

   return $self->_list_dbIDs("rnaproduct");
}


=head2 remove

  Arg [1]    : Bio::EnsEMBL::RNAProduct $rnaproduct
               The RNAProduct to be removed from the database
  Example    : $rpID = $rp_adaptor->remove($rnaproduct, $transcript->dbID());
  Description: Removes a RNAProduct, along with all associated information
               from the database.
  Returntype : none
  Exceptions : throw on incorrect arguments
  Caller     : ?
  Status     : Stable

=cut

sub remove {
  my ($self, $rnaproduct) = @_;

  if (!ref($rnaproduct) || !$rnaproduct->isa('Bio::EnsEMBL::RNAProduct')) {
    throw("$rnaproduct is not a EnsEMBL rnaproduct");
  }

  my $db = $self->db();

  # Do nothing if the object is not stored to begin with
  if (!$rnaproduct->is_stored($db)) {
    return;
  }

  # Remove xrefs
  my $dbe_adaptor = $db->get_DBEntryAdaptor();
  for my $dbe (@{ $rnaproduct->get_all_DBEntries() }) {
    $dbe_adaptor->remove_from_object($dbe, $rnaproduct, 'RNAProduct');
  }

  # Remove attributes
  my $attr_adaptor = $db->get_AttributeAdaptor();
  $attr_adaptor->remove_from_RNAProduct($rnaproduct);

  # Remove rnaproduct itself
  my $sth = $self->prepare("DELETE FROM rnaproduct WHERE rnaproduct_id = ?");
  $sth->bind_param(1, $rnaproduct->dbID(), SQL_INTEGER);
  $sth->execute();
  $sth->finish();

  # Mark the object as local
  $rnaproduct->dbID(undef);
  $rnaproduct->adaptor(undef);

  return;
}


=head2 store

  Arg [1]    : Bio::EnsEMBL::RNAProduct $rnaproduct
               The RNAProduct to be written to the database
  Arg [2]    : Int $transcript_dbID
               The identifier of the transcript that this RNAProduct is
               associated with
  Example    : $rpID = $rp_adaptor->store($rnaproduct, $transcript->dbID());
  Description: Stores a RNAProduct in the database and returns the new
               internal identifier for the stored RNAProduct.
  Returntype : Int
  Exceptions : throw on incorrect arguments
  Caller     : general
  Status     : Stable

=cut

sub store {
  my ($self, $rnaproduct, $transcript_dbID) = @_;

  if (!ref($rnaproduct) || !$rnaproduct->isa('Bio::EnsEMBL::RNAProduct')) {
    throw("$rnaproduct is not a EnsEMBL rnaproduct - not storing");
  }

  my $db = $self->db();

  # Avoid creating duplicate entries
  if ($rnaproduct->is_stored($db)) {
    return $rnaproduct->dbID();
  }

  my ( $start_exon_dbID, $end_exon_dbID );
  if ( $rnaproduct->start_Exon() ) {
    $start_exon_dbID = $rnaproduct->start_Exon()->dbID();
  }
  if ( $rnaproduct->end_Exon() ) {
    $end_exon_dbID = $rnaproduct->end_Exon()->dbID();
  }

  # Store rnaproduct

  my @columns = qw{ transcript_id seq_start start_exon_id seq_end end_exon_id };

  my @canned_columns = ( 'rnaproduct_type_id' );
  my @canned_values
    = ( '(SELECT rnaproduct_type_id FROM rnaproduct_type WHERE code=?)' );

  if (defined($rnaproduct->stable_id())) {
    push @columns, 'stable_id', 'version';

    my $created = $db->dbc->from_seconds_to_date($rnaproduct->created_date());
    if ($created) {
      push @canned_columns, 'created_date';
      push @canned_values,  $created;
    }

    my $modified = $db->dbc->from_seconds_to_date($rnaproduct->modified_date());
    if ($modified) {
      push @canned_columns, 'modified_date';
      push @canned_values,  $modified;
    }
  }

  my $column_string = join(', ', @columns, @canned_columns);
  my $value_string  = join(', ', (q{?}) x @columns, @canned_values);
  my $store_rnaproduct_sql
    = "INSERT INTO rnaproduct ( $column_string ) VALUES ( $value_string )";
  my $rp_st = $self->prepare($store_rnaproduct_sql);

  $rp_st->bind_param(  1, $transcript_dbID,         SQL_INTEGER);
  $rp_st->bind_param(  2, $rnaproduct->start(),     SQL_INTEGER);
  $rp_st->bind_param(  3, $start_exon_dbID,         SQL_INTEGER);
  $rp_st->bind_param(  4, $rnaproduct->end(),       SQL_INTEGER);
  $rp_st->bind_param(  5, $end_exon_dbID,           SQL_INTEGER);
  if (defined($rnaproduct->stable_id())) {
    $rp_st->bind_param(6, $rnaproduct->stable_id(), SQL_VARCHAR);
    $rp_st->bind_param(7, $rnaproduct->version(),   SQL_INTEGER);
  }
  $rp_st->bind_param(  8, $rnaproduct->type_code(), SQL_VARCHAR);
  $rp_st->execute();
  $rp_st->finish();

  # Retrieve the newly assigned dbID
  my $rp_dbID = $self->last_insert_id('rnaproduct_id', undef, 'rnaproduct');

  # Store attributes
  $rnaproduct->synchronise_attributes();
  my $attr_adaptor = $db->get_AttributeAdaptor();
  $attr_adaptor->store_on_RNAProduct($rp_dbID,
                                     $rnaproduct->get_all_Attributes());

  # Store xrefs
  my $dbe_adaptor = $db->get_DBEntryAdaptor();
  for my $dbe (@{ $rnaproduct->get_all_DBEntries() }) {
    $dbe_adaptor->store($dbe, $rp_dbID, 'RNAProduct', 1);
  }

  # Link the rnaproduct object to its data row
  $rnaproduct->dbID($rp_dbID);
  $rnaproduct->adaptor($self);

  return $rp_dbID;
}


=head2 _list_dbIDs

  Arg[1]      : String $table
  Arg[2]      : String $column
  Example     : $rnaproduct_adaptor->_list_dbIDs('rnaproduct', 'rnaproduct_id');
  Description : Local reimplementation to ensure multi-species RNAProducts
                are limited to their species alone
  Returntype  : ArrayRef of specified IDs
  Caller      : Internal
  Status      : Unstable

=cut

sub _list_dbIDs {
  my ($self, $table, $column) = @_;
  my $ids;
  if($self->is_multispecies()) {
    $column ||= "${table}_id";
    my $sql = <<"SQL";
SELECT `rp`.`${column}` FROM rnaproduct rp
JOIN transcript t USING (transcript_id)
JOIN seq_region sr USING (seq_region_id)
JOIN coord_system cs USING (coord_system_id)
WHERE cs.species_id = ?
SQL
    return $self->dbc()->sql_helper()->execute_simple(-SQL => $sql, -PARAMS => [$self->species_id()]);
  }
  else {
    $ids = $self->SUPER::_list_dbIDs($table, $column);
  }
  return $ids;
}


=head2 _fetch_direct_query
  Arg [1]    : reference to an array consisting of:
                - the name of the column to use in the WHERE clause,
                - the value fields from that column are to be equal to,
                - the data type of that column (e.g. SQL_INTEGER)
  Description: PRIVATE internal method shared between public fetch methods
               in order to avoid duplication of SQL-query logic. NOT SAFE
               to be handed directly to users because in its current form
               it can be trivially exploited to inject arbitrary SQL.
  Returntype : ArrayRef of either Bio::EnsEMBL::RNAProducts or undefs
  Exceptions : throws if RNAProduct type is absent or unknown
  Caller     : internal
  Status     : At Risk (In Development)

=cut

sub _fetch_direct_query {
  my ($self, $where_args) = @_;

  my $rp_created_date =
    $self->db()->dbc()->from_date_to_seconds('created_date');
  my $rp_modified_date =
    $self->db()->dbc()->from_date_to_seconds('modified_date');

  my $sql_template = <<'SQL';
SELECT rp.rnaproduct_id, pt.code, rp.transcript_id, rp.seq_start,
  rp.start_exon_id, rp.seq_end, rp.end_exon_id, rp.stable_id, rp.version,
  %s, %s
 FROM rnaproduct rp
 JOIN rnaproduct_type pt ON rp.rnaproduct_type_id = pt.rnaproduct_type_id
 WHERE %s = ?
SQL
  my $sql = sprintf( $sql_template,
                     $rp_created_date, $rp_modified_date, $where_args->[0] );
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $where_args->[1], $where_args->[2]);
  $sth->execute();
  my $results = $self->_obj_from_sth($sth);
  $sth->finish();

  return $results;
}


=head2 _obj_from_sth
  Arg [1]    : DBI statement handle
  Description: PRIVATE internal method shared between public SQL-query
               methods in order to avoid duplication of object-creation
               logic.
  Returntype : ArrayRef of either Bio::EnsEMBL::RNAProducts or undefs
  Exceptions : throws if RNAProduct type is absent or unknown
  Caller     : internal
  Status     : At Risk (In Development)

=cut

sub _obj_from_sth {
  my ($self, $sth) = @_;
  my @return_data;

  my $sql_data = $sth->fetchall_arrayref();
  my $transcript_adaptor = $self->db()->get_TranscriptAdaptor();
  while (my $row_ref = shift @{$sql_data}) {
    my ($rnaproduct_id, $type_code, $transcript_id, $seq_start, $start_exon_id,
        $seq_end, $end_exon_id, $stable_id, $version, $created_date,
        $modified_date) = @{$row_ref};

    if (!defined($rnaproduct_id)) {
      push @return_data, undef;
      next;
    }

    my $exon_adaptor = $self->db()->get_ExonAdaptor();
    my ( $start_exon, $end_exon );
    if (defined $start_exon_id ) {
      $start_exon   = $exon_adaptor->fetch_by_dbID( $start_exon_id );
    }
    if (defined $end_exon_id ) {
      $end_exon     = $exon_adaptor->fetch_by_dbID( $end_exon_id );
    }

    my $class_name = Bio::EnsEMBL::Utils::RNAProductTypeMapper::mapper()
      ->type_code_to_class($type_code);
    my $rnaproduct = $class_name->new_fast( {
      'dbID'          => $rnaproduct_id,
      'type_code'     => $type_code,
      'adaptor'       => $self,
      'start'         => $seq_start,
      'end'           => $seq_end,
      'start_exon'    => $start_exon,
      'end_exon'      => $end_exon,
      'stable_id'     => $stable_id,
      'version'       => $version,
      'created_date'  => $created_date || undef,
      'modified_date' => $modified_date || undef,
    } );

    my $transcript = $transcript_adaptor->fetch_by_dbID($transcript_id);
    $rnaproduct->transcript($transcript);

    push @return_data, $rnaproduct;
  }

  return \@return_data;
}


=head2 _tables

  Arg [1]    : none
  Description: PROTECTED implementation of superclass abstract method.
               Returns the names, aliases of the tables to use for queries.
  Returntype : list of listrefs of strings
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

# This method IS used, at the superclass level
sub _tables { ## no critic (Subroutines::ProhibitUnusedPrivateSubroutines)
  return (['rnaproduct', 'rp']);
}


1;
