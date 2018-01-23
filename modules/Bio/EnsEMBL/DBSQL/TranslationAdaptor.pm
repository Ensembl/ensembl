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

Bio::EnsEMBL::DBSQL::TranslationAdaptor - Provides a means to fetch and store
Translation objects from a database.

=head1 DESCRIPTION

This adaptor provides a means to retrieve and store
Bio::EnsEMBL::Translation objects from/in a database.

Translation objects only truly make sense in the context of their
transcripts so the recommended means to retrieve Translations is
by retrieving the Transcript object first, and then fetching the
Translation.

=head1 SYNOPSIS

  use Bio::EnsEMBL::Registry;

  Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
  );

  $transcript_adaptor =
    Bio::EnsEMBL::Registry->get_adaptor( "human", "core",
    "transcript" );

  $translation_adaptor =
    Bio::EnsEMBL::Registry->get_adaptor( "human", "core",
    "translation" );

  my $transcript = $transcript_adaptor->fetch_by_dbID(131243);
  my $translation =
    $translation_adaptor->fetch_by_Transcript($transcript);

  print("Translation Start Site: "
      . $translation->start_Exon()->stable_id() . " "
      . $translation->start()
      . "\n" );
  print("Translation Stop: "
      . $translation->end_Exon()->stable_id() . " "
      . $translation->end() );

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::TranslationAdaptor;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate );
use Bio::EnsEMBL::Utils::Scalar qw( assert_ref );


@ISA = qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );

=head2 fetch_all_alternative_by_Transcript

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript
  Example    :

    @tl = @{
      $translation_adaptor->fetch_all_alternative_by_Transcript(
                                                            $transcript)
      };

  Description: Retrieves all alternative translations associated with a
               particular transcript.  If no alternative translation is
               found, a reference to an empty list is returned.

  Returntype : listref of Bio::EnsEMBL::Translation
  Exceptions : throw on incorrect argument
  Caller     : Transcript
  Status     : Stable

=cut

sub fetch_all_alternative_by_Transcript {
  my ( $self, $transcript ) = @_;

  assert_ref($transcript, 'Bio::EnsEMBL::Transcript');

  my $tl_created_date =
    $self->db()->dbc()->from_date_to_seconds('tl.created_date');
  my $tl_modified_date =
    $self->db()->dbc()->from_date_to_seconds('tl.modified_date');

  my $sql =
    sprintf( "SELECT tl.translation_id, tl.start_exon_id, "
      . "tl.end_exon_id, tl.seq_start, tl.seq_end, "
      . "tl.stable_id, tl.version, %s, %s "
      . "FROM translation tl "
      . "JOIN transcript t "
      . "ON (t.transcript_id = tl.transcript_id) "
      . "WHERE tl.transcript_id = ? "
      . "AND tl.translation_id != t.canonical_translation_id",
    $tl_created_date, $tl_modified_date );

  my $transcript_id = $transcript->dbID();
  my $sth           = $self->prepare($sql);
  $sth->bind_param( 1, $transcript_id, SQL_INTEGER );

  $sth->execute();

  my (
    $translation_id, $start_exon_id, $end_exon_id,
    $seq_start,      $seq_end,       $stable_id,
    $version,        $created_date,  $modified_date
  );

  $sth->bind_columns(
    \(
      $translation_id, $start_exon_id, $end_exon_id,
      $seq_start,      $seq_end,       $stable_id,
      $version,        $created_date,  $modified_date
    ) );

  # Get all alternative translations.
  my $translations = [];
  while ( $sth->fetch() ) {
    if ( !defined($translation_id) ) { next }

    my ( $start_exon, $end_exon );

    # this will load all the exons whenever we load the translation
    # but I guess thats ok ....

    foreach my $exon ( @{ $transcript->get_all_Exons() } ) {
      if ( $exon->dbID() == $start_exon_id ) { $start_exon = $exon }
      if ( $exon->dbID() == $end_exon_id )   { $end_exon   = $exon }
    }

    if ( !( defined($start_exon) && defined($end_exon) ) ) {
      throw(
        sprintf(
          "Could not find start or end exon in transcript_id=%d\n",
          $transcript->dbID() ) );
    }

    my $translation =
      Bio::EnsEMBL::Translation->new_fast( {
                             'dbID'          => $translation_id,
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

    $translation->transcript($transcript);

    push( @{$translations}, $translation );

  } ## end while ( $sth->fetch() )

  return $translations;
} ## end sub fetch_all_by_Transcript

=head2 fetch_by_Transcript

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript
  Example    : $tl = $translation_adaptor->fetch_by_Transcript($transcript);
  Description: Retrieves a Translation via its associated transcript.
               If the Translation is not found, undef is returned.
  Returntype : Bio::EnsEMBL::Translation
  Exceptions : throw on incorrect argument
  Caller     : Transcript
  Status     : Stable

=cut

sub fetch_by_Transcript {
  my ( $self, $transcript ) = @_;

  assert_ref( $transcript, 'Bio::EnsEMBL::Transcript' );

  my $tl_created_date =
    $self->db()->dbc()->from_date_to_seconds('tl.created_date');
  my $tl_modified_date =
    $self->db()->dbc()->from_date_to_seconds('tl.modified_date');

  my $sql =
    sprintf( "SELECT tl.translation_id, tl.start_exon_id, "
      . "tl.end_exon_id, tl.seq_start, tl.seq_end, "
      . "tl.stable_id, tl.version, %s, %s "
      . "FROM translation tl "
      . "JOIN transcript tr "
      . "ON (tl.translation_id = tr.canonical_translation_id) "
      . "WHERE tr.transcript_id = ?",
    $tl_created_date, $tl_modified_date );

  my $transcript_id = $transcript->dbID();
  my $sth           = $self->prepare($sql);
  $sth->bind_param( 1, $transcript_id, SQL_INTEGER );

  $sth->execute();

  my (
    $translation_id, $start_exon_id, $end_exon_id,
    $seq_start,      $seq_end,       $stable_id,
    $version,        $created_date,  $modified_date
  ) = $sth->fetchrow_array();
  $sth->finish();

  if ( !defined($translation_id) ) { return undef }

  my ( $start_exon, $end_exon );

  # this will load all the exons whenever we load the translation
  # but I guess thats ok ....

  foreach my $exon ( @{ $transcript->get_all_Exons() } ) {
    if ( $exon->dbID() == $start_exon_id ) { $start_exon = $exon }
    if ( $exon->dbID() == $end_exon_id )   { $end_exon   = $exon }
  }

  if ( !( defined($start_exon) && defined($end_exon) ) ) {
    throw(
      sprintf( "Could not find start or end exon in transcript_id=%d\n",
        $transcript->dbID() ) );
  }

  my $translation =
    Bio::EnsEMBL::Translation->new_fast( {
                             'dbID'          => $translation_id,
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

  $translation->transcript($transcript);

  return $translation;
} ## end sub fetch_by_Transcript



=head2 fetch_all_by_external_name

  Arg [1]    : string $external_name
               The external identifier for the translation(s) to be
               obtained.
  Arg [2]    : (optional) string $external_db_name
               The name of the external database from which the
               identifier originates.
  Arg [3]    : Boolean override. Force SQL regex matching for users
               who really do want to find all 'NM%'
  Example    : my @translations =
                  @{ $trl_adaptor->fetch_all_by_external_name('BRCA2') };
               my @many_translations = 
                  @{ $trl_adaptor->fetch_all_by_external_name('BRCA%') };
  Description: Retrieves a list of translations fetched via an
               external identifier.  Note that this may not be a
               particularly useful method, because translations
               do not make much sense out of the context of
               their transcript.  It may be better to use the
               TranscriptAdaptor::fetch_all_by_external_name instead.
               SQL wildcards % and _ are supported in the $external_name
               but their use is somewhat restricted for performance reasons.
               Users that really do want % and _ in the first three characters
               should use argument 3 to prevent optimisations
  Returntype : reference to a list of Translations
  Exceptions : none
  Caller     : general
  Status     : Medium Risk
             :   At some time may be deprecated to instead use 
             :   TranscriptAdaptor::fetch_all_by_external_name 

=cut

sub fetch_all_by_external_name {
  my ( $self, $external_name, $external_db_name, $override ) = @_;

  my $entry_adaptor = $self->db->get_DBEntryAdaptor();

  my @ids = $entry_adaptor->list_translation_ids_by_extids( 
            $external_name, $external_db_name, $override );

  my $transcript_adaptor = $self->db()->get_TranscriptAdaptor();

  my @reference;
  my @non_reference;
  foreach my $id (@ids) {
    my $transcript = $transcript_adaptor->fetch_by_translation_id($id);

    if ( defined($transcript) ) {
      my $translation = $self->fetch_by_Transcript($transcript);
      if($transcript->slice()->is_reference()) {
        push(@reference, $translation);
      }
      else {
        push(@non_reference, $translation);
      }
    }
  }

  return [@reference, @non_reference];
}

=head2 fetch_all_by_GOTerm

  Arg [1]   : Bio::EnsEMBL::OntologyTerm
              The GO term for which translations should be fetched.

  Example:  @translations = @{
              $translation_adaptor->fetch_all_by_GOTerm(
                $go_adaptor->fetch_by_accession('GO:0030326') ) };

  Description   : Retrieves a list of translations that are
                  associated with the given GO term, or with any of
                  its descendent GO terms.

  Return type   : listref of Bio::EnsEMBL::Translation
  Exceptions    : Throws of argument is not a GO term
  Caller        : general
  Status        : Stable

=cut

sub fetch_all_by_GOTerm {
  my ( $self, $term ) = @_;

  assert_ref( $term, 'Bio::EnsEMBL::OntologyTerm' );
  if ( $term->ontology() ne 'GO' ) {
    throw('Argument is not a GO term');
  }

  my $entryAdaptor = $self->db->get_DBEntryAdaptor();

  my %unique_dbIDs;
  foreach my $accession ( map { $_->accession() }
                          ( $term, @{ $term->descendants() } ) )
  {
    my @ids =
      $entryAdaptor->list_translation_ids_by_extids( $accession, 'GO' );
    foreach my $dbID (@ids) { $unique_dbIDs{$dbID} = 1 }
  }

  my @result;
  if ( scalar( keys(%unique_dbIDs) ) > 0 ) {
    my $transcript_adaptor = $self->db()->get_TranscriptAdaptor();

    foreach my $dbID ( sort { $a <=> $b } keys(%unique_dbIDs) ) {
      my $transcript =
        $transcript_adaptor->fetch_by_translation_id($dbID);
      if ( defined($transcript) ) {
        push( @result, $self->fetch_by_Transcript($transcript) );
      }
    }
  }

  return \@result;
} ## end sub fetch_all_by_GOTerm

=head2 fetch_all_by_GOTerm_accession

  Arg [1]   : String
              The GO term accession for which genes should be
              fetched.

  Example   :

    @genes =
      @{ $gene_adaptor->fetch_all_by_GOTerm_accession('GO:0030326') };

  Description   : Retrieves a list of genes that are associated with
                  the given GO term, or with any of its descendent
                  GO terms.  The genes returned are in their native
                  coordinate system, i.e. in the coordinate system
                  in which they are stored in the database.  If
                  another coordinate system is required then the
                  Gene::transfer or Gene::transform method can be
                  used.

  Return type   : listref of Bio::EnsEMBL::Gene
  Exceptions    : Throws of argument is not a GO term accession
  Caller        : general
  Status        : Stable

=cut

sub fetch_all_by_GOTerm_accession {
  my ( $self, $accession ) = @_;

  if ( $accession !~ /^GO:/ ) {
    throw('Argument is not a GO term accession');
  }

  my $goAdaptor =
    Bio::EnsEMBL::Registry->get_adaptor( 'Multi', 'Ontology',
                                         'OntologyTerm' );

  my $term = $goAdaptor->fetch_by_accession($accession);

  return $self->fetch_all_by_GOTerm($term);
}

=head2 store

  Arg [1]    : Bio::EnsEMBL::Translation $translation
               The translation object to be stored in the database 
  Example    : $transl_id = $translation_adaptor->store($translation);
  Description: Stores a translation object in the database
  Returntype : int - the new dbID of the stored translation
  Exceptions : thrown if the dbID of the start_Exon or end_Exon is not 
               defined.
               thrown if only partial stable id information is present (e.g.
               identifier but not version number)
  Caller     : Transcript::store
  Status     : Stable

=cut

sub store {
  my ( $self, $translation, $transcript_id )  = @_;
  
  my $start_exon = $translation->start_Exon();
  my $end_exon   = $translation->end_Exon();
 
  if(!$start_exon) {
    throw("Translation must define a start_Exon to be stored.");
  }
 
  if(!$end_exon) {
    throw("Translation must define an end_Exon to be stored.");
  }
 
  if(!$start_exon->dbID) {
    throw("start_Exon must have a dbID for Translation to be stored.");
  }

  if(!$end_exon->dbID) {
    throw("end_Exon must have a dbID for Translation to be stored.");
  }

  my @columns = qw(
                 seq_start
                 start_exon_id
                 seq_end
                 end_exon_id
                 transcript_id
  );

  my @canned_columns;
  my @canned_values;

  if (defined($translation->stable_id)) {
      push @columns, 'stable_id', 'version';

      my $created = $self->db->dbc->from_seconds_to_date($translation->created_date());
      my $modified = $self->db->dbc->from_seconds_to_date($translation->modified_date());

      if ($created) {
	push @canned_columns, 'created_date';
	push @canned_values,  $created;
      }
      if ($modified) {
	push @canned_columns, 'modified_date';
	push @canned_values,  $modified;
      }
  }

  my $columns = join(', ', @columns, @canned_columns);
  my $values  = join(', ', ('?') x @columns, @canned_values);
  my $store_translation_sql = qq(
        INSERT INTO translation ( $columns ) VALUES ( $values )
  );

  my $sth = $self->prepare($store_translation_sql);
  $sth->bind_param(1,$translation->start,SQL_INTEGER);
  $sth->bind_param(2,$translation->start_Exon->dbID,SQL_INTEGER);
  $sth->bind_param(3,$translation->end,SQL_INTEGER);
  $sth->bind_param(4,$translation->end_Exon->dbID,SQL_INTEGER);
  $sth->bind_param(5,$transcript_id,SQL_INTEGER);


  if (defined($translation->stable_id)) {
 
    $sth->bind_param(6, $translation->stable_id,SQL_VARCHAR);
    $sth->bind_param(7, $translation->version,SQL_VARCHAR);
  }

  $sth->execute();
 
  my $transl_dbID = $self->last_insert_id('translation_id', undef, 'translation');

  #
  # store object xref mappings to translations
  #
 
  my $dbEntryAdaptor = $self->db()->get_DBEntryAdaptor();
  # store each of the xrefs for this translation
  foreach my $dbl ( @{$translation->get_all_DBEntries} ) {
     $dbEntryAdaptor->store( $dbl, $transl_dbID, "Translation", 1 );
  }

  #storing the protein features associated with the translation
  my $pfadaptor = $self->db->get_ProteinFeatureAdaptor();
  foreach my $pf(@{$translation->get_all_ProteinFeatures}){
    $pfadaptor->store($pf, $transl_dbID);
  }

  $translation->get_all_Attributes();

  # store any translation attributes that are defined
  my $attr_adaptor = $self->db->get_AttributeAdaptor();
  $attr_adaptor->store_on_Translation($transl_dbID,
                                      $translation->get_all_Attributes());

  $translation->dbID($transl_dbID);
  $translation->adaptor($self);

  return $transl_dbID;
}



=head2 remove

  Arg [1]    : Bio::EnsEMBL::Translation $translation
  Example    : $translation_adaptor->remove($translation);
  Description: Removes a translation completely from the database, and all
               associated information including protein features etc.
  Returntype : none
  Exceptions : throw on incorrect arguments
               warning if translation is not in this database
  Caller     : TranscriptAdaptor::remove
  Status     : Stable

=cut

sub remove {
  my $self = shift;
  my $translation = shift;

  if(!ref($translation) || !$translation->isa('Bio::EnsEMBL::Translation')) {
    throw("Bio::EnsEMBL::Translation argument expected.");
  }

  if( !$translation->is_stored($self->db()) ) {
    warning("Cannot remove translation " . $translation->dbID() . 
            ". Is not stored in this database.");
    return;
  }

  # remove athe attributes associated with this translation
  my $attrib_adp = $self->db->get_AttributeAdaptor;
  $attrib_adp->remove_from_Translation($translation);

  # remove all xref associations to this translation
  my $dbe_adaptor = $self->db()->get_DBEntryAdaptor();
  foreach my $dbe (@{$translation->get_all_DBEntries()}) {
    $dbe_adaptor->remove_from_object($dbe, $translation, 'Translation');
  }

  # remove all protein_features on this translation
  my $sth = $self->prepare
    ("DELETE FROM protein_feature WHERE translation_id = ?");
  $sth->bind_param(1,$translation->dbID,SQL_INTEGER);
  $sth->execute();
  $sth->finish();

  # remove the translation itself

  $sth = $self->prepare("DELETE FROM translation WHERE translation_id = ?" );
  $sth->bind_param(1,$translation->dbID,SQL_INTEGER);
  $sth->execute();
  $sth->finish();

  $translation->dbID( undef );
  $translation->adaptor(undef);

  return
}


=head2 list_dbIDs

  Arg [1]    : none
  Example    : @translation_ids = @{$translation_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all translations in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub list_dbIDs {
   my ($self) = @_;

   return $self->_list_dbIDs("translation");
}


=head2 list_stable_ids

  Arg [1]    : none
  Example    : @transl_stable_ids = @{$transl_adaptor->list_stable_dbIDs()};
  Description: Gets an array of stable ids for all translations in the current 
               db
  Returntype : reference to a list of strings
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub list_stable_ids {
   my ($self) = @_;

   return $self->_list_dbIDs("translation", "stable_id");
}

=head2 _list_dbIDs

  Arg[1]      : String $table
  Arg[2]      : String $column
  Example     : $transl_adaptor->_list_dbIDs('translation','translation_id');
  Description : Local reimplementation to ensure multi-species translations
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
    my $sql = <<SQL;
select `tr`.`${column}` 
from translation tr
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



=head2 fetch_by_dbID

  Arg [1]    : int $dbID
               The internal identifier of the Translation to obtain
  Example    : $translation = $translation_adaptor->fetch_by_dbID(1234);
  Description: This fetches a Translation object via its internal id.
               This is only debatably useful since translations do
               not make much sense outside of the context of their
               Transcript.  Consider using fetch_by_Transcript instead.
  Returntype : Bio::EnsEMBL::Translation, or undef if the translation is not
               found.
  Exceptions : warning if an additional (old style) Transcript argument is
               provided
  Caller     : ?
  Status     : Stable

=cut

sub fetch_by_dbID {
  my ( $self, $dbID, $transcript ) = @_;

  if ($transcript) {
    deprecate(   "Use of fetch_by_dbID "
               . "with a Transcript argument is deprecated."
               . "Use fetch_by_Transcript instead." );
  }

  if ( !defined($dbID) ) {
    throw("dbID argument is required");
  }

  my $transcript_adaptor = $self->db()->get_TranscriptAdaptor();
  $transcript = $transcript_adaptor->fetch_by_translation_id($dbID);

  if ( defined($transcript) ) {
    my $translation = $self->fetch_by_Transcript($transcript);

    if ( defined($translation) && $translation->dbID()==$dbID ) {
      return $translation;
    }

    my @alt_translations =
      @{ $self->fetch_all_alternative_by_Transcript($transcript) };

    foreach my $alt_translation (@alt_translations) {
      if ( $alt_translation->dbID() == $dbID ) {
        return $alt_translation;
      }
    }
  }

  return undef;
} ## end sub fetch_by_dbID


=head2 fetch_by_stable_id

  Arg [1]    : string $stable_id
               The stable identifier of the Translation to obtain
  Example    : $translation = $translation_adaptor->fetch_by_stable_id("ENSP00001");
  Description: This fetches a Translation object via its stable id.
               This is only debatably useful since translations do
               not make much sense outside of the context of their
               Transcript.  Consider using fetch_by_Transcript instead.
  Returntype : Bio::EnsEMBL::Translation or undef if the translation is not
               found.
  Exceptions : warning if an additional (old style) Transcript argument is
               provided
  Caller     : ?
  Status     : Stable

=cut

sub fetch_by_stable_id {
   my ($self,$stable_id) = @_;

   if(!$stable_id) {
     throw("stable id argument is required");
   }

   my $transcript_adaptor = $self->db()->get_TranscriptAdaptor();
   my $transcript = 
     $transcript_adaptor->fetch_by_translation_stable_id($stable_id);

   return undef if(!$transcript);

   return $self->fetch_by_Transcript($transcript);
}

=head2 fetch_by_stable_id_version

  Arg [1]    : String $id 
               The stable ID of the gene to retrieve
  Arg [2]    : Integer $version
               The version of the stable_id to retrieve
  Example    : $gene = $gene_adaptor->fetch_by_stable_id('ENSG00000148944', 14);
  Description: Retrieves a gene object from the database via its stable id and version.
               The gene will be retrieved in its native coordinate system (i.e.
               in the coordinate system it is stored in the database). It may
               be converted to a different coordinate system through a call to
               transform() or transfer(). If the gene or exon is not found
               undef is returned instead.
  Returntype : Bio::EnsEMBL::Gene or undef
  Exceptions : if we cant get the gene in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_stable_id_version {
   my ($self,$stable_id, $version) = @_;

   if(!$stable_id) {
     throw("stable id argument is required");
   }

   # Enforce that version be numeric
   return unless($version =~ /^\d+$/);

   my $transcript_adaptor = $self->db()->get_TranscriptAdaptor();
   my $transcript = 
       $transcript_adaptor->fetch_by_translation_stable_id_version($stable_id, $version);

   return if(!$transcript);

   return $self->fetch_by_Transcript($transcript);
}

=head2 fetch_all_by_Transcript_list

  Arg [1]    : reference to list of Bio::EnsEMBL::Transcripts $transcripts
               The list of $transcripts to obtain Translation object for.
  Example    : @translations = @{$tla->fetch_all_by_Transcript_list([$t1,$t2]);
  Description: Fetches all translations associated with the list of transcripts
               passed to this method.  The passed transcripts will also have
               their translation set by this method.
  Returntype : Reference to list of Bio::EnsEMBL::Translations
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Transcript_list {
  my ($self,$transcripts) = @_;

  if(!defined($transcripts) || ref($transcripts) ne 'ARRAY') {
    throw("reference to list of Transcripts argument is required");
  }

  return [] if(!@$transcripts);

  my %trans_hash = map {$_->dbID() => $_} @$transcripts;
  my @id_list = keys %trans_hash;

  my @out;

  # mysql is faster and we ensure that we do not exceed the max query size by
  # splitting large queries into smaller queries of 200 ids
  my $max_size = 200;

  my ( $transcript_id, $translation_id, $start_exon_id, $end_exon_id,
       $seq_start, $seq_end, $stable_id, $version, 
       $created_date, $modified_date );

  my %ex_hash;

  while(@id_list) {
    my @ids;
    if(@id_list > $max_size) {
      @ids = splice(@id_list, 0, $max_size);
    } else {
      @ids = splice(@id_list, 0);
    }

    my $id_str;
    if(@ids > 1)  {
      $id_str = " IN (" . join(',', @ids). ")";
    } else {
      $id_str = " = " . $ids[0];
    }
    
    my $canonical_lookup = $self->dbc()->sql_helper()->execute_into_hash(
      -SQL => 'SELECT transcript_id, canonical_translation_id FROM transcript WHERE transcript_id '.$id_str
    );

    my $created_date = $self->db->dbc->from_date_to_seconds("tl.created_date");
    my $modified_date = $self->db->dbc->from_date_to_seconds("tl.modified_date");

    my $sth = $self->prepare
      ("SELECT tl.transcript_id, tl.translation_id, tl.start_exon_id,
           tl.end_exon_id, tl.seq_start, tl.seq_end,
           tl.stable_id, tl.version, " . $created_date . "," .
       $modified_date . 
       " FROM translation tl
         WHERE tl.transcript_id $id_str");

    $sth->execute();

    $sth->bind_columns( \$transcript_id, \$translation_id, \$start_exon_id, \$end_exon_id,
                        \$seq_start, \$seq_end, \$stable_id, \$version,
			\$created_date, \$modified_date );

    while($sth->fetch()) {
      my ($start_exon, $end_exon);

      # this will load all the exons whenever we load the translation
      # but I guess thats ok ....

      my $tr = $trans_hash{$transcript_id};

      foreach my $exon (@{$tr->get_all_Exons()}) {
        if(!$start_exon && $exon->dbID() == $start_exon_id ) {
          $start_exon = $exon;
          last if($end_exon);
        }

        if(!$end_exon && $exon->dbID() == $end_exon_id ) {
          $end_exon = $exon;
          last if($start_exon);
        }
      }

      unless($start_exon && $end_exon) {
        throw("Could not find start or end exon in transcript\n");
      }

      my $tl =  Bio::EnsEMBL::Translation->new
        (-dbID => $translation_id,
         -seq_start => $seq_start,
         -seq_end => $seq_end,
         -start_exon => $start_exon,
         -end_exon => $end_exon,
         -stable_id => $stable_id,
         -version => $version,
	 -created_date => $created_date || undef,
	 -modified_date => $modified_date || undef);
      
      $tl->adaptor($self);
      my $canonical_translation_id = $canonical_lookup->{$transcript_id};
      $tr->translation($tl) if $translation_id == $canonical_translation_id;

      push @out, $tl;
    }
  }

  return \@out;
}


=head2 fetch_all

  Example     : $translations = $translation_adaptor->fetch_all();
  Description : Retrieves all canonical and alternative translations 
                stored in the database.
  Returntype  : listref of Bio::EnsEMBL::Translation
  Caller      : general
  Status      : At Risk

=cut

sub fetch_all {
  my ($self) = @_;
  my $transcript_adaptor = $self->db()->get_TranscriptAdaptor();

  my @translations;
  foreach my $transcript (@{$transcript_adaptor->fetch_all}) {
     my $translation = $self->fetch_by_Transcript($transcript);
     if ($translation) {
	 push @translations, $translation;
     }
     foreach my $alt_translation (@{$self->fetch_all_alternative_by_Transcript($transcript)}) {
	 push @translations, $alt_translation;
     }
  }
  return \@translations;
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
  return (['translation', 'tl']);
}

1;
