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

Bio::EnsEMBL::DBSQL::DBEntryAdaptor -
MySQL Database queries to load and store external object references.

=head1 SYNOPSIS

  $db_entry_adaptor =
    $registry->get_adaptor( 'Human', 'Core', 'DBEntry' );

  $db_entry = $db_entry_adaptor->fetch_by_dbID($id);

  my $gene_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Gene' );

  my $gene = $gene_adaptor->fetch_by_stable_id('ENSG00000101367');

  @db_entries = @{ $db_entry_adaptor->fetch_all_by_Gene($gene) };
  @gene_ids   = $db_entry_adaptor->list_gene_ids_by_extids('BAB15482');

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::DBEntryAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::IdentityXref;
use Bio::EnsEMBL::GoXref;

use Bio::EnsEMBL::Utils::Exception qw(deprecate throw warning);

use vars qw(@ISA);
use strict;

@ISA = qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );

=head2 fetch_by_dbID

  Arg [1]    : int $dbID
               the unique database identifier for the DBEntry to retrieve
  Example    : my $db_entry = $db_entry_adaptor->fetch_by_dbID($dbID);
  Description: Retrieves a dbEntry from the database via its unique
               identifier.
  Returntype : Bio::EnsEMBL::DBEntry
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_dbID {
  my ( $self, $dbID ) = @_;

  my $sth = $self->prepare(
    "SELECT  xref.xref_id,
            xref.dbprimary_acc,
            xref.display_label,
            xref.version,
            exDB.dbprimary_acc_linkable,
            exDB.display_label_linkable,
            exDB.priority,
            exDB.db_name,
            exDB.db_display_name,
            exDB.db_release,
            es.synonym,
            xref.info_type,
            xref.info_text,
            exDB.type,
            exDB.secondary_db_name,
            exDB.secondary_db_table,
            xref.description
    FROM    (xref, external_db exDB)
    LEFT JOIN external_synonym es ON
            es.xref_id = xref.xref_id
    WHERE   xref.xref_id = ?
    AND     xref.external_db_id = exDB.external_db_id" );

  $sth->bind_param( 1, $dbID, SQL_INTEGER );
  $sth->execute();

  my $exDB;

  my $max_rows = 1000;

  while ( my $rowcache = $sth->fetchall_arrayref( undef, $max_rows ) ) {
      #$description refers to the external_db description, while $desc was referring the xref description
    while ( my $arrayref = shift( @{$rowcache} ) ) {
      my ( $refID,               $dbprimaryId,
           $displayid,           $version,
	   $primary_id_linkable,
           $display_id_linkable, $priority,
           $dbname,              $db_display_name,
           $release,             $synonym,
           $info_type,           $info_text,
           $type,                $secondary_db_name,
           $secondary_db_table,  $description
      ) = @$arrayref;

      if ( !defined($exDB) ) {
        $exDB =
          Bio::EnsEMBL::DBEntry->new(
                           -adaptor             => $self,
                           -dbID                => $dbID,
                           -primary_id          => $dbprimaryId,
                           -display_id          => $displayid,
                           -version             => $version,
                           -release             => $release,
                           -dbname              => $dbname,
                           -primary_id_linkable => $primary_id_linkable,
                           -display_id_linkable => $display_id_linkable,
                           -priority            => $priority,
                           -db_display_name     => $db_display_name,
                           -info_type           => $info_type,
                           -info_text           => $info_text,
                           -type                => $type,
                           -secondary_db_name   => $secondary_db_name,
                           -secondary_db_table  => $secondary_db_table,
			   -description         => $description
          );


      }

      if ( defined($synonym) ) { $exDB->add_synonym($synonym) }

    } ## end while ( my $arrayref = shift...
  } ## end while ( my $rowcache = $sth...

  $sth->finish();

  return $exDB;
} ## end sub fetch_by_dbID


=head2 fetch_by_db_accession

  Arg [1]    : string $dbname - The name of the database which the provided
               accession is for.
  Arg [2]    : string $accession - The accesion of the external reference to
               retrieve.
  Example    : my $xref = $dbea->fetch_by_db_accession('Interpro','IPR003439');
               print $xref->description(), "\n" if($xref);
  Description: Retrieves a DBEntry (xref) via the name of the database
               it is from and its primary accession in that database.
               Undef is returned if the xref cannot be found in the
               database.
               NOTE:  In a multi-species database, this method will
               return all the entries matching the search criteria, not
               just the ones associated with the current species.
  Returntype : Bio::EnsEMBL::DBSQL::DBEntry
  Exceptions : thrown if arguments are incorrect
  Caller     : general, domainview
  Status     : Stable

=cut

sub fetch_by_db_accession {
  my ( $self, $dbname, $accession ) = @_;

  my $sth = $self->prepare(
    "SELECT xref.xref_id,
            xref.dbprimary_acc,
            xref.display_label,
            xref.version,
            exDB.dbprimary_acc_linkable,
            exDB.display_label_linkable,
            exDB.priority,
            exDB.db_name,
            exDB.db_display_name,
            exDB.db_release,
            es.synonym,
            xref.info_type,
            xref.info_text,
            exDB.type,
            exDB.secondary_db_name,
            exDB.secondary_db_table,
            xref.description
    FROM    (xref, external_db exDB)
    LEFT JOIN external_synonym es ON
            es.xref_id = xref.xref_id
    WHERE  xref.dbprimary_acc = ?
    AND    exDB.db_name = ?
    AND    xref.external_db_id = exDB.external_db_id" );

  $sth->bind_param( 1, $accession, SQL_VARCHAR );
  $sth->bind_param( 2, $dbname,    SQL_VARCHAR );
  $sth->execute();

  if ( !$sth->rows() && lc($dbname) eq 'interpro' ) {
  # This is a minor hack that means that results still come back even
  # when a mistake was made and no interpro accessions were loaded into
  # the xref table.  This has happened in the past and had the result of
  # breaking domainview

    $sth->finish();
    $sth = $self->prepare(
      "SELECT   NULL,
                i.interpro_ac,
                i.id,
                NULL,
                NULL,
                'Interpro',
                NULL,
                NULL
        FROM    interpro i
        WHERE   i.interpro_ac = ?" );

    $sth->bind_param( 1, $accession, SQL_VARCHAR );
    $sth->execute();
  }

  my $exDB;

  my $max_rows = 1000;

  while ( my $rowcache = $sth->fetchall_arrayref( undef, $max_rows ) ) {
    while ( my $arrayref = shift( @{$rowcache} ) ) {
      my ( $dbID,                $dbprimaryId,
           $displayid,           $version,
           $primary_id_linkable,
           $display_id_linkable, $priority,
           $dbname,              $db_display_name,
           $release,             $synonym,
           $info_type,           $info_text,
           $type,                $secondary_db_name,
           $secondary_db_table,  $description
      ) = @$arrayref;

      if ( !defined($exDB) ) {
        $exDB =
          Bio::EnsEMBL::DBEntry->new(
                           -adaptor             => $self,
                           -dbID                => $dbID,
                           -primary_id          => $dbprimaryId,
                           -display_id          => $displayid,
                           -version             => $version,
                           -release             => $release,
                           -dbname              => $dbname,
                           -primary_id_linkable => $primary_id_linkable,
                           -display_id_linkable => $display_id_linkable,
                           -priority            => $priority,
                           -db_display_name     => $db_display_name,
                           -info_type           => $info_type,
                           -info_text           => $info_text,
                           -type                => $type,
                           -secondary_db_name   => $secondary_db_name,
                           -secondary_db_table  => $secondary_db_table,
			   -description         => $description
          );


      }

      if ($synonym) { $exDB->add_synonym($synonym) }

    } ## end while ( my $arrayref = shift...
  } ## end while ( my $rowcache = $sth...

  $sth->finish();

  return $exDB;
} ## end sub fetch_by_db_accession


=head2 store

  Arg [1]    : Bio::EnsEMBL::DBEntry $exObj
               The DBEntry (xref) to be stored
  Arg [2]    : Int $ensID
               The dbID of an EnsEMBL object to associate with this external
               database entry
  Arg [3]    : string $ensType ('Transcript', 'Translation', 'Gene')
               The type of EnsEMBL object that this external database entry is
               being associated with.
  Arg [4]    : boolean $ignore_release
               If unset or zero, will require that the release string
               of the DBEntry object is identical to the release of the
               external database.  If set and non-zero, will ignore the
               release information.
  Example    : $dbea->store($db_entry, $transcript_id, 'Transcript');
  Description: Stores a reference to an external database (if it is not stored
               already) and associates an EnsEMBL object of a specified type
               with the external identifier.
  Returntype : int - the dbID of the newly created external refernce
  Exceptions : thrown when invalid dbID is passed to this method
  Caller     : scripts which load Xrefs and ObjectXrefs, etc. into Ensembl
  Status     : Stable

=cut

sub store {
  my ( $self, $exObj, $ensID, $ensType, $ignore_release ) = @_;

  my $dbJustInserted;

  #
  # backwards compatibility check:
  # check if $ensID is an object; if so, use $obj->dbID
  #
  my $ensembl_id;

  if ($ensID =~ /^\d+$/) {
    $ensembl_id = $ensID;

  } elsif ( ref($ensID) eq 'Bio::EnsEMBL::Gene' or
            ref($ensID) eq 'Bio::EnsEMBL::Transcript' or
            ref($ensID) eq 'Bio::EnsEMBL::Translation') {

    warning("You should pass DBEntryAdaptor->store() a dbID rather than an ensembl object to store the xref on");

    if ($ensID->dbID) {
      $ensembl_id = $ensID->dbID;
    } else {
      throw("$ensType ".$ensID->display_id." doesn't have a dbID, can't store xref");
    }

  } else {
    throw("Invalid dbID passed to DBEntryAdaptor->store()");
  }

  #
  # Check for the existance of the external_db, throw if it does not exist
  #

  my ($dbRef, $release_clause);

  if ( !$ignore_release ) {

	if(defined $exObj->release()){
	  $release_clause = " AND db_release = ?";
	}
	else{
	  $release_clause = " AND db_release is NULL";
	}


    my $sth = $self->prepare( "
     SELECT external_db_id
       FROM external_db
      WHERE db_name    = ?
        $release_clause" );

    $sth->bind_param( 1, $exObj->dbname(),  SQL_VARCHAR );
    $sth->bind_param( 2, $exObj->release(), SQL_VARCHAR ) if defined $exObj->release();

    $sth->execute();

    ($dbRef) = $sth->fetchrow_array();

    if ( !$dbRef ) {
      throw(
             sprintf( "external_db [%s] release [%s] does not exist",
                      $exObj->dbname(), $exObj->release() ) );
    }
  } else {

    my $sth = $self->prepare( "
     SELECT external_db_id
       FROM external_db
      WHERE db_name = ? " );

    $sth->bind_param( 1, $exObj->dbname(), SQL_VARCHAR );

    $sth->execute();

    ($dbRef) = $sth->fetchrow_array();

    if ( !$dbRef ) {
      throw(
          sprintf( "external_db [%s] does not exist", $exObj->dbname() )
      );
    }
  }

  #
  # Check for the existance of the external reference, add it if not present
  #
  my $sth = $self->prepare( "
       SELECT xref_id
         FROM xref
        WHERE external_db_id = ?
          AND dbprimary_acc  = ?
          AND version        = ?
          AND info_type      = ?
          AND info_text      = ?" );

  $sth->bind_param(1,$dbRef,SQL_INTEGER);
  $sth->bind_param(2,$exObj->primary_id,SQL_VARCHAR);
  $sth->bind_param(3,$exObj->version,SQL_VARCHAR);
  $sth->bind_param(4,$exObj->info_type,SQL_VARCHAR);
  $sth->bind_param(5,$exObj->info_text,SQL_VARCHAR);

  $sth->execute();
  my ($dbX) = $sth->fetchrow_array();

  $sth->finish();
  if(!$dbX) {
    if(!$exObj->primary_id()) {
      throw("DBEntry cannot be stored without a primary_id attribute.");
    }

    #
    # store the new xref
    #
    $sth = $self->prepare( "
       INSERT ignore INTO xref
       SET dbprimary_acc = ?,
           display_label = ?,
           version = ?,
           description = ?,
           external_db_id = ?,
           info_type = ?,
           info_text = ?");
    $sth->bind_param(1, $exObj->primary_id,SQL_VARCHAR);
    $sth->bind_param(2, $exObj->display_id,SQL_VARCHAR);
    $sth->bind_param(3, $exObj->version,SQL_VARCHAR);
    $sth->bind_param(4, $exObj->description,SQL_VARCHAR);
    $sth->bind_param(5, $dbRef,SQL_INTEGER);
    $sth->bind_param(6, $exObj->info_type, SQL_VARCHAR);
    $sth->bind_param(7, $exObj->info_text, SQL_VARCHAR);

    $sth->execute();

    $dbX = $sth->{'mysql_insertid'};
    $sth->finish();
    #
    # store the synonyms for the new xref
    # 
    my $synonym_check_sth = $self->prepare(
              "SELECT xref_id, synonym
               FROM external_synonym
               WHERE xref_id = ?
               AND synonym = ?");

    my $synonym_store_sth = $self->prepare(
        "INSERT ignore INTO external_synonym
         SET xref_id = ?, synonym = ?");

    my $synonyms = $exObj->get_all_synonyms();
    foreach my $syn ( @$synonyms ) {
	$synonym_check_sth->bind_param(1,$dbX,SQL_INTEGER);
	$synonym_check_sth->bind_param(2,$syn,SQL_VARCHAR);
	$synonym_check_sth->execute();
      my ($dbSyn) = $synonym_check_sth->fetchrow_array(); 
	$synonym_store_sth->bind_param(1,$dbX,SQL_INTEGER);
	$synonym_store_sth->bind_param(2,$syn,SQL_VARCHAR);
	$synonym_store_sth->execute() if(!$dbSyn);
    }
    $synonym_check_sth->finish();
    $synonym_store_sth->finish();
  }
  #
  # check if the object mapping was already stored
  #
  $sth = $self->prepare(
    qq(
SELECT  xref_id
FROM    object_xref
WHERE   xref_id = ?
  AND   ensembl_object_type = ?
  AND   ensembl_id = ?
  AND   (   linkage_annotation = ?
    OR      linkage_annotation IS NULL  )) );

  $sth->bind_param( 1, $dbX,                         SQL_INTEGER );
  $sth->bind_param( 2, $ensType,                     SQL_VARCHAR );
  $sth->bind_param( 3, $ensembl_id,                  SQL_INTEGER );
  $sth->bind_param( 4, $exObj->linkage_annotation(), SQL_VARCHAR );

  $sth->execute();

  my ($tst) = $sth->fetchrow_array();

  $sth->finish();

  if(!$tst) {
    #
    # Store the reference to the internal ensembl object
    #
    my $analysis_id;
    if($exObj->analysis()) {
      $analysis_id =
	$self->db()->get_AnalysisAdaptor->store($exObj->analysis());
    } else {
      $analysis_id = undef;
    }

    $sth = $self->prepare(
      qq(
INSERT IGNORE INTO object_xref
  SET   xref_id = ?,
        ensembl_object_type = ?,
        ensembl_id = ?,
        linkage_annotation = ?,
        analysis_id = ? ) );

    $sth->bind_param( 1, $dbX,                         SQL_INTEGER );
    $sth->bind_param( 2, $ensType,                     SQL_VARCHAR );
    $sth->bind_param( 3, $ensembl_id,                  SQL_INTEGER );
    $sth->bind_param( 4, $exObj->linkage_annotation(), SQL_VARCHAR );
    $sth->bind_param( 5, $analysis_id,                 SQL_INTEGER);
 
   #print "stored xref id $dbX in obejct_xref\n";
    $sth->execute();
    $exObj->dbID($dbX);
    $exObj->adaptor($self);
    my $Xidt = $sth->{'mysql_insertid'};

    #
    # If this is an IdentityXref need to store in that table too
    # If its GoXref add the linkage type to go_xref table
    #
    if ($exObj->isa('Bio::EnsEMBL::IdentityXref')) {
      $sth = $self->prepare( "
             INSERT ignore INTO identity_xref
             SET object_xref_id = ?,
             xref_identity = ?,
             ensembl_identity = ?,
             xref_start = ?,
             xref_end   = ?,
             ensembl_start = ?,
             ensembl_end = ?,
             cigar_line = ?,
             score = ?,
             evalue = ?" );
      $sth->bind_param(1,$Xidt,SQL_INTEGER);
      $sth->bind_param(2,$exObj->xref_identity,SQL_INTEGER);
      $sth->bind_param(3,$exObj->ensembl_identity,SQL_INTEGER);
      $sth->bind_param(4,$exObj->xref_start,SQL_INTEGER);
      $sth->bind_param(5,$exObj->xref_end,SQL_INTEGER);
      $sth->bind_param(6,$exObj->ensembl_start,SQL_INTEGER);
      $sth->bind_param(7,$exObj->ensembl_end,SQL_INTEGER);
      $sth->bind_param(8,$exObj->cigar_line,SQL_LONGVARCHAR);
      $sth->bind_param(9,$exObj->score,SQL_DOUBLE);
      $sth->bind_param(10,$exObj->evalue,SQL_DOUBLE);
      $sth->execute();
    } elsif( $exObj->isa( 'Bio::EnsEMBL::GoXref' )) {
      $sth = $self->prepare( "
             INSERT ignore INTO go_xref
                SET object_xref_id = ?,
                    source_xref_id = ?,
                    linkage_type = ? " );
      foreach my $info (@{$exObj->get_all_linkage_info()}) {
	  my( $lt, $sourceXref ) = @{$info};
	  my $sourceXid = undef;
	  if( $sourceXref ){ 
	      $sourceXref->is_stored($self->dbc) || $sourceXref->store;
	      $sourceXid = $sourceXref->dbID;
	  }
	  $sth->bind_param(1,$Xidt,SQL_INTEGER);
	  $sth->bind_param(2,$sourceXid,SQL_INTEGER);
	  $sth->bind_param(3,$lt,SQL_VARCHAR);
	  $sth->execute();
      }
    }
  } 
  return $dbX;
}


=head2 exists

  Arg [1]    : Bio::EnsEMBL::DBEntry $dbe
  Example    : if($dbID = $db_entry_adaptor->exists($dbe)) { do stuff; }
  Description: Returns the db id of this DBEntry if it exists in this database
               otherwise returns undef.  Exists is defined as an entry with 
               the same external_db and display_id
  Returntype : int
  Exceptions : thrown on incorrect args
  Caller     : GeneAdaptor::store, TranscriptAdaptor::store
  Status     : Stable

=cut

sub exists {
  my ($self, $dbe) = @_ ;

  unless($dbe && ref $dbe && $dbe->isa('Bio::EnsEMBL::DBEntry')) {
    throw("arg must be a Bio::EnsEMBL::DBEntry not [$dbe]");
  }

  my $sth = $self->prepare('SELECT x.xref_id 
                            FROM   xref x, external_db xdb
                            WHERE  x.external_db_id = xdb.external_db_id
                            AND    x.display_label = ?
                            AND    xdb.db_name = ?
                            AND    x.dbprimary_acc = ?');

  $sth->bind_param(1,$dbe->display_id,SQL_VARCHAR);
  $sth->bind_param(2,$dbe->dbname,SQL_VARCHAR);
  $sth->bind_param(3,$dbe->primary_id,SQL_VARCHAR);
  $sth->execute();

  my ($dbID) = $sth->fetchrow_array;

  $sth->finish;

  return $dbID;
}


=head2 fetch_all_by_Gene

  Arg [1]    : Bio::EnsEMBL::Gene $gene 
               (The gene to retrieve DBEntries for)
  Arg [2]    : optional external database name
  Arg [3]    : optional external_db type 
  Example    : @db_entries = @{$db_entry_adaptor->fetch_by_Gene($gene)};
  Description: This returns a list of DBEntries associated with this gene.
               Note that this method was changed in release 15.  Previously
               it set the DBLinks attribute of the gene passed in to contain
               all of the gene, transcript, and translation xrefs associated
               with this gene.
  Returntype : listref of Bio::EnsEMBL::DBEntries; may be of type IdentityXref if
               there is mapping data, or GoXref if there is linkage data.
  Exceptions : thows if gene object not passed
  Caller     : Bio::EnsEMBL::Gene
  Status     : Stable

=cut

sub fetch_all_by_Gene {
  my ( $self, $gene, $ex_db_reg, $exdb_type ) = @_;

  if(!ref($gene) || !$gene->isa('Bio::EnsEMBL::Gene')) {
    throw("Bio::EnsEMBL::Gene argument expected.");
  }

  return $self->_fetch_by_object_type($gene->dbID(), 'Gene', $ex_db_reg, $exdb_type);
}


=head2 fetch_all_by_Transcript

  Arg [1]    : Bio::EnsEMBL::Transcript
  Arg [2]    : optional external database name
  Arg [3]    : optional external_db type 
  Example    : @db_entries = @{$db_entry_adaptor->fetch_by_Gene($trans)};
  Description: This returns a list of DBEntries associated with this 
               transcript. Note that this method was changed in release 15.
               Previously it set the DBLinks attribute of the gene passed in
               to contain all of the gene, transcript, and translation xrefs
               associated with this gene.
  Returntype : listref of Bio::EnsEMBL::DBEntries; may be of type IdentityXref if
               there is mapping data, or GoXref if there is linkage data.
  Exceptions : throes if transcript argument not passed
  Caller     : Bio::EnsEMBL::Gene
  Status     : Stable

=cut

sub fetch_all_by_Transcript {
  my ( $self, $trans, $ex_db_reg, $exdb_type ) = @_;

  if(!ref($trans) || !$trans->isa('Bio::EnsEMBL::Transcript')) {
    throw("Bio::EnsEMBL::Transcript argument expected.");
  }

  return $self->_fetch_by_object_type( $trans->dbID(), 'Transcript', $ex_db_reg, $exdb_type);
}


=head2 fetch_all_by_Translation

  Arg [1]    : Bio::EnsEMBL::Translation $trans
               (The translation to fetch database entries for)
  Arg [2]    : optional external database name
  Arg [3]    : optional externaldb type 
  Example    : @db_entries = @{$db_entry_adptr->fetch_by_Translation($trans)};
  Description: Retrieves external database entries for an EnsEMBL translation
  Returntype : listref of Bio::EnsEMBL::DBEntries; may be of type IdentityXref if
               there is mapping data, or GoXref if there is linkage data.
  Exceptions : throws if translation object not passed
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Translation {
  my ( $self, $trans, $ex_db_reg, $exdb_type ) = @_;

  if(!ref($trans) || !$trans->isa('Bio::EnsEMBL::Translation')) {
    throw('Bio::EnsEMBL::Translation argument expected.');
  }
  if( ! $trans->dbID ){ 
    warning( "Cannot fetch_all_by_Translation without a dbID" );
    return [];
  }

  return $self->_fetch_by_object_type( $trans->dbID(), 'Translation', $ex_db_reg, $exdb_type );
}



=head2 remove_from_object

  Arg [1]    : Bio::EnsEMBL::DBEntry $dbe - The external reference which
               is to be disassociated from an ensembl object.
  Arg [2]    : Bio::EnsEMBL::Storable $object - The ensembl object the
               external reference is to be disassociated from
  Arg [3]    : string $object_type - The type of the ensembl object.
               E.g. 'Gene', 'Transcript', 'Translation'
  Example    :
               # remove all dbentries from this translation
               foreach my $dbe (@{$translation->get_all_DBEntries()}) {
                 $dbe_adaptor->remove($dbe, $translation, 'Translation');
               }
  Description: Removes an association between an ensembl object and a
               DBEntry (xref).  This does not remove the actual xref from
               the database, only its linkage to the ensembl object.
  Returntype : none
  Exceptions : Throw on incorrect arguments.
               Warning if object or dbentry is not stored in this database.
  Caller     : TranscriptAdaptor::remove, GeneAdaptor::remove,
               TranslationAdaptor::remove
  Status     : Stable

=cut

sub remove_from_object {
  my $self = shift;
  my $dbe  = shift;
  my $object = shift;
  my $object_type = shift;

  if(!ref($dbe) || !$dbe->isa('Bio::EnsEMBL::DBEntry')) {
    throw("Bio::EnsEMBL::DBEntry argument expected.");
  }

  if(!ref($object) || !$dbe->isa('Bio::EnsEMBL::Storable')) {
    throw("Bio::EnsEMBL::Storable argument expected.");
  }

  if(!$object_type) {
    throw("object_type string argument expected.");
  }

  # make sure both the dbentry and the object it is allegedly linked to
  # are stored in this database

  if(!$object->is_stored($self->db())) {
    warning("Cannot remove DBEntries for $object_type " . $object->dbID() .
            ". Object is not stored in this database.");
    return;
  }

  if(!$dbe->is_stored($self->db())) {
    warning("Cannot remove DBEntry ".$dbe->dbID() . ". Is not stored " .
            "in this database.");
    return;
  }

  # obtain the identifier of the link from the object_xref table
  #No need to compare linkage_annotation here
  my $sth = $self->prepare
    ("SELECT ox.object_xref_id " .
     "FROM   object_xref ox ".
     "WHERE  ox.xref_id = ? " .
     "AND    ox.ensembl_id = ? " .
     "AND    ox.ensembl_object_type = ?");
  $sth->bind_param(1,$dbe->dbID,SQL_INTEGER);
  $sth->bind_param(2,$object->dbID,SQL_INTEGER);
  $sth->bind_param(3,$object_type,SQL_VARCHAR);
  $sth->execute();

  if(!$sth->rows() == 1) {
    $sth->finish();
    return;
  }

  my ($ox_id) = $sth->fetchrow_array();
  $sth->finish();

  # delete from the tables which contain additional linkage information

  $sth = $self->prepare("DELETE FROM go_xref WHERE object_xref_id = ?");
  $sth->bind_param(1,$ox_id,SQL_INTEGER);
  $sth->execute();
  $sth->finish();

  $sth = $self->prepare("DELETE FROM identity_xref WHERE object_xref_id = ?");
  $sth->bind_param(1,$ox_id,SQL_INTEGER);
  $sth->execute();
  $sth->finish();

  # delete the actual linkage itself
  $sth = $self->prepare("DELETE FROM object_xref WHERE object_xref_id = ?");
  $sth->bind_param(1,$ox_id,SQL_INTEGER);
  $sth->execute();
  $sth->finish();

  return;
}


=head2 _fetch_by_object_type

  Arg [1]    : string $ensID
  Arg [2]    : string $ensType
  			   (object type to be returned) 
  Arg [3]    : optional $exdbname (external database name)
  Arf [4]    : optional $exdb_type (external database type)
  Example    : $self->_fetch_by_object_type( $translation_id, 'Translation' )
  Description: Fetches DBEntry by Object type
               NOTE:  In a multi-species database, this method will
               return all the entries matching the search criteria, not
               just the ones associated with the current species.
  Returntype : arrayref of DBEntry objects; may be of type IdentityXref if
               there is mapping data, or GoXref if there is linkage data.
  Exceptions : none
  Caller     : fetch_all_by_Gene
               fetch_all_by_Translation
               fetch_all_by_Transcript
  Status     : Stable

=cut

sub _fetch_by_object_type {
  my ( $self, $ensID, $ensType, $exdbname, $exdb_type ) = @_;

  my @out;

  if ( !defined($ensID) ) {
    throw("Can't fetch_by_EnsObject_type without an object");
  }

  if ( !defined($ensType) ) {
    throw("Can't fetch_by_EnsObject_type without a type");
  }

  #  my $sth = $self->prepare("
  my $sql = (<<SSQL);
    SELECT xref.xref_id, xref.dbprimary_acc, xref.display_label, xref.version,
           exDB.dbprimary_acc_linkable, exDB.display_label_linkable,
           exDB.priority,
           exDB.db_name, exDB.db_release, exDB.status, exDB.db_display_name,
           exDB.secondary_db_name, exDB.secondary_db_table,
           oxr.object_xref_id,
           es.synonym,
           idt.xref_identity, idt.ensembl_identity, idt.xref_start,
           idt.xref_end, idt.ensembl_start, idt.ensembl_end,
           idt.cigar_line, idt.score, idt.evalue, oxr.analysis_id,
           gx.linkage_type,
           xref.info_type, xref.info_text, exDB.type, gx.source_xref_id,
           oxr.linkage_annotation, xref.description
    FROM   (xref xref, external_db exDB, object_xref oxr)
    LEFT JOIN external_synonym es on es.xref_id = xref.xref_id 
    LEFT JOIN identity_xref idt on idt.object_xref_id = oxr.object_xref_id
    LEFT JOIN go_xref gx on gx.object_xref_id = oxr.object_xref_id
    WHERE  xref.xref_id = oxr.xref_id
      AND  xref.external_db_id = exDB.external_db_id 
      AND  oxr.ensembl_id = ?
      AND  oxr.ensembl_object_type = ?
SSQL
  $sql .= " AND exDB.db_name like '" . $exdbname . "' " if ($exdbname);
  $sql .= " AND exDB.type like '" . $exdb_type . "' "   if ($exdb_type);
  my $sth = $self->prepare($sql);

  $sth->bind_param( 1, $ensID,   SQL_INTEGER );
  $sth->bind_param( 2, $ensType, SQL_VARCHAR );
  $sth->execute();

  my ( %seen, %linkage_types, %synonyms );

  my $max_rows = 1000;

  while ( my $rowcache = $sth->fetchall_arrayref( undef, $max_rows ) ) {
    while ( my $arrRef = shift( @{$rowcache} ) ) {
      my ( $refID,                  $dbprimaryId,
           $displayid,              $version,
           $primary_id_linkable,
           $display_id_linkable,    $priority,
           $dbname,                 $release,
           $exDB_status,            $exDB_db_display_name,
           $exDB_secondary_db_name, $exDB_secondary_db_table,
           $objid,                  $synonym,
           $xrefid,                 $ensemblid,
           $xref_start,             $xref_end,
           $ensembl_start,          $ensembl_end,
           $cigar_line,             $score,
           $evalue,                 $analysis_id,
           $linkage_type,           $info_type,
           $info_text,              $type,
           $source_xref_id,         $link_annotation,
	   $description
      ) = @$arrRef;

      my $linkage_key =
        ( $linkage_type || '' ) . ( $source_xref_id || '' );


      my $analysis = undef;
      if ( defined($analysis_id) ) {
	$analysis =
	  $self->db()->get_AnalysisAdaptor()->fetch_by_dbID($analysis_id);
      }

      my %obj_hash = ( 'adaptor'            => $self,
                       'dbID'               => $refID,
                       'primary_id'         => $dbprimaryId,
                       'display_id'         => $displayid,
                       'version'            => $version,
                       'release'            => $release,
                       'info_type'          => $info_type,
                       'info_text'          => $info_text,
                       'type'               => $type,
                       'secondary_db_name'  => $exDB_secondary_db_name,
                       'secondary_db_table' => $exDB_secondary_db_table,
                       'dbname'             => $dbname,
                       'description'        => $description,
                       'linkage_annotation' => $link_annotation,
                       'analysis'           => $analysis);

      # Using an outer join on the synonyms as well as on identity_xref,
      # we now have to filter out the duplicates (see v.1.18 for
      # original). Since there is at most one identity_xref row per
      # xref, this is easy enough; all the 'extra' bits are synonyms.
      my $source_xref;
      if ( !$seen{$refID} ) {
	
	my $exDB;
        if ( ( defined($xrefid) ) ) {  # an xref with similarity scores
          $exDB = Bio::EnsEMBL::IdentityXref->new_fast( \%obj_hash );
          $exDB->xref_identity($xrefid);
          $exDB->ensembl_identity($ensemblid);

          $exDB->cigar_line($cigar_line);
          $exDB->xref_start($xref_start);
          $exDB->xref_end($xref_end); # was not here before 14th Jan 2009 ????
          $exDB->ensembl_start($ensembl_start);
          $exDB->ensembl_end($ensembl_end);
          $exDB->score($score);
          $exDB->evalue($evalue);

        } elsif ( defined $linkage_type && $linkage_type ne "" ) {
          $exDB = Bio::EnsEMBL::GoXref->new_fast( \%obj_hash );
          $source_xref = ( defined($source_xref_id)
                              ? $self->fetch_by_dbID($source_xref_id)
                              : undef );
          $exDB->add_linkage_type( $linkage_type, $source_xref || () );
          $linkage_types{$refID}->{$linkage_key} = 1;

        } else {
          $exDB = Bio::EnsEMBL::DBEntry->new_fast( \%obj_hash );
        }

        if ( defined($exDB_status) ) { $exDB->status($exDB_status) }

        $exDB->primary_id_linkable($primary_id_linkable);
        $exDB->display_id_linkable($display_id_linkable);
        $exDB->priority($priority);
        $exDB->db_display_name($exDB_db_display_name);

        push( @out, $exDB );
        $seen{$refID} = $exDB;

      } ## end if ( !$seen{$refID} )

      # $exDB still points to the same xref, so we can keep adding GO
      # evidence tags or synonyms.

      if ( defined($synonym) && !$synonyms{$refID}->{$synonym} ) {
        if ( defined($synonym) ) {
          $seen{$refID}->add_synonym($synonym);
        }
        $synonyms{$refID}->{$synonym} = 1;
      }

      if (    defined($linkage_type)
           && $linkage_type ne ""
           && !$linkage_types{$refID}->{$linkage_key} )
      {
        $source_xref = ( defined($source_xref_id)
                            ? $self->fetch_by_dbID($source_xref_id)
                            : undef );
        $seen{$refID}
          ->add_linkage_type( $linkage_type, $source_xref || () );
        $linkage_types{$refID}->{$linkage_key} = 1;
      }
    } ## end while ( my $arrRef = shift...
  } ## end while ( my $rowcache = $sth...

  return \@out;
} ## end sub _fetch_by_object_type

=head2 list_gene_ids_by_external_db_id

  Arg [1]    : string $external_id
  Example    : @gene_ids = $dbea->list_gene_ids_by_external_db_id(1020);
  Description: Retrieve a list of geneid by an external identifier that
               is linked to any of the genes transcripts, translations
               or the gene itself.
               NOTE:  If more than one external identifier has the
               same primary accession then genes for each of these is
               returned.
               NOTE:  In a multi-species database, this method will
               return all the entries matching the search criteria, not
               just the ones associated with the current species.
  Returntype : list of ints
  Exceptions : none
  Caller     : unknown
  Status     : Stable

=cut

sub list_gene_ids_by_external_db_id{
   my ($self,$external_db_id) = @_;

   my %T = map { ($_, 1) }
       $self->_type_by_external_db_id( $external_db_id, 'Translation', 'gene' ),
       $self->_type_by_external_db_id( $external_db_id, 'Transcript',  'gene' ),
       $self->_type_by_external_db_id( $external_db_id, 'Gene' );
   return keys %T;
}

=head2 list_gene_ids_by_extids

  Arg [1]    : string $external_name
  Arg [2]    : (optional) string $external_db_name
  Example    : @gene_ids = $dbea->list_gene_ids_by_extids('CDPX');
  Description: Retrieve a list of geneid by an external identifier that is 
               linked to  any of the genes transcripts, translations or the 
               gene itself 
  Returntype : list of ints
  Exceptions : none
  Caller     : unknown
  Status     : Stable

=cut

sub list_gene_ids_by_extids {
  my ( $self, $external_name, $external_db_name ) = @_;

  my %T = map { ( $_, 1 ) }
    $self->_type_by_external_id( $external_name, 'Translation', 'gene',
                                 $external_db_name ),
    $self->_type_by_external_id( $external_name, 'Transcript', 'gene',
                                 $external_db_name ),
    $self->_type_by_external_id( $external_name, 'Gene', undef,
                                 $external_db_name );

  return keys %T;
}


=head2 list_transcript_ids_by_extids

  Arg [1]    : string $external_name
  Arg [2]    : (optional) string $external_db_name
  Example    : @tr_ids = $dbea->list_transcript_ids_by_extids('BCRA2');
  Description: Retrieve a list transcript ids by an external identifier that 
               is linked to any of the genes transcripts, translations or the 
               gene itself 
  Returntype : list of ints
  Exceptions : none
  Caller     : unknown
  Status     : Stable

=cut

sub list_transcript_ids_by_extids {
  my ( $self, $external_name, $external_db_name ) = @_;

  my %T = map { ( $_, 1 ) }
    $self->_type_by_external_id( $external_name, 'Translation',
                                 'transcript',   $external_db_name
    ),
    $self->_type_by_external_id( $external_name, 'Transcript', undef,
                                 $external_db_name );

  return keys %T;
}


=head2 list_translation_ids_by_extids

  Arg [1]    : string $external_name
  Arg [2]    : (optional) string $external_db_name
  Example    : @tr_ids = $dbea->list_translation_ids_by_extids('GO:0004835');
  Description: Gets a list of translation IDs by external display IDs
  Returntype : list of Ints
  Exceptions : none
  Caller     : unknown
  Status     : Stable

=cut

sub list_translation_ids_by_extids {
  my ( $self, $external_name, $external_db_name ) = @_;

  return
    $self->_type_by_external_id( $external_name, 'Translation', undef,
                                 $external_db_name );
}

=head2 _type_by_external_id

  Arg [1]    : string $name - dbprimary_acc
  Arg [2]    : string $ensType - ensembl_object_type
  Arg [3]    : (optional) string $extraType
  Arg [4]    : (optional) string $external_db_name
  	       other object type to be returned
  Example    : $self->_type_by_external_id($name, 'Translation');
               NOTE:  In a multi-species database, this method will
               return all the entries matching the search criteria, not
               just the ones associated with the current species.
  Description: Gets
  Returntype : list of dbIDs (gene_id, transcript_id, etc.)
  Exceptions : none
  Caller     : list_translation_ids_by_extids
               translationids_by_extids
  	       geneids_by_extids
  Status     : Stable

=cut

sub _type_by_external_id {
  my ( $self, $name, $ensType, $extraType, $external_db_name ) = @_;

  my $from_sql  = '';
  my $where_sql = '';
  my $ID_sql    = 'oxr.ensembl_id';

  if ( defined($extraType) ) {
    if ( lc($extraType) eq 'translation' ) {
      $ID_sql = 'tl.translation_id';
    } else {
      $ID_sql = "t.${extraType}_id";
    }

    if ( lc($ensType) eq 'translation' ) {
      $from_sql  = 'transcript t, translation tl, ';
      $where_sql = qq(
          t.transcript_id = tl.transcript_id AND
          tl.translation_id = oxr.ensembl_id AND
          t.is_current = 1 AND
      );
    } else {
      $from_sql  = 'transcript t, ';
      $where_sql = 't.'
        . lc($ensType)
        . '_id = oxr.ensembl_id AND '
        . 't.is_current = 1 AND ';
    }
  }

  if ( lc($ensType) eq 'gene' ) {
    $from_sql  = 'gene g, ';
    $where_sql = 'g.gene_id = oxr.ensembl_id AND g.is_current = 1 AND ';
  } elsif ( lc($ensType) eq 'transcript' ) {
    $from_sql = 'transcript t, ';
    $where_sql =
      't.transcript_id = oxr.ensembl_id AND t.is_current = 1 AND ';
  } elsif ( lc($ensType) eq 'translation' ) {
    $from_sql  = 'transcript t, translation tl, ';
    $where_sql = qq(
        t.transcript_id = tl.transcript_id AND
        tl.translation_id = oxr.ensembl_id AND
        t.is_current = 1 AND
    );
  }

  if ( defined($external_db_name) ) {
    # Involve the 'external_db' table to limit the hits to a particular
    # external database.

    $from_sql .= 'external_db xdb, ';
    $where_sql .=
        'xdb.db_name LIKE '
      . $self->dbc()->db_handle()->quote( $external_db_name . '%' )
      . ' AND xdb.external_db_id = x.external_db_id AND';
  }

  my $query1 = qq(
      SELECT    $ID_sql
      FROM      $from_sql
                xref x,
                object_xref oxr
      WHERE     $where_sql
                ( x.dbprimary_acc = ? OR x.display_label = ? )
      AND         x.xref_id = oxr.xref_id
      AND         oxr.ensembl_object_type = ?);

  my $query2;

  if ( defined($external_db_name) ) {
    # If we are given the name of an external database, we need to join
    # between the 'xref' and the 'object_xref' tables on 'xref_id'.

    $query2 = qq(
      SELECT    $ID_sql
      FROM      $from_sql
                external_synonym syn,
                object_xref oxr,
                xref x
      WHERE     $where_sql
                syn.synonym = ?
      AND       syn.xref_id = oxr.xref_id
      AND       oxr.ensembl_object_type = ?
      AND       x.xref_id = oxr.xref_id);

  } else {
    # If we weren't given an external database name, we can get away
    # with less joins here.

    $query2 = qq(
      SELECT    $ID_sql
      FROM      $from_sql
                external_synonym syn,
                object_xref oxr
      WHERE     $where_sql
                syn.synonym = ?
      AND       syn.xref_id = oxr.xref_id
      AND       oxr.ensembl_object_type = ?);

  }

  my %result;

  my $sth = $self->prepare($query1);

  $sth->bind_param( 1, "$name",  SQL_VARCHAR );
  $sth->bind_param( 2, "$name",  SQL_VARCHAR );
  $sth->bind_param( 3, $ensType, SQL_VARCHAR );
  $sth->execute();
  my $r;
  while ( $r = $sth->fetchrow_array() ) { $result{$r} = 1 }

  $sth = $self->prepare($query2);

  $sth->bind_param( 1, "$name",  SQL_VARCHAR );
  $sth->bind_param( 2, $ensType, SQL_VARCHAR );
  $sth->execute();

  while ( $r = $sth->fetchrow_array() ) { $result{$r} = 1 }

  return keys(%result);

} ## end sub _type_by_external_id

=head2 _type_by_external_db_id

  Arg [1]    : string $type - external_db type
  Arg [2]    : string $ensType - ensembl_object_type
  Arg [3]    : (optional) string $extraType
  	       other object type to be returned
  Example    : $self->_type_by_external_db_id(1030, 'Translation');
  Description: Gets
               NOTE:  In a multi-species database, this method will
               return all the entries matching the search criteria, not
               just the ones associated with the current species.
  Returntype : list of dbIDs (gene_id, transcript_id, etc.)
  Exceptions : none
  Caller     : list_translation_ids_by_extids
               translationids_by_extids
  	       geneids_by_extids
  Status     : Stable

=cut

sub _type_by_external_db_id{
  my ($self, $external_db_id, $ensType, $extraType) = @_;

  my $from_sql = '';
  my $where_sql = '';
  my $ID_sql = "oxr.ensembl_id";

  if (defined $extraType) {
    if (lc($extraType) eq 'translation') {
      $ID_sql = "tl.translation_id";
    } else {
      $ID_sql = "t.${extraType}_id";
    }

    if (lc($ensType) eq 'translation') {
      $from_sql = 'transcript t, translation tl, ';
      $where_sql = qq(
          t.transcript_id = tl.transcript_id AND
          tl.translation_id = oxr.ensembl_id AND
          t.is_current = 1 AND
      );
    } else {
      $from_sql = 'transcript t, ';
      $where_sql = 't.'.lc($ensType).'_id = oxr.ensembl_id AND '.
          't.is_current = 1 AND ';
    }
  }

  if (lc($ensType) eq 'gene') {
    $from_sql = 'gene g, ';
    $where_sql = 'g.gene_id = oxr.ensembl_id AND g.is_current = 1 AND ';
  } elsif (lc($ensType) eq 'transcript') {
    $from_sql = 'transcript t, ';
    $where_sql = 't.transcript_id = oxr.ensembl_id AND t.is_current = 1 AND ';
  } elsif (lc($ensType) eq 'translation') {
    $from_sql = 'transcript t, translation tl, ';
    $where_sql = qq(
        t.transcript_id = tl.transcript_id AND
        tl.translation_id = oxr.ensembl_id AND
        t.is_current = 1 AND
    );
  }

  my $query = 
    "SELECT $ID_sql
       FROM $from_sql xref x, object_xref oxr
      WHERE $where_sql x.external_db_id = ? AND
  	     x.xref_id = oxr.xref_id AND oxr.ensembl_object_type= ?";

  my %result;

  my $sth = $self->prepare($query);

  $sth->bind_param( 1, "$external_db_id", SQL_VARCHAR );
  $sth->bind_param( 2, $ensType,          SQL_VARCHAR );
  $sth->execute();

  while ( my $r = $sth->fetchrow_array() ) { $result{$r} = 1 }

  return keys(%result);
}


=head2 fetch_all_by_description

  Arg [1]    : string description to search for. Include % etc in this string
  Arg [2]    : <optional> string $dbname. Name of the database to search

  Example    : @canc_refs = @{$db_entry_adaptor->fetch_all_by_description("%cancer%")};
               @db_entries = @{$db_entry_adaptor->fetch_all_by_description("%cancer%","MIM_MORBID")};
  Description: Retrieves DBEntries that match the description.
               Optionally you can search on external databases type.
               NOTE:  In a multi-species database, this method will
               return all the entries matching the search criteria, not
               just the ones associated with the current species.
  Returntype : ref to array of Bio::EnsEMBL::DBSQL::DBEntry
  Exceptions : None.
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_description {
  my ( $self, $description, $dbname ) = @_;

  my @results = ();

  my $sql =
    "SELECT xref.xref_id, xref.dbprimary_acc, xref.display_label,
           xref.version,
           exDB.dbprimary_acc_linkable, exDB.display_label_linkable, exDB.priority,
           exDB.db_name, exDB.db_display_name, exDB.db_release, es.synonym,
           xref.info_type, xref.info_text, exDB.type, exDB.secondary_db_name,
           exDB.secondary_db_table, xref.description
    FROM   (xref, external_db exDB)
    LEFT JOIN external_synonym es on es.xref_id = xref.xref_id
    WHERE  xref.description like ?
    AND    xref.external_db_id = exDB.external_db_id";

  if ( defined($dbname) ) { $sql .= " AND exDB.db_name = ? " }

  my $sth = $self->prepare($sql);

  $sth->bind_param( 1, $description, SQL_VARCHAR );

  if ( defined($dbname) ) {
    $sth->bind_param( 2, $dbname, SQL_VARCHAR );
  }

  $sth->execute();

  my $max_rows = 1000;

  while ( my $rowcache = $sth->fetchall_arrayref( undef, $max_rows ) ) {
    while ( my $arrayref = shift( @{$rowcache} ) ) {
      my ( $dbID,                $dbprimaryId,
           $displayid,           $version,
           $primary_id_linkable,
           $display_id_linkable, $priority,
           $ex_dbname,           $db_display_name,
           $release,             $synonym,
           $info_type,           $info_text,
           $type,                $secondary_db_name,
           $secondary_db_table,  $description
      ) = @$arrayref;

      my $exDB =
        Bio::EnsEMBL::DBEntry->new(
                           -adaptor             => $self,
                           -dbID                => $dbID,
                           -primary_id          => $dbprimaryId,
                           -display_id          => $displayid,
                           -version             => $version,
                           -release             => $release,
                           -dbname              => $ex_dbname,
                           -primary_id_linkable => $primary_id_linkable,
                           -display_id_linkable => $display_id_linkable,
                           -priority            => $priority,
                           -db_display_name     => $db_display_name,
                           -info_type           => $info_type,
                           -info_text           => $info_text,
                           -type                => $type,
                           -secondary_db_name   => $secondary_db_name,
                           -secondary_db_table  => $secondary_db_table,
                           -description         => $description
        );

      if ($synonym) { $exDB->add_synonym($synonym) }

      push @results, $exDB;

    } ## end while ( my $arrayref = shift...
  } ## end while ( my $rowcache = $sth...

  $sth->finish();

  return \@results;
} ## end sub fetch_all_by_description


=head2 fetch_all_by_source

  Arg [1]    : string source to search for. Include % etc in this string
               if you want to use SQL patterns

  Example    : @unigene_refs = @{$db_entry_adaptor->fetch_all_by_source("%unigene%")};
  Description: Retrieves DBEntrys that match the source name.               
  Returntype : ref to array of Bio::EnsEMBL::DBSQL::DBEntry
  Exceptions : None.
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_source {
  my ( $self, $source ) = @_;

  my @results = ();

  my $sql =
    "SELECT xref.xref_id, xref.dbprimary_acc, xref.display_label,
           xref.version,
           exDB.dbprimary_acc_linkable, exDB.display_label_linkable, exDB.priority,
           exDB.db_name, exDB.db_display_name, exDB.db_release, es.synonym,
           xref.info_type, xref.info_text, exDB.type, exDB.secondary_db_name,
           exDB.secondary_db_table, xref.description
    FROM   (xref, external_db exDB)
    LEFT JOIN external_synonym es on es.xref_id = xref.xref_id
    WHERE  exDB.db_name like ?
    AND    xref.external_db_id = exDB.external_db_id";


  my $sth = $self->prepare($sql);

  $sth->bind_param( 1, $source, SQL_VARCHAR );

  $sth->execute();

  my $max_rows = 1000;

  while ( my $rowcache = $sth->fetchall_arrayref( undef, $max_rows ) ) {
    while ( my $arrayref = shift( @{$rowcache} ) ) {
      my ( $dbID,                $dbprimaryId,
           $displayid,           $version,
           $primary_id_linkable,
           $display_id_linkable, $priority,
           $dbname,              $db_display_name,
           $release,             $synonym,
           $info_type,           $info_text,
           $type,                $secondary_db_name,
           $secondary_db_table,  $description
      ) = @$arrayref;

      my $exDB =
        Bio::EnsEMBL::DBEntry->new(
                           -adaptor             => $self,
                           -dbID                => $dbID,
                           -primary_id          => $dbprimaryId,
                           -display_id          => $displayid,
                           -version             => $version,
                           -release             => $release,
                           -dbname              => $dbname,
                           -primary_id_linkable => $primary_id_linkable,
                           -display_id_linkable => $display_id_linkable,
                           -priority            => $priority,
                           -db_display_name     => $db_display_name,
                           -info_type           => $info_type,
                           -info_text           => $info_text,
                           -type                => $type,
                           -secondary_db_name   => $secondary_db_name,
                           -secondary_db_table  => $secondary_db_table,
                           -description         => $description
        );

      if ($synonym) { $exDB->add_synonym($synonym) }

      push @results, $exDB;

    } ## end while ( my $arrayref = shift...
  } ## end while ( my $rowcache = $sth...

  $sth->finish();

  return \@results;
} ## end sub fetch_all_by_source


=head2 fetch_all_synonyms

  Arg [1]    : dbID of DBEntry to fetch synonyms for. Used in lazy loading of synonyms.

  Example    : @canc_refs = @{$db_entry_adaptor->fetch_all_synonyms(1234)};
  Description: Fetches the synonyms for a particular DBEntry.
  Returntype : listref of synonyms. List referred to may be empty if there are no synonyms.
  Exceptions : None.
  Caller     : General
  Status     : At Risk

=cut


sub fetch_all_synonyms {

  my ($self, $dbID) = @_;

  my @synonyms;

  my $sth = $self->prepare("SELECT synonym FROM external_synonym WHERE xref_id = ?");

  $sth->bind_param(1, $dbID, SQL_INTEGER);
  $sth->execute();
  while ( my $arrayref = $sth->fetchrow_arrayref()){
    my ($synonym) = @$arrayref;
    push (@synonyms, $synonym);
  }

  @synonyms = () if (!@synonyms);

  return \@synonyms;

}


=head2 get_db_name_from_external_db_id

  Arg [1]    : external_dbid of database to get the database_name
  Example    : my $db_name = $db_entry_adaptor->get_db_name_from_external_db_id(1100);
  Description: Gets the database name for a certain external_db_id
  Returntype : scalar
  Exceptions : None.
  Caller     : General
  Status     : At Risk

=cut

sub get_db_name_from_external_db_id{
    my $self = shift;
    my $external_db_id = shift;

    my $sth = $self->prepare("SELECT db_name FROM external_db WHERE external_db_id = ?");

    $sth->bind_param(1, $external_db_id, SQL_INTEGER);
    $sth->execute();
    my ($db_name) = $sth->fetchrow_array();
    $sth->finish();
    return $db_name;

}

=head2 geneids_by_extids

  Description: DEPRECATED use list_gene_ids_by_extids instead

=cut

sub geneids_by_extids{
   my ($self,$name) = @_;
   deprecate(" use 'list_gene_ids_by_extids instead");
   return $self->list_gene_ids_by_extids( $name );
}


=head2 translationids_by_extids

  Description:  DEPRECATED use list_translation_ids_by_extids instead

=cut

sub translationids_by_extids{
  my ($self,$name) = @_;
  deprecate("Use list_translation_ids_by_extids instead");
  return $self->list_translation_ids_by_extids( $name );
}


=head2 transcriptids_by_extids

  Description: DEPRECATED use transcriptids_by_extids instead

=cut

sub transcriptids_by_extids{
   my ($self,$name) = @_;
   deprecate("Use list_transcript_ids_by_extids instead.");
   return $self->list_transcript_ids_by_extids( $name );
}


1;

