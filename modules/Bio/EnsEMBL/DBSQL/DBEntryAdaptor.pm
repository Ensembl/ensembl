# EnsEMBL External object reference reading writing adaptor for mySQL
#
# Copyright EMBL-EBI 2001
#
# Author: Arne Stabenau
# 
# Date : 06.03.2001
#

=head1 NAME

Bio::EnsEMBL::DBSQL::DBEntryAdaptor - 
MySQL Database queries to load and store external object references.

=head1 SYNOPSIS

$db_entry_adaptor = $db_adaptor->get_DBEntryAdaptor();
$dbEntry = $db_entry_adaptor->fetch_by_dbID($id);

=head1 CONTACT

Post questions to the EnsEMBL developer list <ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::DBSQL::DBEntryAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::IdentityXref;
use Bio::EnsEMBL::GoXref;

use vars qw(@ISA);
use strict;

@ISA = qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );


=head2 fetch_by_dbID

  Arg [1]    : int $dbID
               the unique database identifier for the DBEntry to retrieve
  Example    : my $db_entry = $db_entry_adaptor->fetch_by_dbID($dbID);
  Description: retrieves a dbEntry from the database via its unique identifier 
  Returntype : Bio::EnsEMBL::DBEntry
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_dbID {
  my ($self, $dbID ) = @_;

  my $sth = $self->prepare(
   "SELECT xref.xref_id, xref.dbprimary_acc, xref.display_label,
           xref.version, xref.description,
           exDB.db_name, exDB.release, es.synonym
    FROM   xref, external_db exDB
    LEFT JOIN external_synonym es on es.xref_id = xref.xref_id
    WHERE  xref.xref_id = ?
    AND    xref.external_db_id = exDB.external_db_id");

  $sth->execute($dbID);

  my $exDB;
  my %duplicate;

  while ( my $arrayref = $sth->fetchrow_arrayref()){
    my ( $refID, $dbprimaryId, $displayid, $version, $desc, $dbname, 
	 $release, $synonym) = @$arrayref;
    return undef if( ! defined $refID );

    unless ($duplicate{$refID}){
      $duplicate{$refID} = 1;

      $exDB = Bio::EnsEMBL::DBEntry->new
	( -adaptor => $self,
	  -dbID => $dbID,
	  -primary_id => $dbprimaryId,
	  -display_id => $displayid,
	  -version => $version,
	  -release => $release,
	  -dbname => $dbname );

      $exDB->description( $desc ) if ( $desc );
    } # end duplicate

    $exDB->add_synonym( $synonym )  if ($synonym);
  } # end while

  return $exDB;
}


=head2 store

  Arg [1]    : ?? $exObj
  Arg [2]    : ?? $ensObject
  Arg [3]    : ?? $ensType
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut


sub store {
  my ( $self, $exObj, $ensObject, $ensType ) = @_;
  my $dbJustInserted;

  #
  # Check for the existance of the external_db, throw if it does not exist
  #
  my $sth = $self->prepare( "
     SELECT external_db_id
       FROM external_db
      WHERE db_name = ?
        AND release = ?");

  $sth->execute( $exObj->dbname(), $exObj->release() );

  my ($dbRef) =  $sth->fetchrow_array();

  if(!$dbRef) {
    $self->throw("external_db [" . $exObj->dbname . "] release [" .
		 $exObj->release . "] does not exist");
  }

  #
  # Check for the existance of the external reference, add it if not present
  #
  $sth = $self->prepare( "
       SELECT xref_id
         FROM xref
        WHERE external_db_id = ?
          AND dbprimary_acc = ?
          AND version = ?" );

  $sth->execute( $dbRef, $exObj->primary_id(),$exObj->version() );
  my ($dbX) = $sth->fetchrow_array();
  $sth->finish();

  if(!$dbX) {
      #
      # store the new xref
      #
    $sth = $self->prepare( "
       INSERT ignore INTO xref 
       SET dbprimary_acc = ?,
           display_label = ?,
           version = ?,
           description = ?,
           external_db_id = ?" );
    $sth->execute( $exObj->primary_id(), $exObj->display_id(), 
		   $exObj->version(), $exObj->description(), $dbRef);
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
      $synonym_check_sth->execute($dbX, $syn);
      my ($dbSyn) = $synonym_check_sth->fetchrow_array(); 
      $synonym_store_sth->execute($dbX, $syn) if(!$dbSyn);
    }

    $synonym_check_sth->finish();
    $synonym_store_sth->finish();
  }

  #
  # check if the object mapping was already stored
  #
  $sth = $self->prepare (
           "SELECT xref_id
            FROM object_xref
            WHERE xref_id = ?
            AND   ensembl_object_type = ?
            AND   ensembl_id = ?");

  $sth->execute($dbX, $ensType, $ensObject);
  my ($tst) = $sth->fetchrow_array;
  $sth->finish();

  if(!$tst) {
    #
    # Store the reference to the internal ensembl object
    #
    $sth = $self->prepare(
         "INSERT ignore INTO object_xref
          SET xref_id = ?, ensembl_object_type = ?, ensembl_id = ?");

    $sth->execute( $dbX, $ensType, $ensObject );	
    $exObj->dbID( $dbX );
    $exObj->adaptor( $self );

    my $Xidt = $sth->{'mysql_insertid'};

    #
    # If this is an IdentityXref need to store in that table too
    # If its GoXref add the linkage type to go_xref table
    #
    if ($exObj->isa('Bio::EnsEMBL::IdentityXref')) {
      
      my $analysis_id;
      if( $exObj->analysis() ) {
	$analysis_id = $self->db()->get_AnalysisAdaptor()->
	  store( $exObj->analysis() );
      } else {
	$analysis_id = undef;
      }
      
      $sth = $self->prepare( "
             INSERT ignore INTO identity_xref
             SET object_xref_id = ?,
             query_identity = ?,
             target_identity = ?,
             hit_start = ?,
             hit_end = ?,
             translation_start = ?,
             translation_end = ?,
             cigar_line = ?,
             score = ?,
             evalue = ?,
             analysis_id = ?" );

      $sth->execute($Xidt, $exObj->query_identity, $exObj->target_identity,
		   $exObj->hit_start(), $exObj->hit_end(), 
		   $exObj->translation_start(), $exObj->translation_end(),
		   $exObj->cigar_line(), $exObj->score(), $exObj->evalue(),
		   $analysis_id);
    } elsif( $exObj->isa( 'Bio::EnsEMBL::GoXref' )) {
      $sth = $self->prepare( "
             INSERT ignore INTO go_xref
                SET object_xref_id = ?,
                    linkage_type = ? " );
      foreach my $lt (@{$exObj->get_all_linkage_types()}) {
        $sth->execute( $Xidt, $lt );
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

=cut

sub exists {
  my ($self, $dbe) = @_ ;

  unless($dbe && ref $dbe && $dbe->isa('Bio::EnsEMBL::DBEntry')) {
    $self->throw("arg must be a Bio::EnsEMBL::DBEntry not [$dbe]");
  }

  my $sth = $self->prepare('SELECT x.xref_id 
                            FROM   xref x, external_db xdb
                            WHERE  x.external_db_id = xdb.external_db_id
                            AND    x.display_label = ? 
                            AND    xdb.db_name = ?');

  $sth->execute($dbe->display_id, $dbe->dbname);

  my ($dbID) = $sth->fetchrow_array;

  $sth->finish;

  return $dbID;
}



=head2 fetch_all_by_Gene

  Arg [1]    : Bio::EnsEMBL::Gene $gene 
               (The gene to retrieve DBEntries for)
  Example    : @db_entries = @{$db_entry_adaptor->fetch_by_Gene($gene)};
  Description: This returns a list of DBEntries associated with this gene.
               Note that this method was changed in release 15.  Previously
               it set the DBLinks attribute of the gene passed in to contain
               all of the gene, transcript, and translation xrefs associated 
               with this gene.
  Returntype : listref of Bio::EnsEMBL::DBEntries
  Exceptions : none
  Caller     : Bio::EnsEMBL::Gene

=cut

sub fetch_all_by_Gene {
  my ( $self, $gene ) = @_;

  return $self->_fetch_by_object_type($gene->dbID(), 'Gene');
}


=head2 fetch_all_by_RawContig

  Arg [1]    : Bio::EnsEMBL::RawContig $contig
  Example    : @db_entries = @{$db_entry_adaptor->fetch_by_RawContig($contig)}
  Description: Retrieves a list of DBentries by a contig object. Since 
               xrefs are never stored in EnsEMBL for RawContigs this method
               never returns anything.  It has been left in, because of the
               possibility that it is used externally, and because you could
               in theory store xrefs on rawcontigs if you wanted to.
  Returntype : listref of Bio::EnsEMBL::DBEntries
  Exceptions : none
  Caller     : general

=cut

sub fetch_all_by_RawContig {
  my ( $self, $contig ) = @_;

  return $self->_fetch_by_object_type( $contig->dbID(), 'RawContig' );
}


=head2 fetch_all_by_Transcript

  Arg [1]    : Bio::EnsEMBL::Transcript
  Example    : @db_entries = @{$db_entry_adaptor->fetch_by_Gene($trans)};
  Description: This returns a list of DBEntries associated with this 
               transcript. Note that this method was changed in release 15.  
               Previously it set the DBLinks attribute of the gene passed in 
               to contain all of the gene, transcript, and translation xrefs 
               associated with this gene.
  Returntype : listref of Bio::EnsEMBL::DBEntries
  Exceptions : none
  Caller     : Bio::EnsEMBL::Gene 

=cut

sub fetch_all_by_Transcript {
  my ( $self, $trans ) = @_;

  return $self->_fetch_by_object_type( $trans->dbID(), 'Transcript');
}


=head2 fetch_all_by_Translation

  Arg [1]    : Bio::EnsEMBL::Translation $trans
               (The translation to fetch database entries for)
  Example    : @db_entries = @{$db_entry_adptr->fetch_by_Translation($trans)};
  Description: Retrieves external database entries for an EnsEMBL translation
  Returntype : listref of Bio::EnsEMBL::DBEntries
  Exceptions : none
  Caller     : general

=cut

sub fetch_all_by_Translation {
  my ( $self, $trans ) = @_;

  return $self->_fetch_by_object_type( $trans->dbID(), 'Translation' );
}


=head2 fetch_by_object_type

  Arg [1]    : string $ensObj
  Arg [2]    : string $ensType
  			   (object type to be returned) 
  Example    : $self->_fetch_by_object_type( $translation_id, 'Translation' )
  Description: Fetches DBEntry by Object type
  Returntype : arrayref of DBEntry objects
  Exceptions : none
  Caller     : fetch_all_by_Gene
               fetch_all_by_Translation
               fetch_all_by_Transcript
               fetch_all_by_RawContig

=cut

sub _fetch_by_object_type {
  my ( $self, $ensObj, $ensType ) = @_;
  my @out;

  if (!defined($ensObj)) {
    $self->throw("Can't fetch_by_EnsObject_type without an object");
  }
  if (!defined($ensType)) {
    $self->throw("Can't fetch_by_EnsObject_type without a type");
  }
  my $sth = $self->prepare("
    SELECT xref.xref_id, xref.dbprimary_acc, xref.display_label, xref.version,
           xref.description,
           exDB.db_name, exDB.release, exDB.status, 
           oxr.object_xref_id, 
           es.synonym, 
           idt.query_identity, idt.target_identity, idt.hit_start, idt.hit_end,
           idt.translation_start, idt.translation_end, idt.cigar_line,
           idt.score, idt.evalue, idt.analysis_id,
           gx.linkage_type
    FROM   xref xref, external_db exDB, object_xref oxr 
    LEFT JOIN external_synonym es on es.xref_id = xref.xref_id 
    LEFT JOIN identity_xref idt on idt.object_xref_id = oxr.object_xref_id
    LEFT JOIN go_xref gx on gx.object_xref_id = oxr.object_xref_id
    WHERE  xref.xref_id = oxr.xref_id
      AND  xref.external_db_id = exDB.external_db_id 
      AND  oxr.ensembl_id = ?
      AND  oxr.ensembl_object_type = ?
  ");

  $sth->execute($ensObj, $ensType);
  my (%seen, %linkage_types, %synonyms);
  

  while ( my $arrRef = $sth->fetchrow_arrayref() ) {
    my ( $refID, $dbprimaryId, $displayid, $version, 
         $desc, $dbname, $release, $exDB_status, $objid, 
         $synonym, $queryid, $targetid, $hit_start, $hit_end,
	 $translation_start, $translation_end, $cigar_line,
	 $score, $evalue, $analysis_id, $linkage_type ) = @$arrRef;

    my %obj_hash = ( 
		    _adaptor    => $self,
		    _dbID       => $refID,
		    _primary_id => $dbprimaryId,
		    _display_id => $displayid,
		    _version    => $version,
		    _release    => $release,
		    _dbname     => $dbname);


    # using an outer join on the synonyms as well as on identity_xref, we
    # now have to filter out the duplicates (see v.1.18 for
    # original). Since there is at most one identity_xref row per xref,
    # this is easy enough; all the 'extra' bits are synonyms
    if ( !$seen{$refID} )  {
      my $exDB;

      if ((defined $queryid)) {         # an xref with similarity scores
        $exDB = Bio::EnsEMBL::IdentityXref->new_fast(\%obj_hash);       
        $exDB->query_identity($queryid);
        $exDB->target_identity($targetid);
	if( $analysis_id ) {
	  my $analysis = $self->db()->get_AnalysisAdaptor()->
	    fetch_by_dbID( $analysis_id );
	  if( $analysis ) {
	    $exDB->analysis( $analysis );
	  }
	}
	$exDB->cigar_line( $cigar_line );
	$exDB->hit_start( $hit_start );
	$exDB->hit_end( $hit_end );
	$exDB->translation_start( $translation_start );
	$exDB->translation_end( $translation_end );
	$exDB->score( $score );
	$exDB->evalue( $evalue );
      } elsif( defined $linkage_type ) {
        $exDB = Bio::EnsEMBL::GoXref->new_fast( \%obj_hash );
        $exDB->add_linkage_type( $linkage_type );
        $linkage_types{$refID}->{$linkage_type} = 1;
      } else {
        $exDB = Bio::EnsEMBL::DBEntry->new_fast(\%obj_hash);
      }

      $exDB->description($desc)   if(defined($desc));
      $exDB->status($exDB_status) if(defined($exDB_status));

      push( @out, $exDB );
      $seen{$refID} = $exDB;
    } 

    #
    # $exDB still points to the same xref, so we can keep adding
    # go evidence tags or synonyms
    #
    if(defined($synonym) && !$synonyms{$refID}->{$synonym}) {
      $seen{$refID}->add_synonym($synonym) if(defined($synonym));
      $synonyms{$refID}->{$synonym} = 1;
    }

    if(defined($linkage_type) && !$linkage_types{$refID}->{$linkage_type}) {
      $seen{$refID}->add_linkage_type($linkage_type);
      $linkage_types{$refID}->{$linkage_type} = 1;
    }
  }

  return \@out;
}

=head2 list_gene_ids_by_extids

  Arg [1]    : string $external_id
  Example    : none
  Description: Retrieve a list of geneid by an external identifier that is linked to 
               any of the genes transcripts, translations or the gene itself 
  Returntype : listref of strings
  Exceptions : none
  Caller     : unknown

=cut

sub list_gene_ids_by_extids{
   my ($self,$name) = @_;

   my %T = map { ($_,1) }
       $self->_type_by_external_id( $name, 'Translation', 'gene' ),
       $self->_type_by_external_id( $name, 'Transcript',  'gene' ),
       $self->_type_by_external_id( $name, 'Gene' );
   return keys %T;
}

=head2 geneids_by_extids

  Arg [1]    : string $external_id
  Example    : none
  Description: Retrieve a list of geneid by an external identifier that is linked to 
               any of the genes transcripts, translations or the gene itself 
			   (please not that this call is deprecated)
  Returntype : listref of strings
  Exceptions : none
  Caller     : unknown

=cut

sub geneids_by_extids{
   my ($self,$name) = @_;
   warn ("This method is deprecated please use 'list_gene_ids_by_extids");
   return $self->list_gene_ids_by_extids( $name );
}

=head2 list_transcript_ids_by_extids

  Arg [1]    : string $external_id
  Example    : none
  Description: Retrieve a list transcriptid by an external identifier that is linked to 
               any of the genes transcripts, translations or the gene itself 
  Returntype : listref of strings
  Exceptions : none
  Caller     : unknown

=cut

sub list_transcript_ids_by_extids{
   my ($self,$name) = @_;
   my @transcripts;

   my %T = map { ($_,1) }
       $self->_type_by_external_id( $name, 'Translation', 'transcript' ),
       $self->_type_by_external_id( $name, 'Transcript' );
   return keys %T;
}


=head2 transcriptids_by_extids

  Arg [1]    : string $external_id
  Example    : none
  Description: Retrieve a list transcriptid by an external identifier that is linked to 
               any of the genes transcripts, translations or the gene itself 
  Returntype : listref of strings
  Exceptions : none
  Caller     : unknown

=cut

sub transcriptids_by_extids{
   my ($self,$name) = @_;
   warn ("This method is deprecated please use 'list_transcript_ids_by_extids");
   return $self->list_transcript_ids_by_extids( $name );
}


=head2 translationids_by_extids

  Arg [1]    :  string $name 
  Example    :  none
  Description:  Gets a list of translation IDs by external display IDs 
  				(please note that this call is deprecated)
  Returntype :  list of Ints
  Exceptions :  none
  Caller     :  unknown

=cut

sub translationids_by_extids{
  my ($self,$name) = @_;
  warn ("This method is deprecated please use 'list_translation_ids_by_extids");
  return $self->list_translation_ids_by_extids( $name );
}


=head2 list_translation_ids_by_extids

  Arg [1]    :  string $name 
  Example    :  none
  Description:  Gets a list of translation IDs by external display IDs
  Returntype :  list of Ints
  Exceptions :  none
  Caller     :  unknown

=cut

sub list_translation_ids_by_extids{
  my ($self,$name) = @_;
  return $self->_type_by_external_id( $name, 'Translation' );
}

=head2 _type_by_external_id

  Arg [1]    : string $name
  			   (dbprimary_acc)
  Arg [2]    : string $ensType
  			   (Object_type)
  Arg [3]    : string $extraType
  			   (other object type to be returned) - optional
  Example    : $self->_type_by_external_id( $name, 'Translation' ) 
  Description: Gets
  Returntype : list of ensembl_IDs
  Exceptions : none
  Caller     : list_translation_ids_by_extids
               translationids_by_extids
  			   geneids_by_extids

=cut

sub _type_by_external_id{
  my ($self,$name,$ensType,$extraType) = @_;

  my $from_sql = '';
  my $where_sql = '';
  my $ID_sql = "oxr.ensembl_id";
  if(defined $extraType) {
    $ID_sql = "t.${extraType}_id";
    $from_sql = 'transcript as t, ';
    $where_sql = 't.'.lc($ensType).'_id = oxr.ensembl_id and ';
  }
  my @queries = (
    "select $ID_sql
       from $from_sql xref, object_xref as oxr
      where $where_sql xref.dbprimary_acc = ? and
            xref.xref_id = oxr.xref_id and oxr.ensembl_object_type= ?",
    "select $ID_sql 
       from $from_sql xref, object_xref as oxr
      where $where_sql xref.display_label = ? and
            xref.xref_id = oxr.xref_id and oxr.ensembl_object_type= ?",
    "select $ID_sql
       from $from_sql object_xref as oxr, external_synonym as syn
      where $where_sql syn.synonym = ? and
            syn.xref_id = oxr.xref_id and oxr.ensembl_object_type= ?",
  );

# Increase speed of query by splitting the OR in query into three separate 
# queries. This is because the 'or' statments render the index useless 
# because MySQL can't use any fields in the index.

  my %hash = (); 
  foreach( @queries ) {
    my $sth = $self->prepare( $_ );
    $sth->execute("$name", $ensType);
    while( my $r = $sth->fetchrow_array() ) {
      $hash{$r} = 1;
    }
  }
  return keys %hash;
}


1;
