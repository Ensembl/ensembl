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

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk

=head1 APPENDIX

=cut

;

package Bio::EnsEMBL::DBSQL::DBEntryAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::IdentityXref;

use vars qw(@ISA);
use strict;

@ISA = qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );

sub fetch_by_dbID {
  my ($self, $dbID ) = @_;
  
  my $sth = $self->prepare( "
    SELECT Xref.xrefId, Xref.dbprimary_id, Xref.display_id,
           Xref.version, Xref.description,
           exDB.db_name, exDB.release
      FROM Xref, externalDB exDB
     WHERE Xref.xrefId = $dbID
       AND Xref.externalDBId = exDb.externalDBId 
   " );

  $sth->execute();
  my ( $refID, $dbprimaryId, $displayid, $version, $desc, $dbname, $release) =
    $sth->fetchrow_array();

  if( ! defined $refID ) {
    return undef;
  }

  my $exDB = Bio::EnsEMBL::DBEntry->new
    ( -adaptor => $self,
      -dbID => $dbID,
      -primary_id => $dbprimaryId,
      -display_id => $displayid,
      -version => $version,
      -release => $release,
      -dbname => $dbname );
  
  if( $desc ) {
    $exDB->description( $desc );
  }

  
  my $get_synonym = $self->prepare( "
    SELECT synonym 
      FROM externalSynonym
     WHERE xrefId = $dbID
  " );
  $get_synonym->execute();
  
  while( my ($synonym) = $get_synonym->fetchrow_array() ) {
    $exDB->add_synonym( $synonym );
  }

  return $exDB;
}


sub store {
    my ( $self, $exObj, $ensObject, $ensType ) = @_;
    
    # $self->throw( "Sorry, store not yet supported" );
    my $dbUnknown;
    
    # check if db exists
    # urlPattern dbname release
    my $sth = $self->prepare( "
     SELECT externalDBId
       FROM externalDB
      WHERE db_name = ?
        AND release = ?
    " );
    $sth->execute( $exObj->dbname(), $exObj->release() );
    
    my $dbRef;
    
    if(  ($dbRef) =  $sth->fetchrow_array() ) {
    
    } else {
	# store it, get dbID for that
	$sth = $self->prepare( "
       INSERT INTO externalDB 
       SET db_name = ?,
           release = ?
     " );
	$sth->execute( $exObj->dbname(), $exObj->release());
	
	$dbUnknown = 1;
	$sth = $self->prepare( "
       SELECT LAST_INSERT_ID()
     " );
	$sth->execute();
	( $dbRef ) = $sth->fetchrow_array();
	if( ! defined $dbRef ) {
	    $self->throw( "Database entry failed." );
	}
    }
    
    my $dbX;
    
    if( ! $dbUnknown ) {
	$sth = $self->prepare( "
       SELECT xrefId
         FROM Xref
        WHERE externalDBId = ?
          AND dbprimary_id = ?
          AND version = ?
     " );
	$sth->execute( $dbRef, $exObj->primary_id(), 
		       $exObj->version() );
	( $dbX ) = $sth->fetchrow_array();
    } else {
	# dont check for existence
    }
    
    if( ! defined $dbX ) {
	
	$sth = $self->prepare( "
      INSERT INTO Xref 
       SET dbprimary_id = ?,
           display_id = ?,
           version = ?,
           description = ?,
           externalDBId = $dbRef
     " );
	$sth->execute( $exObj->primary_id(), $exObj->display_id(), $exObj->version(),
		       $exObj->description());
	
	$sth = $self->prepare( "
      SELECT LAST_INSERT_ID()
    " );
	$sth->execute();
	( $dbX ) = $sth->fetchrow_array();
	
	# synonyms
	
	
	
	my @synonyms = $exObj->get_synonyms();
	foreach my $syn ( @synonyms ) {
	    
#Check if this synonym is already in the database for the given primary id
	    my $sth = $self->prepare( "
     SELECT xrefId,
            synonym
       FROM externalSynonym
      WHERE xrefId = '$dbX'
        AND synonym = '$syn'
    " );
	    $sth->execute;
	    
	    my @dbSyn = $sth->fetchrow_array();
	    
	    #print STDERR $dbSyn[0],"\n";
	    
	    if( ! @dbSyn ) {
		
		
		$sth = $self->prepare( "
        INSERT INTO externalSynonym
         SET xrefId = $dbX,
            synonym = '$syn'
      " );
		$sth->execute();
	    }
	}
	
	
	$sth = $self->prepare( "
   INSERT INTO objectXref
     SET xrefId = $dbX,
         ensembl_object_type = ?,
         ensembl_id = ?
  " );
	
	$sth->execute( $ensType, $ensObject );
	
	$exObj->dbID( $dbX );
	$exObj->adaptor( $self );
	
	if ($exObj->isa('Bio::EnsEMBL::IdentityXref')) {
	    $sth = $self->prepare( "
      SELECT LAST_INSERT_ID()
    " );
	    $sth->execute();
	    my ( $Xidt ) = $sth->fetchrow_array();
	    
	    $sth = $self->prepare( "
             INSERT INTO identityXref
             SET objectxrefId = $Xidt,
             query_identity = ?,
             target_identity = ?
    " );
	    $sth->execute( $exObj->query_identity, $exObj->target_identity );
	    
	}
    } else {
	# line is already in Xref table. Need to add to objectXref
	$sth = $self->prepare( "
             INSERT INTO objectXref
               SET xrefId = $dbX,
               ensembl_object_type = ?,
               ensembl_id = ?");
	
	$sth->execute( $ensType, $ensObject );
	
	$exObj->dbID( $dbX );
	$exObj->adaptor( $self );
    }
        
    return $dbX;
    
}

=head2 fetch_by_union

 Title   : fetch_by_union
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_by_union{
   my ($self,$id) = @_;
   my @out = ();

   my @xg = $self->fetch_by_gene($id);
   push(@out,@xg);

   my @xrc = $self->fetch_by_rawContig($id);
   push (@out,@xrc);

   my @xtranscript = $self->fetch_by_transcript($id);
   push(@out,@xtranscript);

   my @xtranslation = $self->fetch_by_translation($id);
   push(@out,@xtranslation);

   return @out;
}


sub fetch_by_gene {
  my ( $self, $geneId ) = @_;
  return $self->_fetch_by_EnsObject_type( $geneId, 'Gene' );
}

sub fetch_by_rawContig {
  my ( $self, $rawContigId ) = @_;
  return $self->_fetch_by_EnsObject_type( $rawContigId, 'RawContig' );
}

sub fetch_by_transcript {
  my ( $self, $trscId ) = @_;
  return $self->_fetch_by_EnsObject_type( $trscId, 'Transcript' );
}

sub fetch_by_translation {
  my ( $self, $trslId ) = @_;
  return $self->_fetch_by_EnsObject_type( $trslId, 'Translation' );
}


sub _fetch_by_EnsObject_type {
  my ( $self, $ensObj, $ensType ) = @_;
  my @out;

  my $sth = $self->prepare( "
    SELECT Xref.xrefId, Xref.dbprimary_id, Xref.display_id,
           Xref.version, Xref.description,
           exDB.db_name, exDB.release, oxr.objectxrefId
      FROM Xref, externalDB exDB, objectXref oxr 
     WHERE Xref.xrefId = oxr.xrefId
       AND Xref.externalDBId = exDB.externalDBId 
       AND oxr.ensembl_id = '$ensObj'
       AND oxr.ensembl_object_type = '$ensType'
   " );

  $sth->execute();

 
  while ( my $arrRef = $sth->fetchrow_arrayref() ) {
    my ( $refID, $dbprimaryId, $displayid, $version, $desc, $dbname, $release, $objid ) =
      @$arrRef;;
    
    my $exDB;

    my $sth1 = $self->prepare( "
      SELECT idt.query_identity , idt.target_identity
      FROM   identityXref idt, objectXref oxr
      WHERE  idt.objectxrefId = '$objid'
    ");
 
    $sth1->execute();
    my ($queryid, $targetid) = $sth1->fetchrow_array();

    if ((defined $queryid) || (defined $targetid)) {
	$exDB = Bio::EnsEMBL::IdentityXref->new
	    ( -adaptor => $self,
	      -dbID => $refID,
	      -primary_id => $dbprimaryId,
	      -display_id => $displayid,
	      -version => $version,
	      -release => $release,
	      -dbname => $dbname);
	      
	$exDB->query_identity($queryid);
	$exDB->target_identity($targetid);

    }
    
    else {
	$exDB = Bio::EnsEMBL::DBEntry->new
	    ( -adaptor => $self,
	      -dbID => $refID,
	      -primary_id => $dbprimaryId,
	      -display_id => $displayid,
	      -version => $version,
	      -release => $release,
	      -dbname => $dbname );
    }
    
    if( $desc ) {
	$exDB->description( $desc );
    }

    
    my $sth = $self->prepare( "
      SELECT synonym 
        FROM externalSynonym
       WHERE xrefId = $refID
    " );
    $sth->execute();
  
    while( my ($synonym) = $sth->fetchrow_array() ) {
      $exDB->add_synonym( $synonym );
    }
    push( @out, $exDB );
  }

  return @out;
}

=head2 geneids_by_extids

 Title   : geneids_by_extids
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub geneids_by_extids{
   my ($self,$name) = @_;
   my @genes;

   my $sth = $self->prepare("select tr.gene from transcript tr,Xref x, objectXref oxr where tr.translation = oxr.ensembl_id and oxr.xrefId = x.xrefId and x.display_id = '$name'");
   $sth->execute();

   while( ($a) = $sth->fetchrow_array ) {
       push(@genes,$a);
   }

   return @genes;
}

=head2 transcriptids_by_extids

 Title   : transcriptids_by_extids
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub transcriptids_by_extids{
   my ($self,$name) = @_;
   my @transcripts;
   my @translations = $self->_type_by_external_id($name,'Translation');

foreach my $t (@translations) {
       my $sth = $self->prepare( "select id from transcript where translation = '$t'");
       $sth->execute();
       my $tr = $sth->fetchrow;
       push (@transcripts,$tr);
   }
   return @transcripts;

}

=head2 translationids_by_extids

 Title   : translationids_by_extids
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub translationids_by_extids{
    my ($self,$name) = @_;
    return $self->_type_by_external_id($name,'Transcript');
}

=head2 rawContigids_by_extids

 Title   : rawContigids_by_extids
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub rawContigids_by_extids{
   my ($self,$name) = @_;
   return $self->_type_by_external_id($name,'rawContig');
}



=head2 _type_by_external_id

 Title   : _fetch_type_by_external_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _type_by_external_id{
   my ($self,$name,$ensType) = @_;
   
my @out;


  my $sth = $self->prepare( "
    SELECT oxr.ensembl_id
    FROM Xref, externalDB exDB, objectXref oxr, externalSynonym syn 
     WHERE (Xref.dbprimary_id = '$name'
            AND Xref.xrefId = oxr.xrefId
            AND oxr.ensembl_object_type = '$ensType')
     OR    (Xref.display_id = '$name'
            AND Xref.xrefId = oxr.xrefId
            AND oxr.ensembl_object_type = '$ensType')
     OR    (syn.synonym = '$name'
            AND syn.xrefId = oxr.xrefId
            AND oxr.ensembl_object_type = '$ensType')
         " );

   $sth->execute();
   while ( my $arrRef = $sth->fetchrow_arrayref() ) {
       my ( $ensID) =
	   @$arrRef;;
       push (@out,$ensID);
     }
   
   return @out;
}


# creates all tables for this adaptor
# if they exist they are emptied and newly created
sub create_tables {
  my $self = shift;

  my $sth = $self->prepare( "drop table if exists objectXref, Xref, externalDescription, externalSynonym, externalDB" );
  $sth->execute();

  $sth = $self->prepare( qq{
     CREATE TABLE objectXref(
       ensembl_id VARCHAR(40) not null, 
       ensembl_object_type ENUM( 'RawContig', 'Transcript', 'Gene', 'Translation' ) not null,
       xrefId INT not null,
       PRIMARY KEY( ensembl_object_type, ensembl_id, xrefId ),
       KEY xrefIdx( xrefId, ensembl_object_type, ensembl_id )
     )
   } );
  $sth->execute();
  $sth = $self->prepare( qq{
     CREATE TABLE Xref(
         xrefId INT not null auto_increment,
         externalDBId int not null,
         dbprimary_id VARCHAR(40) not null,
	 display_id VARCHAR(40) not null,
         version VARCHAR(10),
	 description VARCHAR(255),
         PRIMARY KEY( xrefId ),
         KEY idIdx( dbprimary_id ))
   } );

  $sth->execute();

  $sth = $self->prepare( qq{
     CREATE TABLE externalSynonym(
         xrefId INT not null,
         synonym VARCHAR(40) not null,
         PRIMARY KEY( xrefId, synonym ),
	 KEY nameIdx( synonym )) 
   } );
  $sth->execute();

  $sth = $self->prepare( qq{
     CREATE TABLE externalDB(
         externalDBId INT not null auto_increment,
         db_name VARCHAR(40) not null,
	 release VARCHAR(40),
         PRIMARY KEY( externalDBId ) ) 
   } );
  $sth->execute();
}


1;


__END__


# remove the tables from database
sub delete_tables {
}

# check if tables exist
sub exists_tables {
}

ObjectXref
=============
ensembl_id varchar, later int
ensembl_object_type  enum 
xrefId int
primary key (ensembl_id,ensembl_object_type,xrefId) 


Xref
=================
xrefId int (autogenerated) 
externalDBId int
dbprimary_id  varchar
version varchar

primary key (xrefId)

ExternalDescription
=======================
xrefId int
description varchar (256)

primary key (xrefId)

ExternalSynonym
=================
xrefId int
synonym varchar

primary key (external_id,synonym)


ExternalDB
===================
externalDBId int
db_name varchar
release varchar

