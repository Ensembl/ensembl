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
  ## this function may need revamping in the same way as
  ##  _fetch_by_EnsObject_type, using double outer join to get all stuff in
  ##  one go. PL)

  my ($self, $dbID ) = @_;
  
  my $sth = $self->prepare( "
    SELECT xref.xref_id, xref.dbprimary_acc, xref.display_label,
           xref.version, xref.description,
           exDB.db_name, exDB.release
      FROM xref, external_db exDB
     WHERE xref.xref_id = $dbID
       AND xref.external_db_id = exDB.external_db_id 
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
      FROM external_synonym
     WHERE xref_id = $dbID
  " );
  $get_synonym->execute();
  
  while( my ($synonym) = $get_synonym->fetchrow_array() ) {
    $exDB->add_synonym( $synonym );
  }

  return $exDB;
}


sub store {
    my ( $self, $exObj, $ensObject, $ensType ) = @_;
    #print STDERR"storing DBEntryAdaptor\n";
    # $self->throw( "Sorry, store not yet supported" );
    my $dbJustInserted;
    
    # check if db exists
    # urlPattern dbname release
    my $sth = $self->prepare( "
     SELECT external_db_id
       FROM external_db
      WHERE db_name = ?
        AND release = ?
    " );
    $sth->execute( $exObj->dbname(), $exObj->release() );
    
    my $dbRef;
    
    if(  ($dbRef) =  $sth->fetchrow_array() ) {
        $dbJustInserted = 0;
    } else {
	# store it, get dbID for that
	$sth = $self->prepare( "
       INSERT INTO external_db 
       SET db_name = ?,
           release = ?,
           status  = ?
     " );
	print STDERR "STATUS DBE: ".$exObj->status."\n";

	$sth->execute( $exObj->dbname(), $exObj->release(), $exObj->status);
	
	$dbJustInserted = 1;
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
    
    if(  $dbJustInserted ) {
	# dont have to check for existence; cannnot have been inserted at
	# this point, so $dbX is certainly undefined
        $dbX = undef;
    } else {
	$sth = $self->prepare( "
       SELECT xref_id
         FROM xref
        WHERE external_db_id = ?
          AND dbprimary_acc = ?
          AND version = ?
     " );
	$sth->execute( $dbRef, $exObj->primary_id(), 
		       $exObj->version() );
	( $dbX ) = $sth->fetchrow_array();
    }
    
    if( ! defined $dbX ) {
	
	$sth = $self->prepare( "
      INSERT INTO xref 
       SET dbprimary_acc = ?,
           display_label = ?,
           version = ?,
           description = ?,
           external_db_id = $dbRef
     " );
	$sth->execute( $exObj->primary_id(), $exObj->display_label(), $exObj->version(),
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
     SELECT xref_id,
            synonym
       FROM external_synonym
      WHERE xref_id = '$dbX'
        AND synonym = '$syn'
    " );
	    $sth->execute;
	    
	    my ($dbSyn) = $sth->fetchrow_array();
	    
	    #print STDERR $dbSyn[0],"\n";
	    
	    if( ! $dbSyn ) {
		$sth = $self->prepare( "
        INSERT INTO external_synonym
         SET xref_id = $dbX,
            synonym = '$syn'
      " );
		$sth->execute();
	    }
	}
	
	
	$sth = $self->prepare( "
   INSERT INTO object_xref
     SET xref_id = $dbX,
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
             INSERT INTO identity_xref
             SET object_xref_id = $Xidt,
             query_identity = ?,
             target_identity = ?
    " );
	    $sth->execute( $exObj->query_identity, $exObj->target_identity );
	    
	}
    } else {
	$sth = $self->prepare ( "

              SELECT xref_id
              FROM object_xref
              WHERE xref_id = $dbX
              AND   ensembl_object_type = '$ensType'
              AND   ensembl_id = '$ensObject'");
	
	$sth->execute;
	my ($tst) = $sth->fetchrow_array;


	if (! defined $tst) {
	# line is already in xref table. Need to add to object_xref
	    $sth = $self->prepare( "
             INSERT INTO object_xref
               SET xref_id = $dbX,
               ensembl_object_type = ?,
               ensembl_id = ?");
	
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
             INSERT INTO identity_xref
             SET object_xref_id = $Xidt,
             query_identity = ?,
             target_identity = ?
    " );
		$sth->execute( $exObj->query_identity, $exObj->target_identity );
		
	    }


	}
    }
        
    return $dbX;
    
}


sub fetch_by_gene {
  my ( $self, $gene ) = @_;
  my $query1 = "SELECT t.translation_id 
                FROM transcript t
                WHERE t.gene_id = ?";

  my $sth1 = $self->prepare($query1);
  $sth1->execute( $gene->dbID );

  while (my $transid = $sth1->fetchrow) {

    my @translation_xrefs = $self->_fetch_by_EnsObject_type( $transid, 'Translation' );
    foreach my $translink(@translation_xrefs) {
      $gene->add_DBLink($translink);
    }
  }
  my @genelinks = $self->_fetch_by_EnsObject_type( $gene->stable_id, 'Gene' );
  foreach my $genelink ( @genelinks ) {
    $gene->add_DBLink( $genelink );
  }
}

sub fetch_by_rawContig {
  my ( $self, $rawContigId ) = @_;
  return $self->_fetch_by_EnsObject_type( $rawContigId, 'RawContig' );
}

sub fetch_by_transcript {
  my ( $self, $trans ) = @_;

  my $query1 = "SELECT t.translation_id 
                FROM transcript t
                WHERE t.transcript_id = ?";

  my $sth1 = $self->prepare($query1);
  $sth1->execute( $trans->dbID );

  # 
  # Did this to be consistent with fetch_by_Gene, but don't like
  # it (filling in the object). I think returning the array would
  # be better. Oh well. EB
  #
  
  while (my $transid = $sth1->fetchrow) {
      my @translation_xrefs = $self->_fetch_by_EnsObject_type( $transid, 'Translation' );
      foreach my $translink(@translation_xrefs) {
	  $trans->add_DBLink($translink);
      }
  }


}

sub fetch_by_translation {
  my ( $self, $trslId ) = @_;
  return $self->_fetch_by_EnsObject_type( $trslId, 'Translation' );
}


sub _fetch_by_EnsObject_type {
  my ( $self, $ensObj, $ensType ) = @_;
  my @out;
  
  my $sth = $self->prepare("
    SELECT xref.xref_id, xref.dbprimary_acc, xref.display_label,
          xref.version, xref.description,
          exDB.db_name, exDB.release, oxr.object_xref_id, es.synonym, idt.query_identity, idt.target_identity
    FROM xref, external_db exDB, object_xref oxr LEFT JOIN external_synonym es on es.xref_id = xref.xref_id 
                                               LEFT JOIN identity_xref idt on idt.object_xref_id = oxr.object_xref_id
    WHERE xref.xref_id = oxr.xref_id
      AND xref.external_db_id = exDB.external_db_id 
      AND oxr.ensembl_id = '$ensObj'
      AND oxr.ensembl_object_type = '$ensType'
  ");
  
  $sth->execute();
  
  
  my %seen;
  
  while ( my $arrRef = $sth->fetchrow_arrayref() ) {
    my ( $refID, $dbprimaryId, $displayid, $version, $desc, $dbname, $release, $objid, 
         $synonym, $queryid, $targetid ) =
      @$arrRef;
    
    my $exDB;
    
    # using an outer join on the synonyms as well as on identity_xref, we
    # now have to filter out the duplicates (see v.1.18 for
    # original). Since there is at most one identity_xref row per xref,
    # this is easy enough; all the 'extra' bits are synonyms
    if ( !$seen{$refID} )  {
      $seen{$refID}++;
      
      if ((defined $queryid)) {         # an xref with similarity scores
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
        
      } else {
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
      push( @out, $exDB );
    }                                   # if (!$seen{$refID})

    # $exDB still points to the same xref, so we can keep adding synonyms
    #if ($synonym) {
    #  $exDB->add_synonym( $synonym );
    #}
  }                                     # while <a row from database>
  
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

   my $sth = $self->prepare("SELECT DISTINCT( tr.gene_id ) 
                  FROM transcript tr, 
                       xref x, object_xref oxr
                  WHERE tr.translation_id = oxr.ensembl_id 
                    AND oxr.xref_id = x.xref_id 
                    AND x.display_label = '$name'");
   $sth->execute();

   while( ($a) = $sth->fetchrow_array ) {
       push(@genes,$a);
   }

   unless (scalar @genes){ 
       $sth = $self->prepare("SELECT DISTINCT( tr.gene_id ) 
		      FROM transcript tr, 
			   Xref x, objectXref oxr
		      WHERE tr.translation_id = oxr.ensembl_id 
			AND oxr.xrefId = x.xrefId 
			AND x.dbprimary_id='$name'");
       $sth->execute();
       while( ($a) = $sth->fetchrow_array ) {
	   push(@genes,$a);
       }
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
    FROM xref, external_db exDB, object_xref oxr, external_synonym syn 
     WHERE (xref.dbprimary_acc = '$name'
            AND xref.xref_id = oxr.xref_id
            AND oxr.ensembl_object_type = '$ensType')
     OR    (xref.display_label = '$name'
            AND xref.xref_id = oxr.xref_id
            AND oxr.ensembl_object_type = '$ensType')
     OR    (syn.synonym = '$name'
            AND syn.xref_id = oxr.xref_id
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

  my $sth = $self->prepare( "drop table if exists object_xref, xref, externalDescription, external_synonym, external_db" );
  $sth->execute();

  $sth = $self->prepare( qq{
    
    CREATE TABLE object_xref(
			     object_xref_id INT not null auto_increment,
			     ensembl_id int unsigned not null, 
			     ensembl_object_type ENUM( 'RawContig', 'Transcript', 'Gene', 'Translation' ) not null,
			     xref_id INT unsigned not null,
			     
			     UNIQUE ( ensembl_object_type, ensembl_id, xref_id ),
			     KEY xref_index( object_xref_id, xref_id, ensembl_object_type, ensembl_id )
			    );
  } );
  $sth->execute();
  $sth = $self->prepare( qq{
    CREATE TABLE xref (
		       xref_id INT unsigned not null auto_increment,
		       external_db_id int not null,
		       dbprimary_acc VARCHAR(40) not null,
		       display_label VARCHAR(40) not null,
		       version VARCHAR(10) DEFAULT '' NOT NULL,
		       description VARCHAR(255),
		       
		       PRIMARY KEY( xref_id ),
		       UNIQUE KEY id_index( dbprimary_acc, external_db_id ),
		       KEY display_index ( display_label )
		      );
    
   } );

  $sth->execute();

  $sth = $self->prepare( qq{
     CREATE TABLE external_synonym(
         xref_id INT not null,
         synonym VARCHAR(40) not null,
         PRIMARY KEY( xref_id, synonym ),
	 KEY nameIdx( synonym )) 
   } );
  $sth->execute();

  $sth = $self->prepare( qq{
     CREATE TABLE external_db(
         external_db_id INT not null auto_increment,
         db_name VARCHAR(40) not null,
	 release VARCHAR(40),
         PRIMARY KEY( external_db_id ) ) 
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

Objectxref
=============
ensembl_id varchar, later int
ensembl_object_type  enum 
xref_id int
primary key (ensembl_id,ensembl_object_type,xref_id) 


xref
=================
xref_id int (autogenerated) 
external_db_id int
dbprimary_acc  varchar
version varchar

primary key (xref_id)

ExternalDescription
=======================
xref_id int
description varchar (256)

primary key (xref_id)

ExternalSynonym
=================
xref_id int
synonym varchar

primary key (external_id,synonym)


ExternalDB
===================
external_db_id int
db_name varchar
release varchar

