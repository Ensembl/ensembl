# EnsEMBL Transcript reading writing adaptor for mySQL
#
# Copyright EMBL-EBI 2001
#
# Author: Arne Stabenau
# based on 
# Elia Stupkas Gene_Obj
# 
# Date : 20.02.2001
#

=head1 NAME

Bio::EnsEMBL::DBSQL::TranscriptAdaptor - MySQL Database queries to generate and store transcripts/translations.

=head1 SYNOPSIS

Transcripts and Translations are stored and fetched in this
object. Translations never go alone any more. The database only
accepts them (at the moment) in a transcript.  

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Elia Stupka  : elia@ebi.ac.uk
  Ewan Birney  : 

=cut

package Bio::EnsEMBL::DBSQL::TranscriptAdaptor;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;

@ISA = qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );



=head2 fetch_by_dbID

  Arg [1]    : int $transid
               The unique database id for the transcript to be retrieved.
  Example    : $transcript = $transcript_adaptor->fetch_by_dbID(1234);
  Description: Retreives a transcript from the database via its dbID
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general

=cut
    
sub fetch_by_dbID {
    my ($self, $transid) = @_;

    my $seen = 0;
    my $trans = Bio::EnsEMBL::Transcript->new();
    my $exonAdaptor = $self->db->get_ExonAdaptor();

    my $sth = $self->prepare("SELECT exon_id 
                              FROM   exon_transcript 
                              WHERE  transcript_id = ?
                              ORDER BY rank");
    $sth->execute($transid);

    while( my $rowhash = $sth->fetchrow_hashref) {
	my $exon = $exonAdaptor->fetch_by_dbID($rowhash->{'exon_id'});
	$trans->add_Exon($exon);
	$seen = 1;
    }

    if ($seen == 0 ) {
	$self->throw("transcript $transid is not present in db");
    }
    $trans->dbID($transid);
    $trans->adaptor( $self );

    my $ts = $self->prepare("SELECT translation_id translation_id 
                             FROM   transcript 
                             WHERE  transcript_id = ?");
    $ts->execute($transid);
    my ($val) = $ts->fetchrow_array();
    $trans->_translation_id($val);

    return $trans;
}


=head2 fetch_by_stable_id

  Arg [1]    : string $stable_id 
               The stable id of the transcript to retrieve
  Example    : $trans = $trans_adptr->fetch_by_stable_id('ENST00000309301');
  Description: Retrieves a transcript via its stable id
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_stable_id {
  my ($self, $stable_id) = @_;
  my $sth = $self->prepare( "SELECT transcript_id 
                             FROM   transcript_stable_id 
                             WHERE  stable_id = ?" );
  $sth->execute( $stable_id );

  if( my $arr = $sth->fetchrow_arrayref ) {
    my $transcript = $self->fetch_by_dbID( $arr->[0] );
    return $transcript;
  } else {
    $ self->warn( "No Transcript with this stable id found in the database." );
    return undef;
  }
}



=head2 fetch_by_translation_stable_id

  Arg [1]    : string $transl_stable_id
               The stable identifier of the translation of the transcript to 
               retrieve
  Example    : $t = $tadptr->fetch_by_translation_stable_id('ENSP00000311007');
  Description: Retrieves a Transcript object using the stable identifier of
               its translation.
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_translation_stable_id {
  my ($self, $transl_stable_id) = @_;

  my $sth = $self->prepare( "SELECT t.transcript_id
                             FROM   translation_stable_id tsi, transcript t
                             WHERE  tsi.stable_id = ? 
                             AND    t.translation_id = tsi.translation_id");

  $sth->execute($transl_stable_id);

  my ($id) = $sth->fetchrow_array;
  if ($id){
  	return $self->fetch_by_dbID($id);
  } else {
  	return undef;
  }
} 
=head2 fetch_by_translation_id

=cut

sub fetch_by_translation_id {
    my ($self, $transl_id) = @_;
    my $sth = $self->prepare( 
        "SELECT transcript_id
        FROM transcript
        WHERE translation_id = ?");
    $sth->execute($transl_id);
    my ($id) = $sth->fetchrow_array;
    return $id? $self->fetch_by_dbID($id) : undef;
}

=head2 fetch_all_by_DBEntry

  Arg [1]    : in $external_id
                the external identifier for the transcript to be obtained
  Example    : @trans = @{$trans_adaptor->fetch_all_by_DBEntry($ext_id)}
  Description: retrieves a list of transcripts with an external database
                idenitifier $external_id
  Returntype : listref of Bio::EnsEMBL::DBSQL::Transcript in contig coordinates
  Exceptions : none
  Caller        : ?

=cut

sub fetch_all_by_DBEntry {
  my $self = shift;
  my $external_id = shift;
  my @trans = ();

  my $entryAdaptor = $self->db->get_DBEntryAdaptor();
  my @ids = $entryAdaptor->transcriptids_by_extids($external_id);
  foreach my $trans_id ( @ids ) {
    my $trans = $self->fetch_by_dbID( $trans_id );
    if( $trans ) {
        push( @trans, $trans );
    }
  }
  return \@trans;
}


=head2 fetch_all_by_exon_stable_id

  Arg [1]    : string $stable_id 
               The stable id of an exon in a transcript
  Example    : $trans = $trans_adptr->fetch_all_by_exon_stable_id('ENSE00000309301');
  Description: Retrieves a list of transcripts via an exon stable id
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general

=cut

sub fetch_all_by_exon_stable_id {
  my ($self, $stable_id) = @_;
  my @trans ;
  my $sth = $self->prepare( qq(	SELECT et.transcript_id 
  								FROM exon_transcript as et, 
									 exon_stable_id as esi 
								WHERE esi.exon_id = et.exon_id and 
									  esi.stable_id = "$stable_id"  ));
  $sth->execute(  );

  while( my $id = $sth->fetchrow_array ) {    
		my $transcript = $self->fetch_by_dbID( $id );
    	push(@trans, $transcript) if $transcript;

  } 
  if (!@trans) {
    $self->warn( "No Transcript with this exon stable id found in the database." );
    return undef;
  }
  return \@trans;
}


=head2 store

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript
               The transcript to be written to the database
  Arg [2]    : int $gene_dbID
               The identifier of the gene that this transcript is associated 
               with
  Example    : $transID = $transcriptAdaptor->store($transcript, $gene->dbID);
  Description: Stores a transcript in the database and returns the new
               internal identifier for the stored transcript.
  Returntype : int 
  Exceptions : none
  Caller     : general

=cut

sub store {
   my ($self,$transcript,$gene_dbID) = @_;

   if( ! ref $transcript || !$transcript->isa('Bio::EnsEMBL::Transcript') ) {
       $self->throw("$transcript is not a EnsEMBL transcript - not dumping!");
   }

   # store translation
   # then store the transcript
   # then store the exon_transcript table

   my $translation = $transcript->translation();
   my ( $exon_count, $exons );
   $exons = $transcript->get_all_Exons();
   $exon_count = scalar( @{$exons} );

   my $exonAdaptor = $self->db->get_ExonAdaptor();
   foreach my $exon ( @{$exons} ) {
     $exonAdaptor->store( $exon );
   }

   if( defined $translation ) {
     $self->db->get_TranslationAdaptor()->store( $translation );
   }


   # first store the transcript w/o a display xref
   # the display xref needs to be set after xrefs are stored which needs to 
   # happen after transcript is stored...

   my $xref_id = 0;

   # ok - now load this line in
   my $tst = $self->prepare("
        insert into transcript ( gene_id, translation_id, 
                                 exon_count, display_xref_id )
        values ( ?, ?, ?, 0)
        ");

   if( defined $translation ) {
     $tst->execute( $gene_dbID, $translation->dbID(), $exon_count );
   } else {
     $tst->execute( $gene_dbID, 0, $exon_count );
   }

   my $transc_dbID = $tst->{'mysql_insertid'};

   #store the xrefs/object xref mapping
   my $dbEntryAdaptor = $self->db->get_DBEntryAdaptor();

   foreach my $dbe ( @{$transcript->get_all_DBEntries} ) {
     $dbEntryAdaptor->store( $dbe, $transc_dbID, "Transcript" );
   }

   #
   # Update transcript to point to display xref if it is set 
   #
   if(my $dxref = $transcript->display_xref) {
     if(my $dxref_id = $dbEntryAdaptor->exists($dxref)) {
       my $sth = $self->prepare( "update transcript set display_xref_id = ?".
                                 " where transcript_id = ?");
       $sth->execute($dxref_id, $transc_dbID);
       $dxref->dbID($dxref_id);
       $dxref->adaptor($dbEntryAdaptor);
     }
   }

   my $etst = 
     $self->prepare("insert into exon_transcript (exon_id,transcript_id,rank)"
                    ." values (?,?,?)");
   my $rank = 1;
   foreach my $exon ( @{$transcript->get_all_Exons} ) {
     $etst->execute($exon->dbID,$transc_dbID,$rank);
     $rank++;
   }

   if (defined($transcript->stable_id)) {
     if (!defined($transcript->version)) {
       $self->throw("Trying to store incomplete stable id information for " .
                    "transcript");
     }

     my $statement = 
       "INSERT INTO transcript_stable_id(transcript_id,stable_id,version)" .
         " VALUES(?, ?, ?)";
     my $sth = $self->prepare($statement);
     $sth->execute($transc_dbID, $transcript->stable_id, $transcript->version);
   }

   $transcript->dbID( $transc_dbID );
   $transcript->adaptor( $self );
   return $transc_dbID;
}

=head2 get_stable_entry_info

 Title   : get_stable_entry_info
 Usage   : $transcriptAdptor->get_stable_entry_info($transcript)
 Function: gets stable info for gene and places it into the hash
 Returns : 
 Args    : 


=cut

sub get_stable_entry_info {
  my ($self,$transcript) = @_;

  unless( defined $transcript && ref $transcript && 
	  $transcript->isa('Bio::EnsEMBL::Transcript') ) {
    $self->throw("Needs a Transcript object, not a $transcript");
  }

  my $sth = $self->prepare("SELECT stable_id, version 
                            FROM   transcript_stable_id 
                            WHERE  transcript_id = ?");
  $sth->execute($transcript->dbID());

  my @array = $sth->fetchrow_array();
  $transcript->{'_stable_id'} = $array[0];
  $transcript->{'_version'}   = $array[1];
  

  return 1;
}

=head2 get_Interpro_by_transid

  Arg [1]    : string $trans
               the stable if of the trans to obtain
  Example    : @i = $trans_adaptor->get_Interpro_by_transid($trans->stable_id()); 
  Description: gets interpro accession numbers by transcript stable id.
               A hack really - we should have a much more structured 
               system than this
  Returntype : listref of strings 
  Exceptions : none 
  Caller     : domainview? , GeneView

=cut

sub get_Interpro_by_transid {
   my ($self,$transid) = @_;
   my $sql="
	SELECT	i.interpro_ac, 
		x.description 
        FROM	transcript_stable_id tsi,
                transcript t, 
		protein_feature pf, 
		interpro i
     LEFT JOIN  xref x
            ON  x.dbprimary_acc = i.interpro_ac
	WHERE	tsi.stable_id = '$transid' 
	    AND	t.transcript_id = tsi.transcript_id
	    AND	t.translation_id = pf.translation_id 
	    AND	i.id = pf.hit_id";
   
   my $sth = $self->prepare($sql);
   $sth->execute;

   my @out;
   my %h;
   while( (my $arr = $sth->fetchrow_arrayref()) ) {
       if( $h{$arr->[0]} ) { next; }
       $h{$arr->[0]}=1;
       my $string = $arr->[0] .":".$arr->[1];
       push(@out,$string);
   }


   return \@out;
}




sub remove {
  my $self = shift;
  my $transcript = shift;
  my $gene = shift;

  if( ! defined $transcript->dbID() ) {
    return;
  }

  my $exonAdaptor = $self->db->get_ExonAdaptor();
  my $translationAdaptor = $self->db->get_TranslationAdaptor();

  if( defined $transcript->translation ) {
    $translationAdaptor->remove( $transcript->translation );
  }

  my $sth = $self->prepare( "DELETE FROM exon_transcript 
                             WHERE transcript_id = ?" );
  $sth->execute( $transcript->dbID );
  $sth = $self->prepare( "DELETE FROM transcript_stable_id 
                          WHERE transcript_id = ?" );
  $sth->execute( $transcript->dbID );
  $sth = $self->prepare( "DELETE FROM transcript 
                          WHERE transcript_id = ?" );
  $sth->execute( $transcript->dbID );
  
  foreach my $exon ( @{$transcript->get_all_Exons()} ) {  
    my $sth = $self->prepare( "SELECT count(*) 
                               FROM   exon_transcript 
                               WHERE  exon_id = ?" );
    $sth->execute( $exon->dbID );
    my ($count) = $sth->fetchrow_array;
    if($count == 0){ 
      $exonAdaptor->remove( $exon );
    } else{
      $self->warn("exon " . $exon->dbID . " is not exclusive to transcript " . 
		  $transcript->dbID . "\n");
    }

  }
  
  $transcript->{'dbID'} = undef;
}





=head2 get_display_xref

  Arg [1]    : int $dbID
               the database identifier of the transcript for which the name 
               of external db from which its external name is derived.
  Example    : $external_dbname = $transcript_adaptor->get_display_xref_id(42);
  Description: Retrieves the display_xref_id for a transcript.
  Returntype : int
  Exceptions : thrown if $dbId arg is not defined
  Caller     : general

=cut

sub get_display_xref {
  my ($self, $transcript) = @_;

  if( !defined $transcript ) {
      $self->throw("Must call with a Transcript object");
  }

  my $sth = $self->prepare("SELECT e.db_name,
                                   x.display_label,
                                   x.xref_id
                            FROM   transcript t, 
                                   xref x, 
                                   external_db e
                            WHERE  t.transcript_id = ?
                              AND  t.display_xref_id = x.xref_id
                              AND  x.external_db_id = e.external_db_id
                           ");
  $sth->execute( $transcript->dbID );


  my ($db_name, $display_label, $xref_id ) = $sth->fetchrow_array();
  if( !defined $xref_id ) {
    return undef;
  }

  my $db_entry = Bio::EnsEMBL::DBEntry->new
    (
     -dbid => $xref_id,
     -adaptor => $self->db->get_DBEntryAdaptor(),
     -dbname => $db_name,
     -display_id => $display_label
    );

  return $db_entry;
}


=head2 update

  Arg [1]    : Bio::EnsEMBL::Transcript
  Example    : $transcript_adaptor->update($transcript);
  Description: Updates a transcript in the database
  Returntype : None
  Exceptions : thrown if the $transcript is not a Bio::EnsEMBL::Transcript
               warn if trying to update the number of attached exons.  This
               is a far more complex process and is not yet implemented.
               warn if the method is called on a transcript that does not exist 
               in the database.
  Caller     : general

=cut

sub update {
   my ($self,$transcript) = @_;
   my $update = 0;

   if( !defined $transcript || !ref $transcript || !$transcript->isa('Bio::EnsEMBL::Transcript') ) {
       $self->throw("Must update a transcript object, not a $transcript");
   }

   my $sth = $self->prepare("UPDATE transcript
                          SET    display_xref_id = ?,
                                 exon_count = ?
                          WHERE  transcript_id = ?
                         ");

   my $display_xref = $transcript->display_xref();
   my $display_xref_id;

   if( defined $display_xref && $display_xref->dbID() ) {
     $display_xref_id = $display_xref->dbID();
   } else {
     $display_xref_id = 0;
   }

   my $exon_count = scalar( @{$transcript->get_all_Exons()} );
   $sth->execute( $display_xref_id, $exon_count, $transcript->dbID() );
 }

=head2 list_dbIDs

  Arg [1]    : none
  Example    : @transcript_ids = @{$transcript_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all transcripts in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?

=cut

sub list_dbIDs {
   my ($self) = @_;

   return $self->_list_dbIDs("transcript");
}

=head2 list_stable_dbIDs

  Arg [1]    : none
  Example    : @stable_transcript_ids = @{$transcript_adaptor->list_stable_dbIDs()};
  Description: Gets an array of stable ids for all transcripts in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?

=cut

sub list_stable_ids {
   my ($self) = @_;

   return $self->_list_dbIDs("transcript_stable_id", "stable_id");
}



1;

