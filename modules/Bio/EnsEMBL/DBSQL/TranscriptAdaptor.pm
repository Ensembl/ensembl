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

=head1 APPENDIX

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

  return $self->fetch_by_dbID($id);
} 
                             


=head2 store

 Title   : store
 Usage   : $transcriptAdaptor->store( $transcript )
 Function: writes a particular transcript *but not the exons* into
           the database
 Example :
 Returns : 
 Args    : needs a gene ...


=cut

sub store {
   my ($self,$transcript,$gene) = @_;

   if( ! ref $transcript || !$transcript->isa('Bio::EnsEMBL::Transcript') ) {
       $self->throw("$transcript is not a EnsEMBL transcript - not dumping!");
   }

   if( ! ref $gene || !$gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("$gene is not a EnsEMBL gene - not dumping!");
   }

   # store translation
   # then store the transcript
   # then store the exon_transcript table

   my $translation = $transcript->translation();
   my $exon_count;
   if( defined $translation ) {
     $self->db->get_TranslationAdaptor()->store( $translation );
   }
   $exon_count = scalar(@{$transcript->get_all_Exons()});

   # assuming that the store is used during the Genebuil process, set
   # the display_xref_id to 0.  This ought to get re-set during the protein
   # pipeline run.  This probably update to the gene table has yet to be
   # implemented.
   my $xref_id = 0;

   # ok - now load this line in
   my $tst = $self->prepare("
        insert into transcript ( gene_id, translation_id, exon_count, display_xref_id )
        values ( ?, ?, ?, ?)
        ");

   if( defined $translation ) {
     $tst->execute( $gene->dbID(), $translation->dbID(), $exon_count, $xref_id );
   } else {
     $tst->execute( $gene->dbID(), 0, $exon_count, $xref_id );
   }

   $transcript->dbID( $tst->{'mysql_insertid'});
   $transcript->adaptor( $self );

   #print STDERR "Going to look at gene links\n";
   my $dbEntryAdaptor = $self->db->get_DBEntryAdaptor();

   foreach my $dbl ( @{$transcript->get_all_DBLinks} ) {
     $dbEntryAdaptor->store( $dbl, $transcript->dbID, "Transcript" );
   }

   my $etst = $self->prepare("insert into exon_transcript (exon_id,transcript_id,rank) values (?,?,?)");
   my $rank = 1;
   foreach my $exon ( @{$transcript->get_all_Exons} ) {
     $etst->execute($exon->dbID,$transcript->dbID,$rank);
     $rank++;
   }

   if (defined($transcript->stable_id)) {
     if (!defined($transcript->version)) {
       $self->throw("Trying to store incomplete stable id information for transcript");
     }

     my $statement = "INSERT INTO transcript_stable_id(transcript_id," .
                                   "stable_id,version)".
                      " VALUES(" . $transcript->dbID . "," .
                               "'" . $transcript->stable_id . "'," .
                               $transcript->version . 
                               ")";
     my $sth = $self->prepare($statement);
     $sth->execute();
   }


   return $transcript->dbID;

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
        FROM	transcript t, 
		protein_feature pf, 
		interpro i, 
                xref x,
		transcript_stable_id tsi
	WHERE	tsi.stable_id = '$transid' 
	    AND	t.transcript_id = tsi.transcript_id
	    AND	t.translation_id = pf.translation_id 
	    AND	i.id = pf.hit_id 
	    AND	i.interpro_ac = x.dbprimary_acc";
   
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



=head2 get_external_name

  Arg [1]    : int $dbID
               the database identifier of the transcript whose external name 
               is sought
  Example    : $external_name = $transcript_adaptor->get_external_name(42);
  Description: Retrieves the external name for a transcript.  This is implemented
               by joining across the xref and transcript tables, using the 
               display_xref_id column.
  Returntype : string
  Exceptions : thrown if $dbId arg is not defined
  Caller     : general

=cut

sub get_external_name {
  my ($self, $dbID) = @_;

  if( !defined $dbID ) {
      $self->throw("Must call with a dbID");
  }

  my $sth = $self->prepare("SELECT x.display_label 
                            FROM   transcript t, 
                                   xref x 
                            WHERE  t.transcript_id = ?
                              AND  t.display_xref_id = x.xref_id
                           ");
  $sth->execute($dbID);

  my ($xref) = $sth->fetchrow_array();
  if( !defined $xref ) {
    return undef;
  }

  return $xref;
}


=head2 get_external_dbname

  Arg [1]    : int $dbID
               the database identifier of the transcript for which the name of
               external db from which its external name is derived.
  Example    : $external_dbname = $transcript_adaptor->get_external_dbname(42);
  Description: Retrieves the external db name for a transcript from which its external
               name is derived..  This is implemented by joining across the xref, 
               transcript and external_db tables, using the display_xref_id column.
  Returntype : string
  Exceptions : thrown if $dbId arg is not defined
  Caller     : general

=cut

sub get_external_dbname {
  my ($self, $dbID) = @_;

  if( !defined $dbID ) {
      $self->throw("Must call with a dbID");
  }

  my $sth = $self->prepare("SELECT e.db_name 
                            FROM   transcript t, 
                                   xref x, 
                                   external_db e
                            WHERE  t.transcript_id = ?
                              AND  t.display_xref_id = x.xref_id
                              AND  x.external_db_id = e.external_db_id
                           ");
  $sth->execute($dbID);

  my ($db_name) = $sth->fetchrow_array();
  if( !defined $db_name ) {
    return undef;
  }

  return $db_name;
}


=head2 get_display_xref_id

  Arg [1]    : int $dbID
               the database identifier of the transcript for which the name 
               of external db from which its external name is derived.
  Example    : $external_dbname = $transcript_adaptor->get_display_xref_id(42);
  Description: Retrieves the display_xref_id for a transcript.
  Returntype : int
  Exceptions : thrown if $dbId arg is not defined
  Caller     : general

=cut

sub get_display_xref_id {
  my ($self, $dbID) = @_;

  if( !defined $dbID ) {
      $self->throw("Must call with a dbID");
  }

  my $sth = $self->prepare("SELECT display_xref_id 
                            FROM   transcript 
                            WHERE  transcript_id = ?
                           ");
  $sth->execute($dbID);

  my ($xref_id) = $sth->fetchrow_array();
  if( !defined $xref_id ) {
    return undef;
  }

  return $xref_id;
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

   my $sth = $self->prepare("SELECT exon_count, display_xref_id 
                             FROM   transcript 
                             WHERE  transcript_id = ?
                           ");

   $sth->execute($transcript->dbID);

   if ( my ($exons, $dxref_id) = $sth->fetchrow_array() ){

     if ( $exons != scalar(@{$transcript->get_all_Exons})) {
       $self->warn("Trying to update the number of exons on a stored transcript. Not yet implemented.");
     }

     if ( $dxref_id != $transcript->display_xref ) {

       $sth = $self->prepare("UPDATE transcript
                              SET    display_xref_id = ?
                              WHERE  transcript_id = ?
                            ");

       $sth->execute($transcript->display_xref, $transcript->dbID);
     }
   }
   else {
     $self->warn("Trying to update a transcript that is not in the database.  Try store\'ing first.");
   }
}




1;

