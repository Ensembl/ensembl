# EnsEMBL Translation reading writing adaptor for mySQL
#
# Copyright EMBL-EBI 2001
#
# Author: Arne Stabenau
# 
# Date : 21.07.2001
#

=head1 NAME

Bio::EnsEMBL::DBSQL::TranslationAdaptor - Provides a means to fetch and store
Translation objects from a database.

=head1 SYNOPSIS

  my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(...);

  my $transcript_adaptor = $db->get_TranscriptAdaptor();
  my $translation_adaptor = $db->get_TranslationAdaptor();
  my $transcript = $transcript_adaptor->fetch_by_dbID(131243);
  my $translation = $translation_adaptor->fetch_by_Transcript($transcript);

  print "Translation Start Site: " .
        $translation->start_Exon()->stable_id(), " ", $translation->start();
  print "Translation Stop: " . 
        $translation->end_Exon()->stable_id(), " ", $translation->end();


=head1 CONTACT

  ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::TranslationAdaptor;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate );


@ISA = qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );



=head2 fetch_by_Transcript

  Arg [1]    : none, string, int, Bio::EnsEMBL::Example $formal_parameter_name
    Additional description lines
    list, listref, hashref
  Example    :  ( optional )
  Description: testable description
  Returntype : none, txt, int, float, Bio::EnsEMBL::Example
  Exceptions : none
  Caller     : object::methodname or just methodname

=cut

sub fetch_by_Transcript {
  my ( $self, $transcript ) = @_;

  if(!ref($transcript) || !$transcript->isa('Bio::EnsEMBL::Transcript')) {
    throw("Bio::EnsEMBL::Transcript argument is required.");
  }

  my $sql = "
    SELECT tl.translation_id, tl.start_exon_id,
           tl.end_exon_id, tl.seq_start, tl.seq_end,
           tlsi.stable_id, tlsi.version
      FROM translation tl
 LEFT JOIN translation_stable_id tlsi
        ON tlsi.translation_id = tl.translation_id
     WHERE tl.transcript_id = ?";

  my $transcript_id = $transcript->dbID();
  my $sth = $self->prepare( $sql );
  $sth->execute( $transcript_id );

  my ( $translation_id, $start_exon_id, $end_exon_id,
       $seq_start, $seq_end, $stable_id, $version ) = 
  $sth->fetchrow_array();

  if( ! defined $translation_id ) {
    return undef;
  }

  my ($start_exon, $end_exon);

  # this will load all the exons whenever we load the translation
  # but I guess thats ok ....

  foreach my $exon (@{$transcript->get_all_Exons()}) {
    if($exon->dbID() == $start_exon_id ) {
      $start_exon = $exon;
    }

    if($exon->dbID() == $end_exon_id ) {
      $end_exon = $exon;
    }
  }

  unless($start_exon && $end_exon) {
     throw("Could not find start or end exon in transcript\n");
  }

  my $translation = Bio::EnsEMBL::Translation->new
   (
     -dbID => $translation_id,
     -adaptor => $self,
     -seq_start => $seq_start,
     -seq_end => $seq_end,
     -start_exon => $start_exon,
     -end_exon => $end_exon,
     -stable_id => $stable_id,
     -version => $version
   );

  return $translation;
}




=head2 fetch_all_by_external_name

  Arg [1]    : string $external_id
               The external identifier for the tranlsation(s) to be obtained.
  Example    : @tls = @{$trl_adaptor->fetch_all_by_external_name('BRCA2')};
  Description: Retrieves a list of translations fetched via an external
               identifier.  Note that this may not be a particularly useful
               method, because translations do not make much sense out of the 
               context of their transcript.  It may be better to use the
               TranscriptAdaptor::fetch_all_by_external_name instead.
  Returntype : reference to a list of Translations
  Exceptions : none
  Caller     : general

=cut

sub fetch_all_by_external_name {
  my $self = shift;
  my $external_id = shift;

  my $entry_adaptor = $self->db->get_DBEntryAdaptor();
  my @ids = $entry_adaptor->list_translation_ids_by_extids($external_id);

  my $transcript_adaptor = $self->db()->get_TranscriptAdaptor();

  my @out;

  foreach my $id (@ids) {
    my $transcript = $transcript_adaptor->fetch_by_translation_id($id);
    if($transcript) {
      push @out, $self->fetch_by_Transcript($transcript);
    }
  } 

  return \@out;
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
                   
  my $sth = $self->prepare( "
         INSERT INTO translation( seq_start, start_exon_id,
                                  seq_end, end_exon_id, transcript_id) 
         VALUES ( ?,?,?,?,? )");

  $sth->execute( $translation->start(),
                 $translation->start_Exon()->dbID(),
                 $translation->end(),
                 $translation->end_Exon()->dbID(),
                 $transcript_id );
  
  my $transl_dbID = $sth->{'mysql_insertid'};

  #store object xref mappings to translations

  my $dbEntryAdaptor = $self->db()->get_DBEntryAdaptor();
  #store each of the xrefs for this translation
  foreach my $dbl ( @{$translation->get_all_DBEntries} ) {
     $dbEntryAdaptor->store( $dbl, $transl_dbID, "Translation" );
  }
  
  
  if (defined($translation->stable_id)) {
    if (!defined($translation->version)) {
     throw("Trying to store incomplete stable id information for translation");
    }
    
    my $sth = $self->prepare
      ("INSERT INTO translation_stable_id(translation_id, stable_id, version)".
       "     VALUES (?, ?, ?)");
   
    $sth->execute($transl_dbID, $translation->stable_id(), 
                  $translation->version());
    
    $sth->finish();
  }

  $translation->dbID( $transl_dbID );
  $translation->adaptor( $self );

  return $transl_dbID;
}




sub remove {
  my $self = shift;
  my $translation = shift;

  my $sth = $self->prepare("DELETE FROM translation 
                            WHERE translation_id = ?" );
  $sth->execute( $translation->dbID );
  $sth = $self->prepare("DELETE FROM translation_stable_id 
                         WHERE translation_id = ?" );
  $sth->execute( $translation->dbID );
  $translation->dbID( undef );
  $translation->adaptor(undef);
}


=head2 list_dbIDs

  Arg [1]    : none
  Example    : @translation_ids = @{$translation_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all translations in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?

=cut

sub list_dbIDs {
   my ($self) = @_;

   return $self->_list_dbIDs("translation");
}


=head2 list_stable_dbIDs

  Arg [1]    : none
  Example    : @transl_stable_ids = @{$transl_adaptor->list_stable_dbIDs()};
  Description: Gets an array of stable ids for all translations in the current 
               db
  Returntype : reference to a list of strings
  Exceptions : none
  Caller     : general

=cut

sub list_stable_ids {
   my ($self) = @_;

   return $self->_list_dbIDs("translation_stable_id", "stable_id");
}



=head2 fetch_by_dbID

  Arg [1]    : int $dbID
               The internal identifier of the Translation to obtain
  Example    : $translation = $translation_adaptor->fetch_by_dbID(1234);
  Description: This fetches a Translation object via its internal id.  This
               is only debatably useful since translations do not make much
               sense outside of the context of their Translation.  Consider
               using fetch_by_Transcript instead.
  Returntype : Bio::EnsEMBL::Translation or undef if the translation is not
               found.
  Exceptions : warning if an additional (old style) Transcript argument is
               provided
  Caller     : ?

=cut

sub fetch_by_dbID {
   my ($self,$dbID, $transcript) = @_;

   if($transcript) {
     deprecate("Use of fetch_by_dbID with a Transcript argument is deprecated."
               . "Use fetch_by_Transcript instead." );
   }

   if(!$dbID) {
     throw("dbID argument is required");
   }

   my $transcript_adaptor = $self->db()->get_TranscriptAdaptor();
   my $transcript = $transcript_adaptor->fetch_by_translation_id($dbID);

   return undef if(!$transcript);

   return $self->fetch_by_Transcript($transcript);
}


=head2 fetch_by_dbID

  Arg [1]    : int $dbID
               The internal identifier of the Translation to obtain
  Example    : $translation = $translation_adaptor->fetch_by_dbID(1234);
  Description: This fetches a Translation object via its internal id.  This
               is only debatably useful since translations do not make much
               sense outside of the context of their Translation.  Consider
               using fetch_by_Transcript instead.
  Returntype : Bio::EnsEMBL::Translation or undef if the translation is not
               found.
  Exceptions : warning if an additional (old style) Transcript argument is
               provided
  Caller     : ?

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




=head2 fetch_all_by_DBEntry

  Description: DEPRECATED, this has been renames fetch_all_by_external_name

=cut

sub fetch_all_by_DBEntry {
  my $self = shift;
  deprecate("Use fetch_all_by_external_name instead.");
  return $self->fetch_all_by_external_name(@_);
}

=head2 get_stable_entry_info

 Description: DEPRECATED - This method should no longer be needed. Stable
              id info is fetched when the transcript is.

=cut

sub get_stable_entry_info {
  my ($self,$translation) = @_;

  deprecate( "This method shouldnt be necessary any more" );

  unless(defined $translation && ref $translation && 
	 $translation->isa('Bio::EnsEMBL::Translation') ) {
    throw("Needs a Translation object, not a [$translation]");
  }

  my $sth = $self->prepare("SELECT stable_id, version 
                            FROM   translation_stable_id 
                            WHERE  translation_id = ?");
  $sth->execute($translation->dbID());

  my @array = $sth->fetchrow_array();
  $translation->{'_stable_id'} = $array[0];
  $translation->{'_version'}   = $array[1];
  
  return 1;
}

1;
