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

;

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

 Title   : fetch_by_dbID
 Usage   : $transcriptobj->fetch_by_dbID( $dbid )
 Function: 
 Example : $obj->get( .. )
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : gene id and supporting tag if latter not specified, assumes without
	   Note that it is much faster to get genes without supp.evidence!

=cut

    
sub fetch_by_dbID {
    my ($self,$transid) = @_;

    my $seen = 0;
    my $trans = Bio::EnsEMBL::Transcript->new();
    my $exonAdaptor = $self->db->get_ExonAdaptor();

    my $sth = $self->prepare("select exon_id from exon_transcript where transcript_id = $transid ORDER BY rank");
    $sth->execute();

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

    my $ts = $self->prepare("select translation_id from transcript where transcript_id = $transid");
    $ts->execute;
    my ($val) = $ts->fetchrow_array();
    $trans->_translation_id($val);

    return $trans;
}


sub fetch_by_stable_id {
  my ($self, $stable_id) = @_;
  my $sth = $self->prepare( "select transcript_id from transcript_stable_id where stable_id = ?" );
  $sth->execute( $stable_id );

  if( my $arr = $sth->fetchrow_arrayref ) {
    my $transcript = $self->fetch_by_dbID( $arr->[0] );
    return $transcript;
  } else {
    $ self->warn( "No Transcript with this stable id found in the database." );
    return undef;
  }
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
   my $exonAdaptor = $self->db->get_ExonAdaptor();

   if( ! ref $transcript || !$transcript->isa('Bio::EnsEMBL::Transcript') ) {
       $self->throw("$transcript is not a EnsEMBL transcript - not dumping!");
   }

   if( ! ref $gene || !$gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("$gene is not a EnsEMBL gene - not dumping!");
   }

   # store exons
   # store translation
   # then store the transcript
   # then store the exon_transcript table

   foreach my $exon ( $transcript->get_all_Exons() ) {
     $exonAdaptor->store( $exon );
   }

   my $translation = $transcript->translation();
   if( defined $translation ) {
     $self->db->get_TranslationAdaptor()->store( $translation );
   }
   # ok - now load this line in
   my $tst = $self->prepare("
        insert into transcript ( gene_id, translation_id )
        values ( ?, ?)
        ");

   if( defined $translation ) {
     $tst->execute( $gene->dbID(), $translation->dbID() );
   } else {
     $tst->execute( $gene->dbID(), 0 );
   }

   $transcript->dbID( $tst->{'mysql_insertid'});
   $transcript->adaptor( $self );

   #print STDERR "Going to look at gene links\n";
   my $dbEntryAdaptor = $self->db->get_DBEntryAdaptor();

   foreach my $dbl ( $transcript->each_DBLink ) {
     $dbEntryAdaptor->store( $dbl, $transcript->dbID, "Transcript" );
   }

   my $etst = $self->prepare("insert into exon_transcript (exon_id,transcript_id,rank) values (?,?,?)");
   my $rank = 1;
   foreach my $exon ( $transcript->get_all_Exons ) {
     $etst->execute($exon->dbID,$transcript->dbID,$rank);
     $rank++;
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

  if( !defined $transcript || !ref $transcript || !$transcript->isa('Bio::EnsEMBL::Transcript') ) {
     $self->throw("Needs a Transcript object, not a $transcript");
  }

  my $sth = $self->prepare("select stable_id,version from transcript_stable_id where transcript_id = ".$transcript->dbID);
  $sth->execute();

  my @array = $sth->fetchrow_array();
  $transcript->{'_stable_id'} = $array[0];
  $transcript->{'_version'}   = $array[1];
  

  return 1;
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

  my $sth = $self->prepare( "delete from exon_transcript where transcript_id = ?" );
  $sth->execute( $transcript->dbID );
  $sth = $self->prepare( "delete from transcript_stable_id where transcript_id = ?" );
  $sth->execute( $transcript->dbID );
  $sth = $self->prepare( "delete from transcript where transcript_id = ?" );
  $sth->execute( $transcript->dbID );

  foreach my $exon ( $transcript->get_all_Exons() ) {
    $exonAdaptor->remove( $exon );
  }
  

  $transcript->{'dbID'} = undef;
}

1;

