# EnsEMBL Translation reading writing adaptor for mySQL
#
# Copyright EMBL-EBI 2001
#
# Author: Arne Stabenau
# 
# Date : 21.07.2001
#

=head1 NAME

Bio::EnsEMBL::DBSQL::TranslationAdaptor - MySQL Database queries to generate and store translations.

=head1 SYNOPSIS

Translations are stored and fetched with this
object. 

=head1 CONTACT

  ensembl-dev@ebi.ac.uk


=head1 APPENDIX

=cut

;

package Bio::EnsEMBL::DBSQL::TranslationAdaptor;

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
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_by_dbID {
   my ($self,$dbID,$transcript) = @_;


   if( !defined $transcript ) {
     $self->throw("Translations make no sense outside of their parent Transcript objects. You must retrieve with Transcript parent");
   }

   my $statement = "select translation_id tlid,seq_start,start_exon_id,seq_end,end_exon_id from translation where translation_id = $dbID";
   my $sth     = $self->prepare($statement);
   my $res     = $sth->execute();
   my $rowhash = $sth->fetchrow_hashref;

   if( !defined $rowhash ) {
     # assumme this is a translationless transcript deliberately
     return undef;
   }

   my $out = Bio::EnsEMBL::Translation->new();

   $out->start        ($rowhash->{'seq_start'});
   $out->end          ($rowhash->{'seq_end'});

   my $start_exon = $transcript->get_Exon_by_dbID($rowhash->{'start_exon_id'});
   $out->start_exon($start_exon);
   my $end_exon   = $transcript->get_Exon_by_dbID($rowhash->{'end_exon_id'});
   $out->end_exon($end_exon);

   $out->dbID         ($rowhash->{'tlid'});

   $out->adaptor( $self );

   return $out;
}



sub store {
  my ( $self, $translation )  = @_;
  #print STDERR "storing translation\n";

  if( !defined $translation->start_exon->dbID || !defined $translation->end_exon->dbID ) {
    $self->throw("Attempting to write a translation where the dbIDs to the start and exons are not set. This is most likely to be because you assigned the exons for translation start_exon and translation end_exon to be different in memory objects from your trnascript exons - although it could also be an internal error in the adaptors. For your info the exon memory locations are ".$translation->start_exon." and ".$translation->end_exon());
  }



  my $sth = $self->prepare( "insert into translation( seq_start, start_exon_id, seq_end, end_exon_id) values( ?,?,?,? )");

  $sth->execute( $translation->start(),
		 $translation->start_exon()->dbID(),
		 $translation->end(),
		 $translation->end_exon()->dbID() );
  $translation->dbID( $sth->{'mysql_insertid'} );
  $translation->adaptor( $self );

  return $translation->dbID();
}


=head2 get_stable_entry_info

 Title   : get_stable_entry_info
 Usage   : $translationAdaptor->get_stable_entry_info($translation)
 Function: gets stable info for translation and places it into the hash
 Returns : 
 Args    : 


=cut

sub get_stable_entry_info {
  my ($self,$translation) = @_;

  if( !defined $translation || !ref $translation || !$translation->isa('Bio::EnsEMBL::Translation') ) {
     $self->throw("Needs a Translation object, not a $translation");
  }

  my $sth = $self->prepare("select stable_id,version from translation_stable_id where translation_id = ".$translation->dbID);
  $sth->execute();

  my @array = $sth->fetchrow_array();
  $translation->{'_stable_id'} = $array[0];
  $translation->{'_version'}   = $array[1];
  

  return 1;
}


sub remove {
  my $self = shift;
  my $translation = shift;

  my $sth = $self->prepare( "delete from translation where translation_id = ?" );
  $sth->execute( $translation->dbID );
  $sth = $self->prepare( "delete from translation_stable_id where translation_id = ?" );
  $sth->execute( $translation->dbID );
  $translation->dbID( undef );
}

1;
