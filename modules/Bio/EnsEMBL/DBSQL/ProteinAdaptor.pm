
#
# BioPerl module for ProteinAdaptor
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

ProteinAdaptor - DEPRECATED Use TranslationAdaptor or TranscriptAdaptor instead

=head1 DESCRIPTION

DEPRECATED - Do not use this module.  Use the TranslationAdaptor or
TranscriptAdaptor instead. All Protein object functionality has been replaced
on the Translation object so there is no reason to create Protein objects 
anymore.

=head1 CONTACT

  Post questions to the EnsEMBL development list

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::ProteinAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::TranslationAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(deprecate);

use Bio::EnsEMBL::Protein;

@ISA = qw(Bio::EnsEMBL::DBSQL::TranslationAdaptor);


=head2 fetch_by_transcript_id

  Description: DEPRECATED - use TranslationAdaptor::fetch_by_Transcript or
                            Transcript::translation instead

=cut

sub fetch_by_transcript_id{
   my ($self,$transid) = @_;   
   deprecate("Use Transcript::translation instead\n");
   my $tr_adaptor = $self->db->get_TranscriptAdaptor();
   my $tr = $tr_adaptor->fetch_by_dbID($transid);
   my $tl = $tr->translation();
   $tl->adaptor($self);
   return bless $tl, 'Bio::EnsEMBL::Protein';
}


=head2 fetch_by_translation_stable_id

  Description: DEPRECATED - use TranslationAdaptor::fetch_by_stable_id
               or Transcript::translation instead

=cut

sub fetch_by_translation_stable_id {
  my ($self,$transid) = @_;
  deprecate("Use TranslationAdaptor::fetch_by_stable_id or " .
            "Transcript::translation instead.");
  my $tl_adaptor = $self->db->get_TranslationAdaptor();
  my $tl = $tl_adaptor->fetch_by_stable_id($transid);
  $tl->adaptor($self);
  return bless $tl, 'Bio::EnsEMBL::Protein';
}


=head2 fetch_by_translation_id

  Description: DEPRECATED - use TranslationAdaptor::fetch_by_dbID or

=cut

sub fetch_by_translation_id {
   my ($self, $translation_id) = @_;
   deprecate("Use TranslationAdaptor::fetch_by_dbID or " .
            "Transcript::translation instead.");
   
   my $tl_adaptor = $self->db->get_TranslationAdaptor();
   my $tl = $tl_adaptor->fetch_by_dbID($translation_id);
   $tl->adaptor($self);
   return bless $tl, 'Bio::EnsEMBL::Protein';
}


1;


