#
# Ensembl module for Bio::EnsEMBL::Utils::EMBL::GeneWrapper
#
# Cared for by Ewan Birney <ensembl-dev@ebi.ac.uk>
#
# Copyright EMBL / EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Utils::EMBL::TranscriptWrapper

=head1 SYNOPSIS

   foreach $pt ( $slice->get_all_PredictionTranscripts ) {
     push @seq_feats, new Bio::EnsEMBL::Utils::EMBL::TranscriptWrapper($gene);
   }
       

=head1 DESCRIPTION

Allows prediction transcripts to be EMBL dumped as though they were seq 
features (which in EnsEMBL they are not).  Essentially a workaround to force 
ensembl to fit the Bioperl embl dumping code. This wrapper is an ugly 
replacement for some of the old ugly VirtualGene functionality.

=head1 AUTHOR - Graham McVicker

This modules is part of the Ensembl project http://www.ensembl.org

Email ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal 
methods are usually preceded with a _

=cut


use strict;

package Bio::EnsEMBL::Utils::EMBL::TranscriptWrapper;

use Bio::SeqIO::FTHelper;
use Bio::EnsEMBL::Utils::EMBL::Misc qw(features2join_string);
use Bio::EnsEMBL::Root;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Root Bio::EnsEMBL::Exon);


=head2 new

  Arg [1]    : Bio::EnsEMBL::TranscriptI
  Example    : $tw = new Bio::EnsEMBL::Utils::EMBL::TranscriptWrapper($trans)
  Description: Creates a new TranscriptWrapper 'wrapped' around a TranscriptI 
               object toallow it to be embl dumped. 
  Returntype : Bio::EnsEMBL::Utils::EMBL::TranscriptWrapper
  Exceptions : none
  Caller     : Embl dumping scripts

=cut

sub new {
  my ($class, $transcript) = @_;

  my $self = $class->SUPER::new();

  unless($transcript && ref $transcript && $transcript->isa('Bio::EnsEMBL::TranscriptI')) {
    $self->throw("Transcript argument must be a Bio::EnsEMBL::TranscriptI");
  }

  $self->transcript($transcript);
  $self->{_strict_embl_dumping} = 0;

  return $self;
}


=head2 transcript

  Arg [1]    : (optional) Bio::EnsEMBL::Transcript $trans
  Example    : $trans = $transcript_wrapper->transcript();
  Description: getter/setter fot the transcript contained inside this wrapper
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : thrown if the $trans arg is not a Bio::EnsEMBL::TranscriptI
  Caller     : internal

=cut

sub transcript {
  my ($self, $trans) = @_;

  if($trans) {
    unless($trans->isa('Bio::EnsEMBL::TranscriptI')) {
      $self->throw("[$trans] is not a Bio::EnsEMBL::TranscriptI");
    }
    
    $self->{_transcript} = $trans;
  }

  return $self->{_transcript};
}


=head2 to_FTHelper

  Arg [1]    : none
  Example    : none
  Description: Wrapper function - 
               This is where all the embl dumping magic around genes takes 
               place. Most of the embl output for ensembl genes is defined 
               here.
  Returntype : list of Bio::SeqIO::FTHelper
  Exceptions : none
  Caller     : Bio::SeqIO::embl

=cut

sub to_FTHelper {
  my $self = shift;
  my $trans = $self->transcript;
  
  my @out;

  my @dblinks = ();

  if($trans->can('get_all_DBLinks')) {
    @dblinks = @{$trans->get_all_DBLinks()};
  }

  my $loc = features2join_string($trans->get_all_translateable_Exons());

  my $ft = Bio::SeqIO::FTHelper->new();
  $ft->loc($loc);
  $ft->key('CDS');
  
  $ft->add_field('translation', $trans->translate()->seq());
  if($trans->can('translation')) {
    $ft->add_field('cds', $trans->translation->stable_id());
  }
  $ft->add_field('transcript', $trans->stable_id());
  foreach my $dbl (@dblinks) {
    $ft->add_field('db_xref', $dbl->database().":".$dbl->primary_id());
  }
  push(@out, $ft);
  
  foreach my $exon (@{$trans->get_all_translateable_Exons()}) {
    my $ft = Bio::SeqIO::FTHelper->new();

    $ft->loc(features2join_string([$exon]));
    $ft->key("exon");
    
    if( $self->strict_EMBL_dumping()) {
      $ft->add_field('db_xref', 'ENSEMBL:HUMAN-Exon-'.$exon->stable_id());
    } else {
      $ft->add_field('exon_id', $exon->stable_id());
      $ft->add_field('start_phase', $exon->phase());
      $ft->add_field('end_phase', $exon->phase());
    }
    
    push (@out, $ft);
  }
  
  return @out;
}


  
=head2 strict_EMBL_dumping

  Arg [1]    : optional boolean $newval
  Example    : none
  Description: Getter/ setter for strict embl dumping attribute which 
               determines whether embl dumping will be 'strict'. Default is 0.
  Returntype : boolean
  Exceptions : none
  Caller     : to_FTHelper

=cut

sub strict_EMBL_dumping {
  my ($self, $newval) = @_;
  
  if(defined $newval) {
    $self->{_strict_embl_dumping} = $newval;
  }

  return $self->{_strict_embl_dumping};
}
    

=head2 source_tag

  Arg [1]    : none
  Example    : none
  Description: Wrapper function - returns a source tag for the gene
  Returntype : none
  Exceptions : none
  Caller     : Bio::SeqIO::embl

=cut

sub source_tag {
  my ($self, @args) = @_;

  return "genscan";
}


=head2 primary_tag

  Arg [1]    : none
  Example    : none
  Description: Wrapper function returns a primary tag for thie gene
  Returntype : string
  Exceptions : none
  Caller     : Bio::SeqIO::embl

=cut

sub primary_tag {
  my ($self, @args) = @_;

  return "predicted transcript";
}


=head2 has_tag

  Arg [1]    : none
  Example    : none
  Description: Wrapper function - always returns 0
  Returntype : none
  Exceptions : none
  Caller     : Bio::SeqIO::embl

=cut

sub has_tag {
  my $self = shift;

  return 0;
}


=head2 all_tags

  Arg [1]    : none
  Example    : none
  Description: Wrapper function - does nothing
  Returntype : none
  Exceptions : none
  Caller     : Bio::SeqIO::embl

=cut

sub all_tags {
  my $self = shift;

  return;
}


=head2 each_tag_value

  Arg [1]    : none
  Example    : none
  Description: Wrapper function - returns empty list
  Returntype : empty list
  Exceptions : none
  Caller     : Bio::SeqIO::embl

=cut

sub each_tag_value {
  my $self = shift;

  return ();
}



1;
