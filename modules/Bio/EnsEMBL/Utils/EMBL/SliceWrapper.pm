#
# Ensembl module for Bio::EnsEMBL::Utils::EMBL::SliceWrapper
#
# Cared for by Ewan Birney <ensembl-dev@ebi.ac.uk>
#
# Copyright EMBL/EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Utils::EMBL::SliceWrapper

=head1 SYNOPSIS

$slice = $slice_adaptor->fetch_by_chr_start_end('1', 1, 1000000);

#wrap the slice so it can be embl dumped
$sw = new Bio::EnsEMBL::Utils::EMBL::SliceWrapper($slice);
$sw->skip_SeqFeature('external', 1);

#add comments
&Bio::EnsEMBL::Utils::EMBL::EMBL_Dump::add_ensembl_comments($sw);       

#dump to file in embl format
$outfile = new Bio::SeqIO( '-format' => 'EMBL', -file => ">1:1-1000000.embl");
$outfile->write_seq( $sw );

=head1 DESCRIPTION

This is a wrapper that allows slices to be treated as AnnSeq objects and
dumped in embl format correctly.  It is basically a hack to force the Slice 
object conform to the complexities of Bioperl EMBL Dumping.

=head1 AUTHOR - Graham McVicker

This modules is part of the Ensembl project http://www.ensembl.org

Email ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


use strict;

package Bio::EnsEMBL::Utils::EMBL::SliceWrapper;

use vars qw(@ISA);

use Bio::EnsEMBL::Root;
use Bio::Annotation::Collection;
use Bio::SeqFeature::Generic;
use Bio::SeqI;
use Bio::EnsEMBL::Utils::EMBL::GeneWrapper;
use Bio::EnsEMBL::Utils::EMBL::TranscriptWrapper;

@ISA = qw(Bio::EnsEMBL::Root Bio::SeqI);


=head2 new

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Example    : $sw = new Bio::EnsEMBL::Utils::EMBL::SliceWrapper($slice);
  Description: Creates a new SliceWrapper 'wrapped' around a slice object
  Returntype : Bio::EnsEMBL::Utils::EMBL::SliceWrapper
  Exceptions : none
  Caller     : embl dumping scripts 

=cut

sub new {
  my ($class, $slice) = @_;

  my $self = $class->SUPER::new();

  $self->{_skip_hash} = {};


  $self->slice($slice);
  
  $self->annotation(new Bio::Annotation::Collection);

  return $self;
}


=head2 slice

  Arg [1]    : (optional) Bio::EnsEMBL::Slice
  Example    : $slice = $slice_wrapper->slice();
  Description: getter / setter for the slice this is wrapped around
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : thrown if slice arg is not a Bio::PrimarySeqI
  Caller     : internal

=cut

sub slice {
  my ($self, $slice) = @_;

  if($slice) {
    unless($slice->isa("Bio::PrimarySeqI")) {
      $self->throw("[$slice] is not a Bio::PrimarySeqI");
    }
    
    $self->{_slice} = $slice;
  }

  return $self->{_slice};
}


=head2 primary_seq

  Arg [1]    : none
  Example    : none
  Description: Wrapper method - returns the slice object.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : none
  Caller     : Bio::EnsEMBL::Utils::EMBL::EMBL_Dump

=cut

sub primary_seq {
  my ($self) = @_;

  return $self->slice();
}


=head2 annotation

  Arg [1]    : (optional) Bio::Annotation::Collection $annotation
  Example    : none
  Description: Wrapper method - gets/sets annotation object contained 
               in wrapper
  Returntype : Bio::EnsEMBL::Annotation
  Exceptions : none
  Caller     : Bio::EnsEMBL::Utils::EMBL::EMBL_Dump

=cut

sub annotation {
  my ($self, $annotation) = @_;

  if($annotation) {
    unless($annotation->isa("Bio::Annotation::Collection")) {
      $self->throw("[$annotation] is not a Bio::Annotation::Collection");
    }

    $self->{_annotation} = $annotation;
  }

  return $self->{_annotation};
}


=head2 htg_phase

  Arg [1]    : none
  Example    : none
  Description: wrapper function - returns htg_phase
  Returntype : int
  Exceptions : none
  Caller     : Bio::SeqIO::embl

=cut

sub htg_phase {
  my $self = shift;

  return 4;
}


=head2 sv

  Arg [1]    : none
  Example    : none
  Description: Wrapper function - returns sequence version
  Returntype : string
  Exceptions : none
  Caller     : Bio::SeqIO::embl

=cut

sub sv {
  my $self = shift;

  return '';
}


=head2 ac

  Arg [1]    : none
  Example    : $accession = $slice_wrapper->ac()
  Description: Wrapper function - returns a generated accession for this 
               slice, of the form "$chr_name : $chr_start - $chr_end"
  Returntype : string
  Exceptions : none
  Caller     : Bio::SeqIO::embl

=cut

sub ac {
  my $self = shift;

  my $chr   = $self->slice->chr_name();
  my $start = $self->slice->chr_start();
  my $end   = $self->slice->chr_end();
 
  return "Chromosome $chr $start to $end";
}


=head2 id

  Arg [1]    : none
  Example    : none
  Description: wrapper function - returns value of ac function
  Returntype : string
  Exceptions : none
  Caller     : Bio::SeqIO::embl

=cut

sub id {
  my $self = shift;

  return $self->slice->name;
}


=head2 each_date

  Arg [1]    : none
  Example    : none
  Description: wrapper function - returns current time
  Returntype : none
  Exceptions : none
  Caller     : Bio::SeqIO::embl

=cut

sub each_date {
  my $self = shift;

  return (time);
}


=head2 division

  Arg [1]    : none
  Example    : none
  Description: wrapper function - returns undef
  Returntype : none
  Exceptions : none
  Caller     : Bio::SeqIO::embl

=cut

sub division {
  my $self = shift;
#  ???
  return undef;
}


=head2 moltype

  Arg [1]    : none
  Example    : none
  Description: Wrapper function - returns molecule type of contained slice
  Returntype : string
  Exceptions : none
  Caller     : Bio::SeqIO::embl

=cut

sub moltype {
  my $self = shift;

  return $self->slice()->moltype();
}



=head2 species

  Arg [1]    : none
  Example    : none
  Description: Wrapper function - returns the species of the database the
               slice was created from.
  Returntype : Bio::Species
  Exceptions : none
  Caller     : Bio::SeqIO::embl

=cut

sub species {
  my $self = shift;

  return $self->slice->adaptor->db->get_MetaContainer->get_Species();
}


=head2 skip_SeqFeature

  Arg [1]    : string $name 
               the type of seqfeature to be skipped
  Arg [2]    : (optional) boolean $bool
               True value if the feature type defined by the $name arg should
               be skipped false otherwise.
  Example    : $slice_wrapper->skip_SeqFeature('external', 1);
  Description: Turns on and off which features should be included in the embl
               dump of this slice wrapper.
  Returntype : the boolean value assigned to features of type $name
  Exceptions : none
  Caller     : Embl dumping scripts, top_SeqFeatures

=cut

sub skip_SeqFeature {
  my ($self, $name, $bool) = @_;
    
  if(defined $bool) {
    $self->{_skip_hash}->{$name} = $bool;
  }

  return 0 unless defined $self->{_skip_hash}->{$name};

  return $self->{_skip_hash}->{$name};
}  


=head2 top_SeqFeatures

  Arg [1]    : none
  Example    : none
  Description: Wrapper function - 
               returns SeqFeatures on the slice which have not been switched 
               off via the skip_SeqFeature method.  Genes may be included in 
               the list of SeqFeatures, but they will be wrapped in 
               GeneWrappers so they may be properly embl dumped.
  Returntype : list of Bio::SeqFeature objects
  Exceptions : none
  Caller     : Bio::SeqIO::

=cut

sub top_SeqFeatures {
  my $self = shift;

  my @sfs;

  my $slice_length = $self->slice->length();

  unless( $self->skip_SeqFeature('meta') ) {
    my $sf = new Bio::SeqFeature::Generic();
    $sf->start(1);
    $sf->end($slice_length);
    $sf->strand(1);
    $sf->primary_tag( "source" );
    my $species = $self->species;
    $sf->add_tag_value('organism', $species->common_name());
    $sf->add_tag_value('classification', 
		       join(', ', $species->classification()));
    push @sfs, $sf;
  }

  unless($self->skip_SeqFeature('similarity')) {
    push @sfs, @{$self->slice->get_all_SimilarityFeatures()};
  }
  unless($self->skip_SeqFeature('repeat')) {
    push @sfs, @{$self->slice->get_all_RepeatFeatures()};
  }
  unless($self->skip_SeqFeature('external')) {
    push @sfs, @{$self->slice->get_all_ExternalFeatures()};
  }

  unless($self->skip_SeqFeature('snp')) {
    push @sfs, @{$self->slice->get_all_SNPs};
  }

  #filter out features overlapping slice boundary
  my @out = grep { $_->start >0 && $_->end <= $slice_length } @sfs;

  #transcripts and genes are allowed to overlap boundary
  unless($self->skip_SeqFeature('prediction')) {
    push @out, map { new Bio::EnsEMBL::Utils::EMBL::TranscriptWrapper($_) }
        @{$self->slice->get_all_PredictionTranscripts};
  }
  unless($self->skip_SeqFeature('gene')) {
    push @out, map { new Bio::EnsEMBL::Utils::EMBL::GeneWrapper($_) } 
        @{$self->slice->get_all_Genes()};
  }

  return @out;
}


=head2 seq

  Arg [1]    : none
  Example    : none
  Description: Wrapper function returns the result of calling the seq method
               on the contained slice
  Returntype : string
  Exceptions : none
  Caller     : Bio::SeqIO::embl

=cut

sub seq {
  my ($self, @args) = shift;

  return $self->slice()->seq(@args);
}


=head2 length

  Arg [1]    : none
  Example    : none
  Description: Wrapper function - returns the result of calling length on 
               the contained slice
  Returntype : int
  Exceptions : none
  Caller     : Bio::SeqIO::embl

=cut

sub length {
  my $self = shift;

  return $self->slice()->length();
}


=head2 desc

  Arg [1]    : none
  Example    : none
  Description: Wrapper function - returns a suitable embl description for the
               contained slice.
  Returntype : string
  Exceptions : none
  Caller     : Bio::SeqIO::embl 

=cut

sub desc {
  my $self = shift;

  return 'Reannotated sequence via Ensembl';
}

1;   
