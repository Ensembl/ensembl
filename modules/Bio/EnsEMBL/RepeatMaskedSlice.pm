#
# Ensembl module for Bio::EnsEMBL::RepeatMaskedSlice
#
#
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::RepeatMaskedSlice - Arbitary Slice of a genome

=head1 SYNOPSIS

   $sa = $db->get_SliceAdaptor;

   $slice = $sa->fetch_by_region('chromosome', 'X', 1_000_000, 2_000_000);

   $repeat_masked_slice = $sa->get_repeatmasked_seq();

   #get repeat masked sequence:
   my $dna = $repeat_masked_slice->seq();
   $dna    = $repeat_masked_slice->subseq(1, 1000);


=head1 DESCRIPTION

This is a specialised B<Bio::EnsEMBL::Slice> class that is used to 
retrieve repeat masked genomic sequence rather than normal genomic 
sequence.

=head1 AUTHOR - Graham McVicker

=head1 CONTACT

This modules is part of the Ensembl project http://www.ensembl.org

Questions can be posted to the ensembl-dev mailing list:
ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut


package Bio::EnsEMBL::RepeatMaskedSlice;

use strict;
use warnings;

use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);


use vars qw(@ISA);

@ISA = ('Bio::EnsEMBL::Slice');



=head2 new

  Arg [-REPEAT_MASK] : The logic name of the repeats to be used for masking.
                      If not provided, all repeats in the database are used.
  Arg [...]  : Named superclass arguments. See B<Bio::EnsEMBL::Slice>.
  Example    : my $slice = Bio::EnsEMBL::RepeatMaskedSlice->new
                  (-START  => $start,
                   -END    => $end,
                   -STRAND => $strand,
                   -SEQ_REGION_NAME => $seq_region,
                   -COORD_SYSTEM  => $cs,
                   -ADAPTOR => $adaptor,
                   -REPEAT_MASK => ['repeat_masker'],
                   -SOFT_MASK => 1);
  Description: Creates a Slice which behaves exactly as a normal slice but
               that returns repeat masked sequence from the seq method.
  Returntype : Bio::EnsEMBL::RepeatMaskedSlice
  Exceptions : none
  Caller     : RawComputes (PredictionTranscript creation code).

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($logic_names, $soft_mask) = rearrange(['REPEAT_MASK',
                                            'SOFT_MASK'], @_);

  my $self = $class->SUPER::new(@_);


  $logic_names ||= [''];
  if(ref($logic_names) ne 'ARRAY') {
    throw("Reference to list of logic names argument expected.");
  }

  $self->{'repeat_mask_logic_names'} = $logic_names;
  $self->{'soft_mask'} = $soft_mask;

  return $self;
}


=head2 repeat_mask_logic_names

  Arg [1]    : reference to list of strings $logic_names (optional)
  Example    : $rm_slice->repeat_mask_logic_name(['repeat_masker']);
  Description: Getter/Setter for the logic_names of the repeats that are used
               to mask this slices sequence.
  Returntype : reference to list of strings
  Exceptions : none
  Caller     : seq() method

=cut

sub repeat_mask_logic_names {
  my $self = shift;

  if(@_) {
    my $array = shift;
    if(ref($array) ne 'ARRAY') {
      throw('Reference to list of logic names argument expected.');
    }
  }
  
  return $self->{'repeat_mask_logic_names'};
}


=head2 soft_mask

  Arg [1]    : boolean $soft_mask (optional)
  Example    : $rm_slice->soft_mask(0);
  Description: Getter/Setter which is used to turn on/off softmasking of the
               sequence returned by seq.
  Returntype : boolean
  Exceptions : none
  Caller     : seq() method

=cut

sub soft_mask {
  my $self = shift;
  $self->{'soft_mask'} = shift if(@_);
  return $self->{'soft_mask'} || 0;
}



=head2 seq

  Arg [1]    : none
  Example    : print $rmslice->seq(), "\n";
  Description: Retrieves the entire repeat masked sequence for this slice.
               See also the B<Bio::EnsEMBL::Slice> implementation of this 
               method.
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub seq {
  my $self = shift;
  #
  # get all the features
  #
  my $logic_names = $self->repeat_mask_logic_names();
  my $soft_mask   = $self->soft_mask();

  my $repeats = [];

  foreach my $l (@$logic_names) {
    push @{$repeats}, @{$self->get_all_RepeatFeatures($l)};
  }

  #
  # get the dna
  #
  my $dna = $self->SUPER::seq(@_);

  #
  # mask the dna
  #
  $self->_mask_features(\$dna,$repeats,$soft_mask);
  return $dna;
}



=head2 subseq

  Arg [1]    : none
  Example    : print $rmslice->subseq(1, 1000);
  Description: Retrieves a repeat masked sequence from a specified subregion
               of this slice.  See also the B<Bio::EnsEMBL::Slice> 
               implementation of this method.
  Returntype : string
  Exceptions : none
  Caller     : general

=cut


sub subseq {
  my $self = shift;
  #
  # get all the features
  #
  my $logic_names = $self->repeat_mask_logic_names();
  my $soft_mask   = $self->soft_mask();

  my $repeats = [];

  foreach my $l (@$logic_names) {
    push @{$repeats}, @{$self->get_all_RepeatFeatures($l)};
  }

  #
  # get the dna
  #
  my $dna = $self->SUPER::subseq(@_);

  #
  # mask the dna
  #
  $self->_mask_features(\$dna,$repeats,$soft_mask);
  return $dna;
}


1;
