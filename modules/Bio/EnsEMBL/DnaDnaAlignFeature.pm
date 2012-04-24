=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::DnaDnaAlignFeature - Ensembl specific dna-dna pairwise
alignment feature

=head1 SYNOPSIS

  See BaseAlignFeature

=cut


package Bio::EnsEMBL::DnaDnaAlignFeature;

use strict;

use Bio::EnsEMBL::BaseAlignFeature;

use vars qw(@ISA);
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw( Bio::EnsEMBL::BaseAlignFeature );


=head2 new

  Arg [..]   : List of named arguments. (-pair_dna_align_feature_id) defined
               in this constructor, others defined in BaseFeaturePair and 
               SeqFeature superclasses.  
  Example    : $daf = new DnaDnaAlignFeature(-cigar_string => '3M3I12M');
  Description: Creates a new DnaDnaAlignFeature using either a cigarstring or
               a list of ungapped features.  
  Returntype : Bio::EnsEMBL::DnaDnaAlignFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my ($pair_dna_align_feature_id) = rearrange([qw(PAIR_DNA_ALIGN_FEATURE_ID)], @_);
  if (defined $pair_dna_align_feature_id){
      $self->{'pair_dna_align_feature_id'} = $pair_dna_align_feature_id;
  }
  return $self;
}


=head2 pair_dna_align_feature_id

  Arg[1]     : (optional) String $arg - value to set
  Example    : $self->pair_dna_align_feature_id($pair_feature_id);
  Description: Getter/setter for attribute 'pair_dna_align_feature_id'
               The id of the dna feature aligned
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub pair_dna_align_feature_id{
    my ($self, $arg) = @_;
    if (defined $arg){
	$self->{pair_dna_align_feature_id} = $arg;
    }
    return $self->{pair_dna_align_feature_id};
}

=head2 _hit_unit

  Arg [1]    : none
  Description: PRIVATE implementation of abstract superclass method.  Returns
               1 as the 'unit' used for the hit sequence. 
  Returntype : int
  Exceptions : none
  Caller     : Bio::EnsEMBL::BaseAlignFeature
  Status     : Stable

=cut

sub _hit_unit {
  return 1;
}



=head2 _query_unit

  Arg [1]    : none
  Description: PRIVATE implementation of abstract superclass method Returns
               1 as the 'unit' used for the hit sequence.
  Returntype : int
  Exceptions : none
  Caller     : Bio::EnsEMBL::BaseAlignFeature
  Status     : Stable

=cut

sub _query_unit {
  return 1;
}

=head2 restrict_between_positions

  Arg [1]    : int $start
  Arg [2]    : int $end
  Arg [3]    : string $flags
               SEQ = $start and $end apply to the seq sequence
                     i.e. start and end methods
               HSEQ = $start and $end apply to the hseq sequence
                      i.e. hstart and hend methods
  Example    : $daf->restrict_between_positions(150,543,"SEQ")
  Description: Build a new DnaDnaAlignFeature object that fits within
               the new specified coordinates and sequence reference, cutting
               any pieces hanging upstream and downstream.
  Returntype : Bio::EnsEMBL::DnaDnaAlignFeature object
  Exceptions : 
  Caller     : 
  Status     : Stable

=cut

sub restrict_between_positions {
  my ($self,$start,$end,$seqref) = @_;

  unless (defined $start && $start =~ /^\d+$/) {
    $self->throw("The first argument is not defined or is not an integer");
  }
  unless (defined $end && $end =~ /^\d+$/) {
    $self->throw("The second argument is not defined or is not an integer");
  }
  unless (defined $seqref &&
          ($seqref eq "SEQ" || $seqref eq "HSEQ")) {
    $self->throw("The third argument is not defined or is not equal to 'SEQ' or 'HSEQ'");
  }

# symbolic method references should be forbidden!
# need to be rewrite at some stage.

  my ($start_method1,$end_method1,$strand_method1,$start_method2,$end_method2,$strand_method2) =
    qw(start end strand hstart hend hstrand);

  if ($seqref eq "HSEQ") {
    ($start_method1,$end_method1,$strand_method1,$start_method2,$end_method2,$strand_method2) =
    qw(hstart hend hstrand start end strand);
  }

  my @restricted_features;

  foreach my $ungapped_feature ($self->ungapped_features) {

    if ($ungapped_feature->$start_method1() > $end ||
        $ungapped_feature->$end_method1() < $start) {

      next;

    } elsif ($ungapped_feature->$end_method1() <= $end &&
             $ungapped_feature->$start_method1() >= $start) {

      push @restricted_features, $ungapped_feature;

    } else {

      if ($ungapped_feature->$strand_method1() eq $ungapped_feature->$strand_method2()) {

        if ($ungapped_feature->$start_method1() < $start) {

          my $offset = $start - $ungapped_feature->$start_method1();
          $ungapped_feature->$start_method1($start);
          $ungapped_feature->$start_method2($ungapped_feature->$start_method2() + $offset);

        }
        if ($ungapped_feature->$end_method1() > $end) {

          my $offset = $ungapped_feature->$end_method1() - $end;
          $ungapped_feature->$end_method1($end);
          $ungapped_feature->$end_method2($ungapped_feature->$end_method2() - $offset);

        }
      } else {

        if ($ungapped_feature->$start_method1() < $start) {

          my $offset = $start - $ungapped_feature->$start_method1();
          $ungapped_feature->$start_method1($start);
          $ungapped_feature->$end_method2($ungapped_feature->$end_method2() - $offset);

        }
        if ($ungapped_feature->$end_method1() > $end) {

          my $offset = $ungapped_feature->$end_method1() - $end;
          $ungapped_feature->$end_method1($end);
          $ungapped_feature->$start_method2($ungapped_feature->$start_method2() + $offset);

        }
      }
      
      push @restricted_features, $ungapped_feature;
    }
  }

  if (scalar @restricted_features) {
    my $DnaDnaAlignFeature = new Bio::EnsEMBL::DnaDnaAlignFeature('-features' =>\@restricted_features);
    if (defined $self->slice) {
      $DnaDnaAlignFeature->slice($self->slice);
    }
    if (defined $self->hslice) {
      $DnaDnaAlignFeature->hslice($self->hslice);
    }
    return $DnaDnaAlignFeature;
  } else {
    return undef;
  }
}

=head2 alignment_strings

  Arg [1]    : list of string $flags
               FIX_SEQ = does not introduce gaps (dashes) in seq aligned sequence
                         and delete the corresponding insertions in hseq aligned sequence
               FIX_HSEQ = does not introduce gaps (dashes) in hseq aligned sequence
                         and delete the corresponding insertions in seq aligned sequence
               NO_SEQ = return the seq aligned sequence as an empty string
               NO_HSEQ = return the hseq aligned sequence as an empty string
               This 2 last flags would save a bit of time as doing so no querying to the core
               database in done to get the sequence.
  Example    : $daf->alignment_strings or
               $daf->alignment_strings("FIX_HSEQ") or
               $daf->alignment_strings("NO_SEQ","FIX_SEQ")
  Description: Allows to rebuild the alignment string of both the seq and hseq sequence
               using the cigar_string information and the slice and hslice objects
  Returntype : array reference containing 2 strings
               the first corresponds to seq
               the second corresponds to hseq
  Exceptions : 
  Caller     : 
  Status     : Stable

=cut


sub alignment_strings {
  my ( $self, @flags ) = @_;

  # set the flags
  my $seq_flag = 1;
  my $hseq_flag = 1;
  my $fix_seq_flag = 0;
  my $fix_hseq_flag = 0;

  for my $flag ( @flags ) {
    $seq_flag = 0 if ($flag eq "NO_SEQ");
    $hseq_flag = 0 if ($flag eq "NO_HSEQ");
    $fix_seq_flag = 1 if ($flag eq "FIX_SEQ");
    $fix_hseq_flag = 1 if ($flag eq "FIX_HSEQ");
  } 

  my ($seq, $hseq);
  $seq = $self->slice->subseq($self->start, $self->end, $self->strand) if ($seq_flag || $fix_seq_flag);  
  $hseq = $self->hslice->subseq($self->hstart, $self->hend, $self->hstrand) if ($hseq_flag || $fix_hseq_flag);
  
  my $rseq= "";
  # rseq - result sequence
  my $rhseq= "";
  # rhseq - result hsequence

  my $seq_pos = 0;
  my $hseq_pos = 0;

  my @cig = ( $self->cigar_string =~ /(\d*[DIM])/g );

  for my $cigElem ( @cig ) {
    my $cigType = substr( $cigElem, -1, 1 );
    my $cigCount = substr( $cigElem, 0 ,-1 );
    $cigCount = 1 unless $cigCount;

    if( $cigType eq "M" ) {
        $rseq .= substr( $seq, $seq_pos, $cigCount ) if ($seq_flag);
        $rhseq .= substr( $hseq, $hseq_pos, $cigCount ) if ($hseq_flag);
      $seq_pos += $cigCount;
      $hseq_pos += $cigCount;
    } elsif( $cigType eq "D" ) {
      if( ! $fix_seq_flag ) {
        $rseq .=  "-" x $cigCount if ($seq_flag);
        $rhseq .= substr( $hseq, $hseq_pos, $cigCount ) if ($hseq_flag);
      }
      $hseq_pos += $cigCount;
    } elsif( $cigType eq "I" ) {
      if( ! $fix_hseq_flag ) {
        $rseq .= substr( $seq, $seq_pos, $cigCount ) if ($seq_flag);
        $rhseq .= "-" x $cigCount if ($hseq_flag);
      }
      $seq_pos += $cigCount;
    }
  }
  return [ $rseq,$rhseq ];
}

=head2 get_SimpleAlign

  Arg [1]    : list of string $flags
               translated = by default, the sequence alignment will be on nucleotide. With translated flag
                            the aligned sequences are translated.
               uc = by default aligned sequences are given in lower cases. With uc flag, the aligned 
                    sequences are given in upper cases.
  Example    : $daf->get_SimpleAlign or
               $daf->get_SimpleAlign("translated") or
               $daf->get_SimpleAlign("translated","uc")
  Description: Allows to rebuild the alignment string of both the seq and hseq sequence
               using the cigar_string information and the slice and hslice objects
  Returntype : a Bio::SimpleAlign object
  Exceptions : 
  Caller     : 
  Status     : Stable

=cut

sub get_SimpleAlign {
  my ( $self, @flags ) = @_;

  # setting the flags
  my $uc = 0;
  my $translated = 0;

  for my $flag ( @flags ) {
    $uc = 1 if ($flag =~ /^uc$/i);
    $translated = 1 if ($flag =~ /^translated$/i);
  }

  my $sa = Bio::SimpleAlign->new();

  #Hack to try to work with both bioperl 0.7 and 1.2:
  #Check to see if the method is called 'addSeq' or 'add_seq'
  my $bio07 = 0;
  if(!$sa->can('add_seq')) {
    $bio07 = 1;
  }

  my ($sb_seq,$qy_seq) = @{$self->alignment_strings};

  my $loc_sb_seq = Bio::LocatableSeq->new(-SEQ    => $uc ? uc $sb_seq : lc $sb_seq,
                                          -START  => $self->seq_region_start,
                                          -END    => $self->seq_region_end,
                                          -ID     => $self->seqname,
                                          -STRAND => $self->strand);

  $loc_sb_seq->seq($uc ? uc $loc_sb_seq->translate->seq
                   : lc $loc_sb_seq->translate->seq) if ($translated);

  my $loc_qy_seq = Bio::LocatableSeq->new(-SEQ    => $uc ? uc $qy_seq : lc $qy_seq,
                                          -START  => $self->hseq_region_start,
                                          -END    => $self->hseq_region_end,
                                          -ID     => $self->hseqname,
                                          -STRAND => $self->hstrand);

  $loc_qy_seq->seq($uc ? uc  $loc_qy_seq->translate->seq
                   : lc $loc_qy_seq->translate->seq) if ($translated);

  if($bio07) {
    $sa->addSeq($loc_sb_seq);
    $sa->addSeq($loc_qy_seq);
  } else {
    $sa->add_seq($loc_sb_seq);
    $sa->add_seq($loc_qy_seq);
  }

  return $sa;
}

1;
