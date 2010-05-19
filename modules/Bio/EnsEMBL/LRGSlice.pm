=head1 LICENSE

  Copyright (c) 1999-2010 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::LRGSlice - Arbitary Slice of a genome

=head1 SYNOPSIS

  $sa = $db->get_SliceAdaptor;

  $slice =
    $sa->fetch_by_region( 'LRG', 'LRG3');

  # get some attributes of the slice
  my $seqname = $slice->seq_region_name();
  my $start   = $slice->start();
  my $end     = $slice->end();

  # get the sequence from the slice
  my $seq = $slice->seq();

  # get some features from the slice
  foreach $gene ( @{ $slice->get_all_Genes } ) {
    # do something with a gene
  }

  foreach my $feature ( @{ $slice->get_all_DnaAlignFeatures } ) {
    # do something with dna-dna alignments
  }

=head1 DESCRIPTION

A LRG Slice object represents a region of a genome.  It can be used to retrieve
sequence or features from an area of interest.

=head1 METHODS

=cut

package Bio::EnsEMBL::LRGSlice;
use vars qw(@ISA);
use strict;

use Bio::PrimarySeqI;

my $reg = "Bio::EnsEMBL::Registry";

use Bio::EnsEMBL::Slice;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Slice);

sub new{
  my $class = shift;

  my $self = bless {}, $class ;

  my $slice = $self = $class->SUPER::new( @_);

 return $self;
}

sub stable_id {
    my $self = shift;
    return $self->seq_region_name;
}


sub display_xref {
    my $self = shift;
    return $self->seq_region_name;
}

sub feature_Slice {
  my $self = shift;
  return $self->{_chrom_slice} if defined($self->{_chrom_slice});

  my $max=-99999999999;
  my $min=9999999999;
  my $chrom;
  my $strand;

#  print STDERR "working out feature slcie\n";
  foreach my $segment (@{$self->project('chromosome')}) {
    my $from_start = $segment->from_start();
    my $from_end    = $segment->from_end();
    my $to_name    = $segment->to_Slice->seq_region_name();
    $chrom = $to_name;

    my $to_start    = $segment->to_Slice->start();
    my $to_end    = $segment->to_Slice->end();
    if($to_start > $max){
      $max = $to_start;
    }
    if($to_start < $min){
      $min = $to_start;
    }
    if($to_end > $max){
      $max = $to_end;
    }
    if($to_end <  $min){
      $min = $to_end;
    }
    my $ori        = $segment->to_Slice->strand();
    $strand = $ori;
  }
  if(!defined($chrom)){
    warn "Could not project to chromosome for ".$self->name."??\n";
    return undef;
  }
  my $chrom_slice = $self->adaptor->fetch_by_region("chromosome",$chrom, $min, $max, $strand);
  $self->{_chrom_slice} = $chrom_slice;
  return $self->{_chrom_slice};
}

sub DESTROY{
}


1;
