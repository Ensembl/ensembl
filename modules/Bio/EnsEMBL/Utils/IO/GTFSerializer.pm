=pod

=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 NAME

Bio::EnsEMBL::Utils::IO::GTFSerializer - Transcript to GTF converter

=head1 SYNOPSIS

use Bio::EnsEMBL::Utils::IO::GTFSerializer;

my $serializer = Bio::EnsEMBL::Utils::IO::GTFSerializer->new($output_fh);

=head1 DESCRIPTION

Subclass of Serializer that can turn a transcript into a series of lines 
for the GTF format.

=cut

package Bio::EnsEMBL::Utils::IO::GTFSerializer;

use strict;
use warnings;

use  Bio::Tools::CodonTable;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::Utils::IO::FeatureSerializer;
use Bio::EnsEMBL::Utils::Scalar qw/check_ref/;

use base qw(Bio::EnsEMBL::Utils::IO::FeatureSerializer);

my %strand_conversion = ( '1' => '+', '0' => '?', '-1' => '-');

=head2 print_feature

    Arg [1]    : Bio::EnsEMBL::Transcript
    Example    : $serializer->print_feature($transcript)
    Description: 
    Returntype : none

=cut

sub print_feature {
  my $self = shift;
  my $feature = shift;

  throw(sprintf "Feature is of type %s. Cannot write non transcripts to GTF", ref($feature))
    unless check_ref($feature, "Bio::EnsEMBL::Transcript");

  my $text_buffer = "";

  #filehandle is inherited
  my $fh = $self->{'filehandle'};
  print $fh $text_buffer;

}

=head2 _make_start_codon_features

    Arg [1]    : Bio::EnsEMBL::Transcript
    Example    : 
    Description: 
    Returntype : Array

=cut

sub _make_start_codon_features {
  my ($self, $trans) = @_;
  defined $trans or
    throw("Transcript object not defined");

  return (()) unless $trans->translation;

  my @translateable = @{$trans->get_all_translateable_Exons};

  my @pepgencoords = $trans->pep2genomic(1,1);

  # cdna can come padded these days so allow gap at the start
  if($pepgencoords[0]->isa('Bio::EnsEMBL::Mapper::Gap')){
    shift @pepgencoords;
  }

  if(scalar(@pepgencoords) > 3) {
    throw(sprintf "Pep start for transcript %s does not map cleanly", $trans->display_id);
  }

  unless($pepgencoords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
    throw(sprintf "Pep start for transcript %s maps to gap", $trans->display_id);
  }
  unless($pepgencoords[$#pepgencoords]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
    throw(sprintf "Pep start (end of) for transcript %s maps to gap", $trans->display_id);
  }

  @translateable = @{$trans->get_all_translateable_Exons};
  my @startc_feat;
  my $phase = 0;
  foreach my $pepgencoord (@pepgencoords) {
    push @startc_feat, 
      new Bio::EnsEMBL::SeqFeature(-seqname => $trans->stable_id,
				   -source_tag => 'starttrans',
				   -primary_tag => 'similarity',
				   -start => $pepgencoord->start,
				   -end   => $pepgencoord->end,
				   -phase => $phase,
				   -strand => $translateable[0]->strand);
    $phase = 3 - ($pepgencoord->end - $pepgencoord->start + 1);
  }
  if ($translateable[0]->strand == 1) {
    @startc_feat = sort {$a->start <=> $b->start } @startc_feat;
  } else {
    @startc_feat = sort {$b->start <=> $a->start } @startc_feat;
  }

  return @startc_feat;

}

=head2 _make_stop_codon_features

    Arg [1]    : Bio::EnsEMBL::Transcript
    Example    : 
    Description: 
    Returntype : Array

=cut

sub _make_stop_codon_features {
  my ($self, $trans) = @_;

  defined $trans or
    throw("Transcript object not defined");

  return (()) unless $trans->translation;

  my @translateable = @{$trans->get_all_translateable_Exons};

  my $cdna_endpos = $trans->cdna_coding_end;

  my @pepgencoords = $trans->cdna2genomic($cdna_endpos-2,$cdna_endpos);

  if(scalar(@pepgencoords) > 3) {
    throw(sprintf "Pep end for transcript %s does not map cleanly", $trans->display_id);
  }
  unless($pepgencoords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
    throw(sprintf "Pep end for transcript %s maps to gap", $trans->display_id);
  }
  unless($pepgencoords[$#pepgencoords]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
    throw(sprintf "Pep end (end of) for transcript %s maps to gap", $trans->display_id);
  }

  my @stopc_feat;
  my $phase = 0;
  foreach my $pepgencoord (@pepgencoords) {
    push @stopc_feat, 
      new Bio::EnsEMBL::SeqFeature(-seqname => $trans->display_id,
				   -source_tag => 'endtrans',
				   -primary_tag => 'similarity',
				   -start => $pepgencoord->start,
				   -end   => $pepgencoord->end,
				   -phase => $phase,
				   -strand => $translateable[0]->strand);
    $phase = 3 - ($pepgencoord->end-$pepgencoord->start+1);
  }

  if ($translateable[0]->strand == 1) {
    @stopc_feat = sort {$a->start <=> $b->start } @stopc_feat;
  } else {
    @stopc_feat = sort {$b->start <=> $a->start } @stopc_feat;
  }

  return @stopc_feat;

}

=head2 _check_start_and_stop

    Arg [1]    : Bio::EnsEMBL::Transcript
    Example    : 
    Description: 
    Returntype : Array

=cut

sub _check_start_and_stop {
  my ($self, ,$trans) = @_;

  return (0,0) unless defined $trans->translation;

  my $tln = $trans->translation;

  my $coding_start = $trans->cdna_coding_start;
  my $coding_end   = $trans->cdna_coding_end;
  my $cdna_seq     = uc($trans->spliced_seq);

  my $startseq     = substr($cdna_seq,$coding_start-1,3);
  my $endseq       = substr($cdna_seq,$coding_end-3,3);

  my $has_start = 1;
  my $has_end = 1;

  # reimplemented because verterbrate specific
  # $has_start = 0  if ($startseq ne "ATG");
  # $has_end = 0 if ($endseq ne "TAG" && $endseq ne "TGA" && $endseq ne "TAA");

  my ($attrib) = @{ $self->slice()->get_all_Attributes('codon_table') };

  my $codon_table_id = $attrib->value()
    if defined $attrib;
  $codon_table_id ||= 1; # default vertebrate codon table

  my $codon_table = Bio::Tools::CodonTable->new( -id => $codon_table_id );

  $has_start = 0 unless $codon_table->is_start_codon($startseq);
  $has_end = 0 unless $codon_table->is_ter_codon($endseq);

  return ($has_start, $has_end);

}

1;
