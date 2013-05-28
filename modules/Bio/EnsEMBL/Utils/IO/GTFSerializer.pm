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

my %strand_conversion = ( '1' => '+', '0' => '.', '-1' => '-');

=head2 print_feature

    Arg [1]    : Bio::EnsEMBL::Transcript
    Example    : $serializer->print_feature($transcript)
    Description: 
    Returntype : none

=cut

sub print_feature {
  my $self = shift;
  my $transcript = shift;

  throw(sprintf "Feature is of type %s. Cannot write non transcripts to GTF", ref($transcript))
    unless check_ref($transcript, "Bio::EnsEMBL::Transcript");

  # filehandle is inherited
  my $fh = $self->{'filehandle'};

  my $slice = $transcript->slice();
  my $sliceoffset = $slice->start-1;
  my $idstr = $slice->seq_region_name; 

  my @startcs =  $self->_make_start_codon_features($transcript);
  my @endcs   =  $self->_make_stop_codon_features($transcript);
  my ($hasstart, $hasend) = $self->_check_start_and_stop($transcript);
  #
  # production does not use this option
  #
  # if (!$include_codons) {
  #   $hasstart = $hasend = 0;
  # }

  # TODO: ask Andy if this is safe
  my $dbname = $transcript->adaptor()->dbc()->dbname();
  my $vegadb = $dbname =~ /vega/;
  my $gene = $transcript->get_Gene();

  my ($biotype_display, $transcript_biotype);
  {
    no warnings 'uninitialized';
    $biotype_display = $vegadb ? $gene->status . '_' . $gene->biotype : $gene->biotype;
    $transcript_biotype = $vegadb ? $transcript->status . '_' . $transcript->biotype : $transcript->biotype;
  }

  my @translateable_exons = 
    @{$transcript->get_all_translateable_Exons} 
      if $transcript->translation;

  my ($count, $intrans, $instop) = (1, 0, 0);

  foreach my $exon (@{$transcript->get_all_Exons}) {
    my $strand = $strand_conversion{$exon->strand};

    print $fh 
      # Column 1 - seqname, the name of the sequence/chromosome the feature is on. Landmark for start below
      $idstr . "\t" . 
      # Column 2 - source
      $transcript_biotype . "\t" . 
      # Column 3 - feature type name
      'exon' . "\t" . 
      # Column 4 - start, the start coordinate of the feature
      ($exon->start+$sliceoffset) . "\t". 
      # Column 5 - end, coordinates (absolute) for the end of this feature
      ($exon->end+$sliceoffset) . "\t". 
      # Column 6 - score, for variations only
      "." . "\t". 
      # Column 7 - strand, forward (+) or reverse (-)
      $strand . "\t". 
      # Column 8 - frame (reading phase), what base of this feature is the first base of a codon
      "." . "\t";
    # Column 9 - attribute, a ;-separated list of key-value pairs (additional feature info)
    $self->_print_attribs($gene, $biotype_display, $transcript, $count, 'exon', $exon, $vegadb);
    print $fh "\n";

    $intrans = 1 
      if $transcript->translation and 
	$exon == $transcript->translation->start_Exon;

    if ($intrans) {
      # print the CDS of this exon

      my $cdsexon = shift @translateable_exons;
      #
      # Here is computing the value of the GTF frame (taking into
      # account the Ensembl convention), but it's misleadingly called phase
      #
      my $phase = $cdsexon->phase;
      if ($cdsexon->phase == 1) {
        $phase = 2;
      } elsif ($cdsexon->phase == 2) {
        $phase = 1;
      } elsif ($cdsexon->phase == -1) {
        $phase = 0;
      }

      my $exon_start = $cdsexon->start;
      my $exon_end   = $cdsexon->end;
      if ($transcript->translation && 
          $hasend && 
          ($exon->end >= $endcs[0]->start && $exon->start <= $endcs[0]->end)) {

        if ($cdsexon->strand == 1) {
          $exon_end = $cdsexon->end - $endcs[0]->length;
        } else {
          $exon_start = $cdsexon->start + $endcs[0]->length;
        }
      }

      if ($exon_start <= $cdsexon->end &&
          $exon_end >= $cdsexon->start &&
          !$instop) {
        print $fh $idstr . "\t" . 
	  $transcript_biotype . "\t" . 
          'CDS' . "\t" . 
          ($exon_start + $sliceoffset) . "\t". 
          ($exon_end + $sliceoffset) . "\t". 
          "." . "\t". 
          $strand . "\t". 
          $phase . "\t";
        $self->_print_attribs($gene, $biotype_display, $transcript, $count, 'CDS');
        print $fh "\n";
      }
    }
    if ($transcript->translation && 
        $exon == $transcript->translation->start_Exon && $hasstart) {
      my $tmpcnt = $count;
      foreach my $startc (@startcs) {

        print $fh $idstr . "\t" . 
	  $transcript_biotype . "\t" . 
          'start_codon' . "\t" . 
          ($startc->start+$sliceoffset) . "\t". 
          ($startc->end+$sliceoffset) . "\t". 
          "." . "\t". 
          $strand . "\t". 
          $startc->phase . "\t";

        $self->_print_attribs($gene, $biotype_display, $transcript, $tmpcnt++, 'start_codon');
        print $fh "\n";
      }
    }
    if ($transcript->translation && 
        ($exon == $transcript->translation->end_Exon)) {
      if ($hasend) {
        my $tmpcnt = $count - $#endcs;

        foreach my $endc (@endcs) {

          print $fh $idstr . "\t" . 
	    $transcript_biotype . "\t" . 
            'stop_codon' . "\t" . 
            ($endc->start+$sliceoffset) . "\t". 
            ($endc->end+$sliceoffset) . "\t". 
            "." . "\t". 
            $strand . "\t". 
            $endc->phase . "\t";

          $self->_print_attribs($gene, $biotype_display, $transcript, $tmpcnt++, 'stop_codon');
          print $fh "\n";
        }
      }
      $intrans = 0;
    }

    if (scalar(@endcs) && 
        ($exon->end >= $endcs[0]->start && $exon->start <= $endcs[0]->end)) {
      $instop = 1;
    }

    $count++;
  }

}

=head2 _print_attribs

    Arg []     : 
    Example    : 
    Description: 
    Returntype : None

=cut

sub _print_attribs {
  my ($self, $gene, $gene_biotype, $transcript, $count, $type, $exon, $vegadb) = @_;

  my $gene_name;
  $gene_name = $gene->external_name;
  $gene_name =~ s/^[A-Z]{1,3}:// if $vegadb;

  my $trans_name;
  $trans_name = $transcript->external_name;
  $trans_name =~ s/^[A-Z]{1,3}:// if $vegadb;

  my $fh = $self->{'filehandle'};

  print $fh "\tgene_id \"" .  get_id_from_obj($gene) . "\";" .
            " transcript_id \"" . get_id_from_obj($transcript) . "\";";
  print $fh " exon_number \"$count\";";
  print $fh " gene_name \"" . $gene_name . "\";" if ($gene_name);
  print $fh " gene_biotype \"" . $gene_biotype ."\";";
  print $fh " transcript_name \"" . $trans_name . "\";" if ($trans_name);
  if ($type eq 'CDS') {
    print $fh ' protein_id "' . get_id_from_obj($transcript->translation) . '";';
  }
  if($exon) {
    printf $fh ' exon_id "%s";', get_id_from_obj($exon);
  }
  return;
}

sub get_id_from_obj {
  my ($obj) = @_;
  my $id = $obj->stable_id();
  $id = $obj->dbID() unless defined $id;
  return $id;
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

  my ($attrib) = @{ $trans->slice()->get_all_Attributes('codon_table') };

  my $codon_table_id = $attrib->value()
    if defined $attrib;
  $codon_table_id ||= 1; # default vertebrate codon table

  my $codon_table = Bio::Tools::CodonTable->new( -id => $codon_table_id );

  $has_start = 0 unless $codon_table->is_start_codon($startseq);
  $has_end = 0 unless $codon_table->is_ter_codon($endseq);

  return ($has_start, $has_end);

}

1;
