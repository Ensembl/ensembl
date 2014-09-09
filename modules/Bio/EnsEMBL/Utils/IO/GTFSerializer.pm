
=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

=pod


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

use Bio::Tools::CodonTable;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::Utils::IO::FeatureSerializer;
use Bio::EnsEMBL::Utils::Scalar qw/check_ref/;
use Bio::EnsEMBL::Feature;

use base qw(Bio::EnsEMBL::Utils::IO::FeatureSerializer);

my %strand_conversion = ( '1' => '+', '0' => '.', '-1' => '-' );

sub print_main_header {
  my $self = shift;
  my $dba = shift;
  my $fh = $self->{'filehandle'};
  
  my $mc = $dba->get_MetaContainer();
  my $gc = $dba->get_GenomeContainer();

  # Get the build. name gives us GRCh37.p1 where as default gives us GRCh37
  my $assembly_name = $gc->get_assembly_name();
  print $fh "#!genome-build ${assembly_name}\n" if $assembly_name;

  # Get the build default
  my $version = $gc->get_version();
  print $fh "#!genome-version ${version}\n" if $version;

  # Get the date of the genome build
  my $assembly_date = $gc->get_assembly_date();
  print $fh "#!genome-date ${assembly_date}\n" if $assembly_date;

  # Get accession and only print if it is there
  my $accession = $gc->get_accession();
  if($accession) {
    my $accession_source = $mc->single_value_by_key('assembly.web_accession_source');
    print $fh "#!genome-build-accession ${accession_source}:${accession}\n";
  }

  # Genebuild last updated
  my $genebuild_last_date = $gc->get_genebuild_last_geneset_update();
  print $fh "#!genebuild-last-updated ${genebuild_last_date}\n" if $genebuild_last_date;

  return;
}

sub print_Gene {
  my ($self, $gene) = @_;

  #Print Gene summary
  my $slice           = $gene->slice();
  my $idstr           = $slice->seq_region_name;
  my $sliceoffset     = $slice->start - 1;
  my $vegadb          = 0;
  my $fh              = $self->{'filehandle'};
  my $biotype_display = q{};
  if($gene->adaptor()) {
    my $dbname = $gene->adaptor()->dbc()->dbname();
    $vegadb = $dbname =~ /vega/;
  }

  {
    no warnings 'uninitialized';
    $biotype_display = $vegadb ? $gene->status . '_' . $gene->biotype : $gene->biotype;
  }

  print $fh sprintf(qq{%s\t%s\tgene\t%d\t%d\t.\t%s\t.\t}, 
        $idstr, $gene->source, 
        ($gene->start()+$sliceoffset), ($gene->end()+$sliceoffset),
        ($strand_conversion{$gene->strand}));
  $self->_print_attribs($gene, $biotype_display, $gene, $biotype_display, 0, 'gene');
  print $fh "\n";

  # Now print all transcripts
  foreach my $t (@{$gene->get_all_Transcripts()}) {
    $self->print_feature($t, $gene);
  }
  return;
}

=head2 print_feature

    Arg [1]    : Bio::EnsEMBL::Transcript
    Example    : $serializer->print_feature($transcript)
    Description: 
    Returntype : none

=cut

sub print_feature {
  my $self       = shift;
  my $transcript = shift;
  my $gene = shift;

  throw( sprintf
           "Feature is of type %s. Cannot write non transcripts to GTF",
         ref($transcript)
  ) unless check_ref( $transcript, "Bio::EnsEMBL::Transcript" );

  # filehandle is inherited
  my $fh = $self->{'filehandle'};

  my $slice       = $transcript->slice();
  my $sliceoffset = $slice->start - 1;
  my $idstr       = $slice->seq_region_name;

  my @startcs = $self->_make_start_codon_features($transcript);
  my @endcs   = $self->_make_stop_codon_features($transcript);
  my ( $hasstart, $hasend ) = $self->_check_start_and_stop($transcript);

  my $vegadb = 0;
  if($transcript->adaptor()) {
    my $dbname = $transcript->adaptor()->dbc()->dbname();
    $vegadb = $dbname =~ /vega/;
  }
  $gene   ||= $transcript->get_Gene();
  my $translation = $transcript->translation();

  my ( $biotype_display, $transcript_biotype );
  {
    no warnings 'uninitialized';
    $biotype_display =
      $vegadb ? $gene->status . '_' . $gene->biotype : $gene->biotype;
    $transcript_biotype =
      $vegadb ? $transcript->status . '_' . $transcript->biotype :
      $transcript->biotype;
  }

  # Find selenocysteine
  my $has_selenocysteine = 0;
  my $selenocysteines = [];
  if($translation) {
    $selenocysteines = $translation->get_all_selenocysteine_SeqEdits();
    $has_selenocysteine = 1 if @{$selenocysteines};
  }

  #Print Transcript summary
  print $fh sprintf(qq{%s\t%s\ttranscript\t%d\t%d\t.\t%s\t.\t}, 
        $idstr, $transcript->source, 
        ($transcript->start()+$sliceoffset), ($transcript->end()+$sliceoffset),
        ($strand_conversion{$transcript->strand}));
  $self->_print_attribs($gene, $biotype_display, $transcript, $transcript_biotype, 0, 'transcript', undef, undef, $has_selenocysteine);
  print $fh "\n";

  #Process any selenocystines we may have
  if($has_selenocysteine) {
    foreach my $edit (@{$selenocysteines}) {
      my $edit_start = $edit->start();
      my @projections = $transcript->pep2genomic($edit_start, $edit_start);
      foreach my $projection (@projections) {
        my $strand = $strand_conversion{$projection->strand()};
        my $start = $projection->start();
        my $end = $projection->end();
        print $fh sprintf(qq{%s\t%s\tSelenocysteine\t%d\t%d\t.\t%s\t.\t}, 
          $idstr, $transcript->source, ($start+$sliceoffset), ($end+$sliceoffset), $strand);
        $self->_print_attribs($gene, $biotype_display, $transcript, $transcript_biotype, 0, 'Selenocystine', undef, undef, $has_selenocysteine);
        print $fh "\n";
      }
    }
  }

  my @translateable_exons;
  @translateable_exons = @{ $transcript->get_all_translateable_Exons }
    if $translation;
    
  my ( $count, $intrans, $instop ) = ( 1, 0, 0 );

  foreach my $exon ( @{ $transcript->get_all_Exons } ) {
    my $strand = $strand_conversion{ $exon->strand };

    print $fh
# Column 1 - seqname, the name of the sequence/chromosome the feature is on. Landmark for start below
      $idstr . "\t" .
      # Column 2 - source
      $transcript->source. "\t" .
      # Column 3 - feature type name
      'exon' . "\t" .
      # Column 4 - start, the start coordinate of the feature
      ( $exon->start + $sliceoffset ) . "\t" .
    # Column 5 - end, coordinates (absolute) for the end of this feature
      ( $exon->end + $sliceoffset ) . "\t" .
      # Column 6 - score, for variations only
      "." . "\t" .
      # Column 7 - strand, forward (+) or reverse (-)
      $strand . "\t" .
# Column 8 - frame (reading phase), what base of this feature is the first base of a codon
      "." . "\t";
# Column 9 - attribute, a ;-separated list of key-value pairs (additional feature info)
    $self->_print_attribs( $gene, $biotype_display, $transcript, $transcript_biotype,
                           $count, 'exon', $exon, $vegadb, $has_selenocysteine );
    print $fh "\n";

    $intrans = 1
      if $translation and
      $exon == $translation->start_Exon;

    my $last_used_coding_exon;
    if ($intrans) {
      # print the CDS of this exon

      my $cdsexon = shift @translateable_exons;
      $last_used_coding_exon = $cdsexon;
   #
   # Here is computing the value of the GTF frame (taking into
   # account the Ensembl convention), but it's misleadingly called phase
   #
      my $phase = $cdsexon->phase;
      if ( $cdsexon->phase == 1 ) {
        $phase = 2;
      }
      elsif ( $cdsexon->phase == 2 ) {
        $phase = 1;
      }
      elsif ( $cdsexon->phase == -1 ) {
        $phase = 0;
      }

      my $exon_start = $cdsexon->start;
      my $exon_end   = $cdsexon->end;
      if ( $translation &&
           $hasend &&
           ( $exon->end >= $endcs[0]->start &&
             $exon->start <= $endcs[0]->end ) )
      {

        if ( $cdsexon->strand == 1 ) {
          $exon_end = $cdsexon->end - $endcs[0]->length;
        }
        else {
          $exon_start = $cdsexon->start + $endcs[0]->length;
        }
      }

      # Before logic was the first 3 conditions. Added the length of
      # exon is less than 3 as the GTF dumper missed off
      # 1bp CDS exons. Seems that $instop was set to true meaning
      # we bailed out. The other conditions are fine. Not sure
      # if this test is the right one to do but it does work
      if ( $exon_start <= $cdsexon->end &&
           $exon_end >= $cdsexon->start &&
           (!$instop || $cdsexon->length() < 3) )
      {
        print $fh $idstr . "\t" . $transcript->source .
          "\t" . 'CDS' . "\t" . ( $exon_start + $sliceoffset ) .
          "\t" . ( $exon_end + $sliceoffset ) .
          "\t" . "." . "\t" . $strand . "\t" . $phase . "\t";
        $self->_print_attribs( $gene, $biotype_display, $transcript, $transcript_biotype,
                               $count, 'CDS', undef, undef, $has_selenocysteine );
        print $fh "\n";
      }
    } ## end if ($intrans)

    # The alternative is that this region could be described as UTR but *only* if
    # the transcript was coding.


    if ( $translation &&
         $exon == $translation->start_Exon &&
         $hasstart )
    {
      my $tmpcnt = $count;
      foreach my $startc (@startcs) {
        # here we should check the start codon covers 3 bases
        print $fh $idstr . "\t" . $transcript->source . "\t" .
          'start_codon' . "\t" . ( $startc->start + $sliceoffset ) .
          "\t" . ( $startc->end + $sliceoffset ) .
          "\t" . "." . "\t" . $strand . "\t" . $startc->phase . "\t";

        $self->_print_attribs( $gene, $biotype_display, $transcript, $transcript_biotype,
                               $tmpcnt++, 'start_codon', undef, undef, $has_selenocysteine );
        print $fh "\n";
      }
    }
    if ( $translation &&
         ( $exon == $translation->end_Exon ) )
    {
      if ($hasend) {
        my $tmpcnt = $count - $#endcs;

        foreach my $endc (@endcs) {
          # here we should check the stop codon covers 3 bases
          print $fh $idstr . "\t" . $transcript->source . "\t" .
            'stop_codon' . "\t" . ( $endc->start + $sliceoffset ) .
            "\t" . ( $endc->end + $sliceoffset ) .
            "\t" . "." . "\t" . $strand . "\t" . $endc->phase . "\t";

          $self->_print_attribs( $gene, $biotype_display, $transcript, $transcript_biotype,
                                 $tmpcnt++, 'stop_codon', undef, undef, $has_selenocysteine );
          print $fh "\n";
        }
      }
      $intrans = 0;
    }

    if ( scalar(@endcs) &&
         ( $exon->end >= $endcs[0]->start &&
           $exon->start <= $endcs[0]->end ) )
    {
      $instop = 1;
    }

    $count++;
  } ## end foreach my $exon ( @{ $transcript...})

  my $utrs = $self->get_all_UTR_features($transcript);
  foreach my $utr (@{$utrs}) {
    my $strand = $strand_conversion{$utr->strand()};
    print $fh sprintf(qq{%s\t%s\tUTR\t%d\t%d\t.\t%s\t.\t}, 
        $idstr, $transcript->source, ($utr->start()+$sliceoffset), ($utr->end+$sliceoffset), $strand);
    $self->_print_attribs($gene, $biotype_display, $transcript, $transcript_biotype, 0, 'UTR', undef, undef, $has_selenocysteine);
    print $fh "\n";
  }

  # my $three_prime_utr = $transcript->three_prime_utr_Feature();
  # if($three_prime_utr) {
  #   my $strand = $strand_conversion{$three_prime_utr->strand()};
  #   print $fh sprintf(qq{%s\t%s\tUTR\t%d\t%d\t.\t%s\t.\t}, 
  #       $idstr, $transcript_biotype, ($three_prime_utr->start()+$sliceoffset), ($three_prime_utr->end+$sliceoffset), $strand);
  #   $self->_print_attribs($gene, $biotype_display, $transcript, 0, 'UTR', undef, undef, $has_selenocysteine);
  #   print $fh "\n";
  # }

  return;
} ## end sub print_feature

=head2 _print_attribs

    Arg []     : 
    Example    : 
    Description: 
    Returntype : None

=cut

sub _print_attribs {
  my ( $self, $gene, $gene_biotype, $transcript, $trans_biotype, $count, $type, $exon,
       $vegadb, $has_selenocysteine )
    = @_;

  my ( $gene_name, $gene_source, $gene_version );
  $gene_name = $gene->external_name;
  $gene_name =~ s/^[A-Z]{1,3}:// if $vegadb;
  $gene_source = $gene->source;
  $gene_version = $gene->version;

  my ( $trans_name, $trans_source, $trans_version );
  $trans_name = $transcript->external_name;
  $trans_name =~ s/^[A-Z]{1,3}:// if $vegadb;
  $trans_source = $transcript->source;
  $trans_version = $transcript->version;

  my $fh = $self->{'filehandle'};

  print $fh "gene_id \"" . get_id_from_obj($gene) ."\";";
  print $fh " gene_version \"" . $gene_version . "\";" if $gene_version;
  if($type ne 'gene') {
    print $fh " transcript_id \"" . get_id_from_obj($transcript) . "\";";
    print $fh " transcript_version \"" . $trans_version . "\";" if $trans_version;
    print $fh " exon_number \"$count\";" if $count > 0;
  }
  print $fh " gene_name \"" . $gene_name . "\";"     if ($gene_name);
  print $fh " gene_source \"" . $gene_source . "\";" if ($gene_source);
  print $fh " gene_biotype \"" . $gene_biotype . "\";";

  if($type ne 'gene') {
    print $fh " transcript_name \"" . $trans_name . "\";"
      if ($trans_name);
    print $fh " transcript_source \"" . $trans_source . "\";"
      if ($trans_source);
    print $fh " transcript_biotype \"" . $trans_biotype . "\";"
      if ($trans_biotype);
    my $ccds_entries = $transcript->get_all_DBEntries('CCDS');
    if(@{$ccds_entries}) {
      print $fh qq{ tag "CCDS";};
      foreach my $ccds (sort { $a->primary_id() cmp $b->primary_id() } @{$ccds_entries}) {
        my $primary_ccds_id = $ccds->primary_id();
        print $fh qq{ ccds_id "$primary_ccds_id";};
      }
    }
  }

  if ( $type eq 'CDS' ) {
    print $fh ' protein_id "' .
      get_id_from_obj( $transcript->translation ) . '";';
    print $fh ' protein_version "' .
      $transcript->translation->version . "\";" if $transcript->translation->version;
  }
  if ($exon) {
    printf $fh ' exon_id "%s";', get_id_from_obj($exon);
    print $fh ' exon_version "' .
      $exon->version . "\";" if $exon->version;
  }

  if($has_selenocysteine) {
    print $fh qq{ tag "seleno";};
  }

  if($transcript && $transcript->isa('Bio::EnsEMBL::Transcript')) {
    foreach my $tag (qw/cds_end_NF cds_start_NF mRNA_end_NF mRNA_start_NF/) {
      my $attributes = $transcript->get_all_Attributes($tag);
      if(@{$attributes}) {
        print $fh qq{ tag "${tag}";};
      }
    }
  }

  return;
} ## end sub _print_attribs

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
  my ( $self, $trans ) = @_;
  defined $trans or throw("Transcript object not defined");

  return ( () ) unless $trans->translation;

  my @pepgencoords = $trans->pep2genomic( 1, 1 );

  # cdna can come padded these days so allow gap at the start
  if ( $pepgencoords[0]->isa('Bio::EnsEMBL::Mapper::Gap') ) {
    shift @pepgencoords;
  }

  if ( scalar(@pepgencoords) > 3 ) {
    throw( sprintf "Pep start for transcript %s does not map cleanly",
           $trans->display_id );
  }

  unless ( $pepgencoords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate') ) {
    throw( sprintf "Pep start for transcript %s maps to gap",
           $trans->display_id );
  }
  unless ( $pepgencoords[$#pepgencoords]
           ->isa('Bio::EnsEMBL::Mapper::Coordinate') )
  {
    throw( sprintf "Pep start (end of) for transcript %s maps to gap",
           $trans->display_id );
  }

  my @translateable = @{ $trans->get_all_translateable_Exons };
  my @startc_feat;
  my $phase = 0;
  foreach my $pepgencoord (@pepgencoords) {
    push @startc_feat,
      new Bio::EnsEMBL::SeqFeature( -seqname     => $trans->stable_id,
                                    -source_tag  => 'starttrans',
                                    -primary_tag => 'similarity',
                                    -start       => $pepgencoord->start,
                                    -end         => $pepgencoord->end,
                                    -phase       => $phase,
                                    -strand => $translateable[0]->strand
      );
    $phase = 3 - ( $pepgencoord->end - $pepgencoord->start + 1 );
  }
  if ( $translateable[0]->strand == 1 ) {
    @startc_feat = sort { $a->start <=> $b->start } @startc_feat;
  }
  else {
    @startc_feat = sort { $b->start <=> $a->start } @startc_feat;
  }

  return @startc_feat;

} ## end sub _make_start_codon_features

=head2 _make_stop_codon_features

    Arg [1]    : Bio::EnsEMBL::Transcript
    Example    : 
    Description: 
    Returntype : Array

=cut

sub _make_stop_codon_features {
  my ( $self, $trans ) = @_;
  defined $trans or throw("Transcript object not defined");

  return ( () ) unless $trans->translation;

  my @translateable = @{ $trans->get_all_translateable_Exons };

  my $cdna_endpos = $trans->cdna_coding_end;

  my @pepgencoords =
    $trans->cdna2genomic( $cdna_endpos - 2, $cdna_endpos );

  if ( scalar(@pepgencoords) > 3 ) {
    throw( sprintf "Pep end for transcript %s does not map cleanly",
           $trans->display_id );
  }
  unless ( $pepgencoords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate') ) {
    throw( sprintf "Pep end for transcript %s maps to gap",
           $trans->display_id );
  }
  unless ( $pepgencoords[$#pepgencoords]
           ->isa('Bio::EnsEMBL::Mapper::Coordinate') )
  {
    throw( sprintf "Pep end (end of) for transcript %s maps to gap",
           $trans->display_id );
  }

  my @stopc_feat;
  my $phase = 0;
  foreach my $pepgencoord (@pepgencoords) {
    push @stopc_feat,
      new Bio::EnsEMBL::SeqFeature( -seqname     => $trans->display_id,
                                    -source_tag  => 'endtrans',
                                    -primary_tag => 'similarity',
                                    -start       => $pepgencoord->start,
                                    -end         => $pepgencoord->end,
                                    -phase       => $phase,
                                    -strand => $translateable[0]->strand
      );
    $phase = 3 - ( $pepgencoord->end - $pepgencoord->start + 1 );
  }

  if ( $translateable[0]->strand == 1 ) {
    @stopc_feat = sort { $a->start <=> $b->start } @stopc_feat;
  }
  else {
    @stopc_feat = sort { $b->start <=> $a->start } @stopc_feat;
  }

  return @stopc_feat;

} ## end sub _make_stop_codon_features

=head2 _check_start_and_stop

    Arg [1]    : Bio::EnsEMBL::Transcript
    Example    : 
    Description: 
    Returntype : Array

=cut

sub _check_start_and_stop {
  my ( $self,, $trans ) = @_;

  return ( 0, 0 ) unless defined $trans->translation;
  my ( $has_start, $has_end );

  # transcript could be annotated has having incomplete
  # CDS at either 5', 3' end or both
  my @attrib = @{ $trans->get_all_Attributes('cds_start_NF') };
  $has_start =
    ( scalar @attrib == 1 and $attrib[0]->value() == 1 ) ? 0 : 1;
  @attrib = @{ $trans->get_all_Attributes('cds_end_NF') };
  $has_end =
    ( scalar @attrib == 1 and $attrib[0]->value() == 1 ) ? 0 : 1;
  return ( 0, 0 ) unless $has_start and $has_end;

#
# even if the transcript is not annotated with incomplete start/end
# CDS, we need to test whether the extracted start/stop codons are valid
#
# use translateable_seq (CDS) instead of spliced_seq (CDNA) which is
# not padded for non-triplet issues
#
# $has_start = $has_end = 1;
# my $coding_start = $trans->cdna_coding_start;
# my $coding_end = $trans->cdna_coding_end;
# my $cdna_seq = uc($trans->spliced_seq);
# my $startseq = substr($cdna_seq, $coding_start-1, 3);
# my $endseq = substr($cdna_seq, $coding_end-3, 3);
#
  my $cds_seq  = uc( $trans->translateable_seq );
  my $startseq = substr( $cds_seq, 0, 3 );
  my $endseq   = substr( $cds_seq, -3 );

# reimplemented since there are alternatively valid codon tables
# $has_start = 0  if ($startseq ne "ATG");
# $has_end = 0 if ($endseq ne "TAG" && $endseq ne "TGA" && $endseq ne "TAA");

  my ($attrib) =
    @{ $trans->slice()->get_all_Attributes('codon_table') };

  my $codon_table_id;
  $codon_table_id = $attrib->value() if defined $attrib;
  $codon_table_id ||= 1;    # default codon table (vertebrate)
  my $codon_table =
    Bio::Tools::CodonTable->new( -id => $codon_table_id );

  $has_start = 0 unless $codon_table->is_start_codon($startseq);
  $has_end   = 0 unless $codon_table->is_ter_codon($endseq);

  return ( $has_start, $has_end );

} ## end sub _check_start_and_stop

sub get_all_UTR_features {
  my ($self, $transcript) = @_;
  my $translation = $transcript->translation();
  return [] if ! $translation;
  
  my @utrs;

  my $cdna_coding_start = $transcript->cdna_coding_start();
  # if it is greater than 1 then it must have UTR
  if($cdna_coding_start > 1) {
    my @projections = $transcript->cdna2genomic(1, ($cdna_coding_start-1));
    foreach my $projection (@projections) {
      next if $projection->isa('Bio::EnsEMBL::Mapper::Gap');
      my $f = Bio::EnsEMBL::Feature->new(
        -START => $projection->start, 
        -END => $projection->end,
        -STRAND => $projection->strand,
      );
      push(@utrs, $f);
    }
  }

  my $cdna_coding_end = $transcript->cdna_coding_end();
  if($cdna_coding_end < $transcript->length()) {
    my @projections = $transcript->cdna2genomic(($cdna_coding_end+1), $transcript->length());
    foreach my $projection (@projections) {
      next if $projection->isa('Bio::EnsEMBL::Mapper::Gap');
      my $f = Bio::EnsEMBL::Feature->new(
        -START => $projection->start, 
        -END => $projection->end,
        -STRAND => $projection->strand,
      );
      push(@utrs, $f);
    } 
  }

  return \@utrs;
}

1;
