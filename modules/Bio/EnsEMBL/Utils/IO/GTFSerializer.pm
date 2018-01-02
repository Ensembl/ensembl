
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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
     my $string = "#!genome-build-accession ";
     $string .= "$accession_source:" if $accession_source;
     $string .= "$accession"; 

     print $fh "$string\n";
     #my $accession_source = $mc->single_value_by_key('assembly.web_accession_source');
     #print $fh "#!genome-build-accession ${accession_source}:${accession}\n";
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
  my $fh              = $self->{'filehandle'};
  my $biotype_display = q{};
  if($gene->adaptor()) {
    my $dbname = $gene->adaptor()->dbc()->dbname();
  }

  {
    no warnings 'uninitialized';
    $biotype_display = $gene->biotype;
  }

  # Skip all trans-splicing transcripts, they can't be dumped in GTF
  # format.
  my @transcripts = @{$gene->get_all_Transcripts()};
  my $single_trans_gene = scalar(@transcripts);
  my %trans_splicing_ids;
  foreach my $t (@transcripts) {
    if (map { $_->value } @{ $t->get_all_Attributes('trans_spliced') }) {
      $trans_splicing_ids{$t->stable_id()} = 1;
    }
  }
  if ((values (%trans_splicing_ids)) && $single_trans_gene == 1) {
    return;
  } else {
    print $fh sprintf(qq{%s\t%s\tgene\t%d\t%d\t.\t%s\t.\t}, 
      $idstr, $gene->source, 
      ($gene->start()+$sliceoffset), ($gene->end()+$sliceoffset),
      ($strand_conversion{$gene->strand}));
    $self->_print_attribs($gene, $biotype_display, $gene, $biotype_display, 0, 'gene');
    print $fh "\n";
  }

  # Now print all transcripts
  foreach my $t (@transcripts) {
    next if (exists ($trans_splicing_ids{$t->stable_id()}));
    $self->print_feature($t, $gene);
  }
  return;
}

sub print_Prediction {
  my ($self, $prediction) = @_;

  #Print prediction transcript summary
  my $slice           = $prediction->slice();
  my $idstr           = $slice->seq_region_name;
  my $sliceoffset     = $slice->start - 1;
  my $fh              = $self->{'filehandle'};

  my $source = $prediction->analysis->gff_source || 'ensembl';

  print $fh sprintf(qq{%s\t%s\ttranscript\t%d\t%d\t.\t%s\t.\t},
        $idstr, $source,
        ($prediction->start()+$sliceoffset), ($prediction->end()+$sliceoffset),
        ($strand_conversion{$prediction->strand}));
  print $fh "\n";

  # Now print all transcripts
  foreach my $exon (@{$prediction->get_all_Exons()}) {
    print $fh sprintf(qq{%s\t%s\texon\t%d\t%d\t.\t%s\t.\t},
          $idstr, $source,
          ($exon->start()+$sliceoffset), ($exon->end()+$sliceoffset),
          ($strand_conversion{$exon->strand}));
    print $fh "\n";
  }
  return;
}

=head2 print_feature

    Arg [1]    : Bio::EnsEMBL::Transcript
    Example    : $serializer->print_feature($transcript)
    Description: Generate all the GTF features required for this transcript
    Returntype : none

=cut

sub print_feature {
  my $self       = shift;
  my $transcript = shift;
  my $gene = shift;

  throw( sprintf "Feature is of type %s. Cannot write non transcripts to GTF", ref($transcript) ) 
    unless check_ref( $transcript, "Bio::EnsEMBL::Transcript" );

  # filehandle is inherited
  my $fh = $self->{'filehandle'};

  my $slice       = $transcript->slice();
  my $sliceoffset = $slice->start - 1;
  my $idstr       = $slice->seq_region_name;

  ## Seq edits are needed for sequence, not coordinates
  $transcript->edits_enabled(0);

  # Create start codon features. Multiple due to outcomes of projection, but they will 
  # overlap the first exon
  my @startcodons = $self->_make_start_codon_features($transcript);

  # Stop codons (multiple "features" from projection), are NOT included in the CDS in GTF, so the
  # CDS coordinates get adjusted in the process of creating the stop.
  my @endcodons   = $self->_make_stop_codon_features($transcript);

  # We then verify the positions really are proper stops and starts, by looking for annotation
  # and checking the bases constitute a stop. These flags are then used to decide whether to
  # actually print the features created by the two calls above.
  my ( $hasstart, $hasend ) = $self->_check_start_and_stop($transcript);

  if($transcript->adaptor()) {
    my $dbname = $transcript->adaptor()->dbc()->dbname();
  }
  $gene   ||= $transcript->get_Gene();
  my $translation = $transcript->translation();

  my ( $biotype_display, $transcript_biotype );
  {
    no warnings 'uninitialized';
    $biotype_display = $gene->biotype;
    $transcript_biotype = $transcript->biotype;
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
  $self->_print_attribs($gene, $biotype_display, $transcript, $transcript_biotype, 0, 'transcript', undef, $has_selenocysteine);
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
        $self->_print_attribs($gene, $biotype_display, $transcript, $transcript_biotype, 0, 'Selenocysteine', undef, $has_selenocysteine);
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
                             $count, 'exon', $exon, $has_selenocysteine );
    print $fh "\n";

    $intrans = 1 if $translation and $exon == $translation->start_Exon;

    my $last_used_coding_exon;
    if ($intrans) {
      # print the CDS of this exon

      my $cdsexon = shift @translateable_exons;
      $last_used_coding_exon = $cdsexon;
   #
   # Computing the value of the GTF frame (taking into
   # account the Ensembl convention where it is called phase )
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
      foreach my $endcodon (@endcodons) {
        if ( $translation &&
             $hasend &&
             ( $cdsexon->seq_region_end >= $endcodon->start &&
               $cdsexon->seq_region_start <= $endcodon->end ) )
        {
          if ( $cdsexon->strand == 1 ) {
            $exon_end = $cdsexon->end - $endcodon->length;
          }
          else {
            $exon_start = $cdsexon->start + $endcodon->length;
          }
        }
      }

      # Exons with length less than 3 are included to avoid missing 1bp CDS exons.

      if ( $exon_start <= $cdsexon->end &&
           $exon_end >= $cdsexon->start &&
           (!$instop || $cdsexon->length() < 3) )
      {
        print $fh $idstr . "\t" . $transcript->source .
          "\t" . 'CDS' . "\t" . ( $exon_start + $sliceoffset ) .
          "\t" . ( $exon_end + $sliceoffset ) .
          "\t" . "." . "\t" . $strand . "\t" . $phase . "\t";
        $self->_print_attribs( $gene, $biotype_display, $transcript, $transcript_biotype,
                               $count, 'CDS', undef, $has_selenocysteine );
        print $fh "\n";
      }
    } ## end if ($intrans)

    # The alternative is that this region could be described as UTR but *only* if
    # the transcript was coding.


    if ( $translation 
         && $exon == $translation->start_Exon 
         && $hasstart ) {
      
      my $tmpcnt = $count; # used to order GTF features in the file
      foreach my $startc (@startcodons) {
        # here we should check the start codon covers 3 bases
        print $fh $idstr . "\t" . $transcript->source . "\t" .
          'start_codon' . "\t" . ( $startc->start ) .
          "\t" . ( $startc->end ) .
          "\t" . "." . "\t" . $strand . "\t" . $startc->phase . "\t";

        $self->_print_attribs( $gene, $biotype_display, $transcript, $transcript_biotype,
                               $tmpcnt++, 'start_codon', undef, $has_selenocysteine );
        print $fh "\n";
      }
    }
    if ( $translation && ( $exon == $translation->end_Exon ) ) {
      if ($hasend) {
        my $tmpcnt = $count - $#endcodons;

        foreach my $endc (@endcodons) {
          # here we should check the stop codon covers 3 bases
          print $fh $idstr . "\t" . $transcript->source . "\t" .
            'stop_codon' . "\t" . ( $endc->start ) .
            "\t" . ( $endc->end ) .
            "\t" . "." . "\t" . $strand . "\t" . $endc->phase . "\t";

          $self->_print_attribs( $gene, $biotype_display, $transcript, $transcript_biotype,
                                 $tmpcnt++, 'stop_codon', undef, $has_selenocysteine );
          print $fh "\n";
        }
      }
      $intrans = 0;
    }

    if ( scalar(@endcodons) &&
         ( $exon->end >= $endcodons[0]->start &&
           $exon->start <= $endcodons[0]->end ) )
    {
      $instop = 1;
    }

    $count++;
  } ## end foreach my $exon ( @{ $transcript...})


  my $utrs = $transcript->get_all_five_prime_UTRs();
  push @$utrs, @{$transcript->get_all_three_prime_UTRs()};
  foreach my $utr (@{$utrs}) {
    my $strand = $strand_conversion{$utr->strand()};
    print $fh sprintf(qq{%s\t%s\t%s\t%d\t%d\t.\t%s\t.\t}, 
        $idstr, $transcript->source, $utr->type, ($utr->seq_region_start()), ($utr->seq_region_end()), $strand);
    $self->_print_attribs($gene, $biotype_display, $transcript, $transcript_biotype, 0, 'UTR', undef, $has_selenocysteine);
    print $fh "\n";
  }

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
       $has_selenocysteine )
    = @_;

  my ( $gene_name, $gene_source, $gene_version );
  $gene_name = $gene->external_name;
  $gene_source = $gene->source;
  $gene_version = $gene->version;

  my ( $trans_name, $trans_source, $trans_version );
  $trans_name = $transcript->external_name;
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
  if ($gene->havana_gene()) {
    print $fh " havana_gene \"" . $gene->havana_gene->display_id() . "\";";
    print $fh " havana_gene_version \"" . $gene->havana_gene->version() . "\";";
  }

  #add projection parent
  my $proj_parent_attributes = $gene->get_all_Attributes("proj_parent_g");
    if (@{$proj_parent_attributes}) {
      my $value = $proj_parent_attributes->[0]->value;
      print $fh qq{ projection_parent_gene "${value}";};
    }

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
    if ($transcript->havana_transcript()) {
      print $fh " havana_transcript \"" . $transcript->havana_transcript->display_id() . "\";";
      print $fh " havana_transcript_version \"" . $transcript->havana_transcript->version() . "\";";
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
    foreach my $tag (qw/cds_end_NF cds_start_NF mRNA_end_NF mRNA_start_NF gencode_basic/) {
      my $attributes = $transcript->get_all_Attributes($tag);
      if(@{$attributes}) {
        my $value = $tag;
        $value = "basic" if $tag eq "gencode_basic";
        print $fh qq{ tag "${value}";};
      }
    }
    my $attributes = $transcript->get_all_Attributes("TSL");
    if (@{$attributes}) {
      my $value = $attributes->[0]->value;
      $value =~ s/tsl//;
      print $fh qq{ transcript_support_level "${value}";};
    }
    my $proj_parent_attributes = $transcript->get_all_Attributes("proj_parent_t");
    if (@{$proj_parent_attributes}) {
      my $value = $proj_parent_attributes->[0]->value;
      print $fh qq{ projection_parent_transcript "${value}";};
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
    Description: The start codon can lay across several seq_regions, hence an array of them
                 is returned.
    Returntype : Array of Bio::EnsEMBL::SeqFeature representing the start codon

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
  unless ( $pepgencoords[$#pepgencoords]->isa('Bio::EnsEMBL::Mapper::Coordinate') ){
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
    Description: This method is only called after the sequence has already been checked
                 As a result, the assumption of a 3-base stop codon is valid.
                 The stop can lay across several seq_regions, hence the need to generate
                 several stop features.
    Returntype : Array of Bio::EnsEMBL::SeqFeature representing the stop codon

=cut

sub _make_stop_codon_features {
  my ( $self, $trans ) = @_;
  defined $trans or throw("Transcript object not defined");

  return ( () ) unless $trans->translation;

  my @translateable = @{ $trans->get_all_translateable_Exons };
  my @stopc_feat;
  my $cdna_endpos = $trans->cdna_coding_end;
  my @pepgencoords = $trans->cdna2genomic( $cdna_endpos - 2, $cdna_endpos );

  if ( scalar(@pepgencoords) > 3 ) {
    throw( sprintf "Pep end for transcript %s does not map cleanly", $trans->display_id );
  }
  unless ( $pepgencoords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate') ) {
    throw( sprintf "Pep end for transcript %s maps to gap", $trans->display_id );
  }
  unless ( $pepgencoords[$#pepgencoords]->isa('Bio::EnsEMBL::Mapper::Coordinate') ) {
    throw( sprintf "Pep end (end of) for transcript %s maps to gap", $trans->display_id );
  }

  # In event of an irregular end phase, return no stop codon features
  my $last_exon = $translateable[$#translateable];
  if ($last_exon->end_phase != 0 && $last_exon->end_phase != -1) {
    return @stopc_feat;
  }
  my $phase = 0;
  foreach my $pepgencoord (@pepgencoords) {
    push @stopc_feat,
      new Bio::EnsEMBL::SeqFeature( -seqname     => $trans->display_id,
                                    -source_tag  => 'endtrans',
                                    -primary_tag => 'similarity',
                                    -start       => $pepgencoord->start,
                                    -end         => $pepgencoord->end,
                                    -phase       => $phase,
                                    -strand      => $translateable[0]->strand
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
    Description: Perform a two step verification of start and stop codons.
                 The results are used later to decide whether to generate the GTF 
                 start and stop features
    Returntype : Array

=cut

sub _check_start_and_stop {
  my ( $self, $trans ) = @_;

  return ( 0, 0 ) unless defined $trans->translation;
  my ( $has_start, $has_end );

  # transcript could be annotated has having incomplete
  # CDS at either 5', 3' end or both. Reject in the event of 
  # no identifiable start or stop from annotators.
  my @attrib = @{ $trans->get_all_Attributes('cds_start_NF') };
  $has_start =
    ( scalar @attrib == 1 and $attrib[0]->value() == 1 ) ? 0 : 1;
  @attrib = @{ $trans->get_all_Attributes('cds_end_NF') };
  $has_end =
    ( scalar @attrib == 1 and $attrib[0]->value() == 1 ) ? 0 : 1;
  return ( 0, 0 ) unless $has_start or $has_end;

#
# even if the transcript is not annotated with incomplete start/end
# CDS, we need to test whether the extracted start/stop codons are valid
#
# use translateable_seq (CDS) instead of spliced_seq (CDNA) which is
# not padded for non-triplet issues

  my $cds_seq  = uc( $trans->translateable_seq );
  my $startseq = substr( $cds_seq, 0, 3 );

  # check last exon phase to verify stop codon length, before the stop codon is pulled directly
  # from the CDS (which contains the stop codon usually). IG genes and other unusual cases can
  # be missing a stop.
  my @exons = @{ $trans->get_all_translateable_Exons };
  my $last_exon = $exons[$#exons];
  my $phase = $last_exon->end_phase;
  $phase = 0 if $phase == -1;
  my $endseq = substr( $cds_seq, -(3-$phase) );

# reimplemented since there are alternatively valid codon tables

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

1;
