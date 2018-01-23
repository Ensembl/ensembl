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

=head1 NAME

BEDSerializer - Feature to BED converter

=head1 SYNOPSIS

use Bio::EnsEMBL::Utils::IO::BEDSerializer;

=head1 DESCRIPTION

Subclass of Serializer that can turn a feature into a line for the BED format.

=cut

package Bio::EnsEMBL::Utils::IO::BEDSerializer;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception;
use URI::Escape;
use Bio::EnsEMBL::Utils::IO::FeatureSerializer;
use Bio::EnsEMBL::Utils::Scalar qw/check_ref/;

use base qw(Bio::EnsEMBL::Utils::IO::FeatureSerializer);

=head2 new

    Constructor
    Arg [1]    : Optional File handle
    Arg [2]    : Default source of the features. Defaults to .
    Arg [3]    : RGB colour to emit. Defaults to black (0,0,0)

    Returntype : Bio::EnsEMBL::Utils::IO::BEDSerializer

=cut

sub new {
    my $class = shift;
    my $self = {
        filehandle => shift,
        source => shift,
        rgb => shift
    };
    bless $self, $class;

    if (!defined ($self->{'filehandle'})) {
        # no file handle, let the handle point to a copy of STDOUT instead
        open $self->{'filehandle'}, ">&STDOUT";
        $self->{'stdout'} = 1;
    }

    $self->{rgb} = '0,0,0' if ! defined $self->{rgb};

    return $self;
}


=head2 print_feature

    Arg [1]    : Bio::EnsEMBL::Feature, subclass or related pseudo-feature
    Example    : $reporter->print_feature($feature,$slice_start_coordinate,"X")
    Description: Asks a feature for its summary, and generates a BED
                 compliant entry to hand back again
    Returntype : none
=cut

sub print_feature {
  my ($self, $feature, $cache) = @_;
  return unless $feature;
  my $bed_array;
  if($feature->isa('Bio::EnsEMBL::Transcript')) {
    $bed_array = $self->_write_Transcript($feature, $cache);
  }
  else {
    $bed_array = $self->_feature_to_bed_array($feature, $cache);
  }
  my $bed_line = join("\t", @{$bed_array});
  my $fh = $self->{'filehandle'};
  print $fh $bed_line, "\n";
  return 1;
}

sub rgb {
  my ($self, $rgb) = @_;
  $self->{rgb} = $rgb if defined $rgb;
  return $self->{rgb};
}

sub _write_Transcript {
  my ($self, $transcript, $cache) = @_;

  # Not liking this. If we are in this situation we need to re-fetch the transcript
  # just so the thing ends up on the right Slice!
  my $new_transcript = $transcript->transfer($transcript->slice()->seq_region_Slice());
  $new_transcript->get_all_Exons(); # force exon loading
  my $bed_array = $self->_feature_to_bed_array($transcript, $cache);
  my $bed_genomic_start = $bed_array->[1]; #remember this is in 0 coords
  # 12 column BED
  my ($coding_start, $coding_end, $exon_starts_string, $exon_lengths_string, $exon_count) = (0,0,q{},q{},0);
  my $rgb = $self->rgb();
  # genePred BED extensions
  my ($second_name, $cds_start_status, $cds_end_status, $exon_frames, $type, $gene_name, $second_gene_name, $gene_type) = (q{},q{none},q{none},q{},q{},q{},q{},q{});
  
  # Set the remaining transcript attributes
  $type = $transcript->biotype();
  $second_name = $transcript->external_name() || q{none};

  # If we have a translation then we do some maths to calc the start of 
  # the thick sections. Otherwise we must have a ncRNA or pseudogene
  # and that thick section is just set to the transcript's end
  if($new_transcript->translation()) {
    my ($cdna_start, $cdna_end) = ($new_transcript->cdna_coding_start(), $new_transcript->cdna_coding_end);
    if($new_transcript->strand() == -1) {
      ($cdna_start, $cdna_end) = ($cdna_end, $cdna_start);
    }
    # Rules are if it's got a coding start we will use it; if not we use the cDNA
    $coding_start = $self->_cdna_to_genome($new_transcript, $cdna_start);
    $coding_start--; # convert to 0 based coords

    #Same again but for the end
    $coding_end = $self->_cdna_to_genome($new_transcript, $cdna_end);

    #Also figure out complete 5' and 3' CDS tags
    my ($five_prime_nc, $three_prime_nc) =  map { $_->[0] } 
                                            grep { defined $_ && @{$_} } 
                                            map { $transcript->get_all_Attributes($_) } 
                                            qw/cds_start_NF cds_end_NF/; 

    
    $cds_start_status = ($five_prime_nc) ? 'incmpl' : 'cmpl';
    $cds_end_status = ($three_prime_nc) ? 'incmpl' : 'cmpl';
    # reverse if we are on the negative strand
    if($transcript->strand() == -1) {
      ($cds_start_status, $cds_end_status) = ($cds_end_status, $cds_start_status);
    }
  }
  else {
    # When we represent non-coding transcripts we set the coding start/end (thickStart/thickEnd)
    # in a BED file to the start of the transcript region (lowest coordinate).
    $coding_start = $bed_genomic_start;
    $coding_start--;
    $coding_end = $coding_start;
  }

  # Now for the interesting bit. Exons are given relative to the bed start
  # so we need to calculate the offset. Exons are also sorted by start otherwise 
  # offset calcs are wrong
  foreach my $exon (sort { $a->seq_region_start() <=> $b->seq_region_start() } @{$new_transcript->get_all_Exons()}) {
    my $exon_start = $exon->seq_region_start();
    $exon_start--; #move into 0 coords
    my $offset = $exon_start - $bed_genomic_start; # just have to minus current start from the genomic start
    $exon_starts_string .= $offset.',';
    $exon_lengths_string .= $exon->length().',';
    
    # We have to re-interpret a phase -1 as 0 if we are on a coding exon. Otherwise we just 
    # leave it (non-coding exons are always -1 it seems from UCSC's examples)
    my $phase = $exon->phase();
    my $exon_coding_start = $exon->coding_region_start($transcript);
    if(defined $exon_coding_start) {
      $phase = 0 if $phase == -1;
    }
    $exon_frames .= $phase.',';
    $exon_count++;
  }

  # Get the gene and populate its fields
  my $gene = $transcript->get_Gene();
  $gene_name = $gene->stable_id();
  $second_gene_name = $gene->external_name() || q{none};
  $gene_type = $gene->biotype();

  # Finally recreate 12 column BED format
  push(@{$bed_array}, $coding_start, $coding_end, $rgb, $exon_count, $exon_lengths_string, $exon_starts_string);
  # And then the gene pred
  push(@{$bed_array}, $second_name, $cds_start_status, $cds_end_status, $exon_frames, $type, $gene_name, $second_gene_name, $gene_type);

  return $bed_array;
}

sub _cdna_to_genome {
  my ($self, $transcript, $coord) = @_;
  my @mapped = $transcript->cdna2genomic($coord, $coord);
  my $genomic_coord = $mapped[0]->start();
  return $genomic_coord;
}

=head2 _feature_to_bed_array

  Arg [1]     : Bio::EnsEMBL::Feature
                The Feature to encode
  Arg [2]     : HashRef
                Cache of retrieved Slice names.
  Description : Takes a feature and returns an extended BED record consumed up to field
                6 of the format (strand). Score is left as 0
  Returntype  : ArrayRef of fields encoded as
                [chr, start, end, display, 0, strand]

=cut

sub _feature_to_bed_array {
  my ($self, $feature, $cache) = @_;
  my $chr_name = $self->_feature_to_UCSC_name($feature, $cache);
  my $start = $feature->seq_region_start() - 1;
  my $end = $feature->seq_region_end();
  my $strand = ($feature->seq_region_strand() == -1) ? '-' : '+'; 
  my $display_id = $feature->display_id();
  my $score = 1000;
  return [ $chr_name, $start, $end, $display_id, $score, $strand ];
}

=head2 _feature_to_UCSC_name

  Arg [1]     : Bio::EnsEMBL::Feature
                Feature object whose seq_region_name must be converted
  Arg [2]     : HashRef
                Cache of retrieved names. Here because we are unaware of the
                lifetime of a View in Catalyst (is it persistent between requests)
  Description : Attempts to figure out what UCSC calls this Slice. First we
                consult the synonyms attached to the Slice, then ask if it is a chromsome
                (accounting for MT -> chrM) and adding a chr to the name if it was a reference
                Slice. If it was none of these the name is just passed through
  Returntype  : String of the UCSC name

=cut

sub _feature_to_UCSC_name {
  my ($self, $feature, $cache) = @_;
  my $seq_region_name = $feature->seq_region_name();

  #Return if the name was already defined (we assume we work within a single species)
  return $cache->{$seq_region_name} if exists $cache->{$seq_region_name};

  if(! $feature->can('slice')) {
    return $cache->{$seq_region_name} = $seq_region_name;
  }
  my $slice = $feature->slice();
  my $ucsc_name;
  my $has_adaptor = ($slice->adaptor()) ? 1 : 0;
  if($has_adaptor) { # if it's got an adaptor we can lookup synonyms
    my $synonyms = $slice->get_all_synonyms('UCSC');
    if(@{$synonyms}) {
      $ucsc_name = $synonyms->[0]->name();
    }
  }
  if(! defined $ucsc_name) {
    # if it's a chromosome then we can test a few more things
    if($slice->is_chromosome()) {
      #MT is a special case; it's chrM
      if($seq_region_name eq 'MT' ) {
        $ucsc_name = 'chrM';
      }
      # If it was a ref region add chr onto it (only check if we have an adaptor)
      elsif($has_adaptor && $slice->is_reference()) {
        $ucsc_name = 'chr'.$seq_region_name;
      }
    }
  }
  #Leave it as the seq region name otherwise
  $ucsc_name = $seq_region_name if ! defined $ucsc_name;
  $cache->{$seq_region_name} = $ucsc_name;
  return $ucsc_name;
}

1;
