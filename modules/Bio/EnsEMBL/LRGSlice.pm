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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

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

  my $max;
  my $min;
  my $chrom;
  my $strand;
  my %segments;

  foreach my $segment (@{$self->project('chromosome')}) {
    my $from_start = $segment->from_start();
    my $from_end   = $segment->from_end();
    my $slice      = $segment->to_Slice;
    
    my $seq_region_id = $slice->adaptor->get_seq_region_id($slice);
    my $to_name    = $slice->seq_region_name();
    
    next if (!$seq_region_id);
    
    if (!$segments{$seq_region_id}) {
      $segments{$seq_region_id}{'chr_name'} = $to_name;
      $segments{$seq_region_id}{'strand'}   = $slice->strand;
      $max=-99999999999;
      $min=9999999999;
    } else {
      $max = $segments{$seq_region_id}{'max'};
      $min = $segments{$seq_region_id}{'min'};
    }

    my $to_start = $slice->start();
    my $to_end   = $slice->end();
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
    
    $segments{$seq_region_id}{'max'} = $max;
    $segments{$seq_region_id}{'min'} = $min;  
  }
  
  my @seq_region_ids = sort {$a <=> $b} keys (%segments);
  
  # Keep the chromosome mapping instead of the patch, if possible
  if (scalar @seq_region_ids) {
    my $seq_id = $seq_region_ids[0];
    $chrom  = $segments{$seq_id}{'chr_name'};
    $max    = $segments{$seq_id}{'max'};
    $min    = $segments{$seq_id}{'min'};
    $strand = $segments{$seq_id}{'strand'};
  }

  if(!defined($chrom)){
    warn "Could not project to chromosome for ".$self->name." ??\n";
    return undef;
  }
  my $chrom_slice = $self->adaptor->fetch_by_region("chromosome",$chrom, $min, $max, $strand);
  $self->{_chrom_slice} = $chrom_slice;
  return $self->{_chrom_slice};
}

sub DESTROY{
}

sub get_all_differences{
  my $self = shift;
  
  my @results;
  
  # get seq_region_attrib diffs (always same-length substitutions)
  ################################################################
  
  my $sth = $self->adaptor->prepare(qq{
    SELECT sra.value
    FROM seq_region_attrib sra, attrib_type at
    WHERE at.code = '_rna_edit'
    AND at.attrib_type_id = sra.attrib_type_id
    AND sra.seq_region_id = ?
  });
  
  $sth->execute($self->get_seq_region_id);
  
  my $edit_string;
  
  $sth->bind_columns(\$edit_string);
  
  while($sth->fetch()) {
    my ($start, $end, $edit) = split " ", $edit_string;
    
    my $slice = $self->sub_Slice($start, $end);
    my $chr_proj = $slice->project("chromosome");
    my $ref_seq = '-';
    if(scalar @$chr_proj == 1) {
      $ref_seq = $chr_proj->[0]->[2]->seq;
    }
    
    
    my $diff = {
      'start' => $start,
      'end'   => $end,
      'type'  => 'substitution',
      'seq'   => $edit,
      'ref'   => $ref_seq,
    };
    
    push @results, $diff;
  }
  
  # get more complex differences via projections
  ##############################################
  
  # project the LRG slice to contig coordinates
  my @segs = @{$self->project("contig")};
  
  # if the LRG projects into more than one segment
  if(scalar @segs > 1) {
    
    my ($prev_end, $prev_chr_start, $prev_chr_end, $prev_was_chr);
    
    foreach my $seg(@segs) {
      
      # is this a novel LRG contig, or does it project to a chromosome?
      my @chr_proj = @{$seg->[2]->project("chromosome")};
      
      # if it is a normal contig
      if(scalar @chr_proj) {
        
        # check if there has been a deletion in LRG
        if($prev_was_chr && $prev_end == $seg->[0] - 1) {
          
          # check it's not just a break in contigs
          unless(
             ($chr_proj[0]->[2]->strand != $self->strand && $prev_chr_start == $chr_proj[0]->[2]->end + 1) ||
             ($chr_proj[0]->[2]->strand != $self->strand && $prev_chr_end == $chr_proj[0]->[2]->start - 1)
          ) {
            
            # now get deleted slice coords, depends on the strand rel to LRG
            my ($s, $e);
            
            # opposite strand
            if($chr_proj[0]->[2]->strand != $self->strand) {
              ($s, $e) = ($prev_chr_start - 1, $chr_proj[0]->[2]->end + 1);
            }
            
            # same strand
            else {
              ($s, $e) = ($prev_chr_end + 1, $chr_proj[0]->[2]->start - 1);
            }
            
            if($s > $e) {
              warn "Oops, trying to create a slice from $s to $e (could have been ", $prev_chr_start - 1, "-", $chr_proj[0]->[2]->end + 1, " or ", $prev_chr_end + 1, "-", $chr_proj[0]->[2]->start - 1, ")";
            }
            
            else {
              # get a slice representing the sequence that was deleted
              my $deleted_slice = $self->adaptor->fetch_by_region("chromosome", $chr_proj[0]->[2]->seq_region_name, $s, $e, $chr_proj[0]->[2]->strand);
              
              my $diff = {
                'start' => $seg->[0],
                'end'   => $prev_end,
                'type'  => 'deletion',
                'seq'   => '-',
                'ref'   => $deleted_slice->seq." ".$deleted_slice->start.'-'.$deleted_slice->end,
              };
              
              push @results, $diff;
            }
          }
        }
        
        $prev_was_chr = 1;
        
        $prev_chr_start = $chr_proj[0]->[2]->start;
        $prev_chr_end = $chr_proj[0]->[2]->end;
      }
      
      # if it is an LRG made-up contig for an insertion
      else {
        $prev_was_chr = 0;
        
        my $diff = {
          'start' => $seg->[0],
          'end'   => $seg->[1],
          'type'  => 'insertion',
          'seq'   => substr($self->seq, $seg->[0] - 1, $seg->[1] - $seg->[0] + 1),
          'ref'   => '-',
        };
        
        push @results, $diff;
      }
      
      $prev_end = $seg->[1];
    }
  }
  
  # return results sorted by start, then end position
  return [sort {$a->{start} <=> $b->{start} || $a->{end} <=> $b->{end}} @results];
}

1;
