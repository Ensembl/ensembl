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

Bio::EnsEMBL::Upstream - Object that defines an upstream region

=head1 SYNOPSIS

  use Bio::EnsEMBL::Upstream;

  my $upstream = Bio::EnsEMBL::Upstream->new(
    -transcript => $transcript,
    -length     => 2000           # bp
  );

  # Retrieve coordinates of upstream region
  my $upstream_region_start = $upstream->upstart;
  my $upstream_region_end   = $upstream->upend;

  # Retrieve coordinates in 'downstream' first intron
  my $intron_region_start = $upstream->downstart;
  my $intron_region_end   = $upstream->downend;

  # Coordinates are returned in the same scheme as the input transcript.
  # However, the coordinates of an upstream region can be transformed to
  # any other scheme using a slice

  $upstream->transform($slice);

  # Coordinates can be retrieved in scheme in the same manner as the
  # above.

=head1 DESCRIPTION

An object that determines the upstream region of a transcript.  Such a
region is non-coding and ensures that other genes or transcripts are
not present.  Ultimately, these objects can be used to looking for
promoter elements.  To this end, it is also possible to derive a region
downstream of the first exon, within the first intron and where promoter
elements sometimes are found.

=head1 METHODS

=cut

package Bio::EnsEMBL::Upstream;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);

@ISA = qw(Bio::EnsEMBL::Feature);


=head2 new

  Arg [transcript] : (optional) Bio::EnsEMBL::Transcript
  Arg [length]     : (optional) int $length
  Example    : $upstream = Bio::EnsEMBL::Upstream->new(-transcript => $transcript, 
						       -length => 2000);
  Description: Creates a new upstream object
  Returntype : Bio::EnsEMBL::Upstream
  Exceptions : none
  Caller     : Bio::EnsEMBL::Transcript, general
  Status     : Stable

=cut

sub new {
  my ($class, @args) = @_;
  my $self = {};
  
  bless $self, $class;
  
  my ($transcript, 
      $length) = rearrange([qw(TRANSCRIPT
			       LENGTH
			      )],@args);
	
  $self->transcript($transcript) 	if defined $transcript;
  $self->length($length) 		if $length;
  
  return $self
}

=head2 transcript

  Arg        : (optional) Bio::EnsEMBL::Transcript
  Example    : $self->transcript($transcript);
  Description: Getter/setter for transcript object
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : Throws if argument is not undefined 
               or a Bio::EnsEMBL::Transcript
  Caller     : $self->new, $self->_derive_coords, 
               $self->_first_coding_Exon
  Status     : Stable

=cut


sub transcript {
  my $self = shift;
  
  if (@_){
    $self->{_transcript} = shift;
    
    if (defined $self->{_transcript}) {
      throw("Transcript is not a Bio::EnsEMBL::Transcript") 
	if (! $self->{_transcript}->isa("Bio::EnsEMBL::Transcript"));
      $self->_flush_cache;
    }
  }
  
  return $self->{_transcript}
}

=head2 length

  Arg        : (optional) int $length
  Example    : $self->length(2000); # bp
  Description: Getter/setter for upstream region length.
  Returntype : int
  Exceptions : Throws if length is requested before it has been set.
  Caller     : $self->new, $self->_derive_coords
  Status     : Stable

=cut

sub length {
  my $self = shift;
  
  if (@_){
    $self->{_length} = shift;
    $self->_flush_cache;
  }
  
  throw("Region length has not been set.")
    unless $self->{_length};
  
  return $self->{_length}
}

=head2 _flush_cache

  Arg        : none
  Example    : $self->_flush_cache;
  Description: Empties cached coordinates (called when 
	       coordinate scheme or region length has changed).
  Returntype : none
  Exceptions : none
  Caller     : $self->length, $self->transform
  Status     : Stable

=cut

sub _flush_cache {
  my $self = shift;
  
  $self->upstart(undef);
  $self->upend(undef);
  $self->downstart(undef);
  $self->downend(undef);
}

=head2 upstart

  Arg        : none
  Example    : $self->upstart;
  Description: Returns the start coordinate of the region 
               upstream of the transcript.  This coordinate 
               is always the furthest from the translation 
               initiation codon, whereas upend always abutts 
               the translation initiation codon.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub upstart {
  my $self = shift;
  
  if (@_) {
    $self->{_upstart} = shift @_;
    return
  }

  if (! defined $self->{_upstart}) {
    $self->_derive_coords('up');
  }

  return $self->{_upstart}
}

=head2 upend

  Arg        : none
  Example    : $self->upend;
  Description: Returns the end coordinate of the region 
               upstream of the transcript.  This coordinate 
               always always abutts the translation 
               initiation codon, whereas upstart always 
               returns the coorindate furthest from the 
               translation initiation codon.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub upend {
  my $self = shift;
	
  if (@_) {
    $self->{_upend} = shift @_;
    return
  }

  if (! defined $self->{_upend}) {
    $self->_derive_coords('up');
  }

  return $self->{_upend}
}

=head2 downstart

  Arg        : none
  Example    : $self->downstart;
  Description: Returns the start coordinate of the region 
               in the first intron of the transcript.  This 
               coordinate is always closest to the first 
               exon (irregardless of strand).
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub downstart {
  my $self = shift;
	
  if (@_) {
    $self->{_downstart} = shift @_;
    return
  }

  if (! defined $self->{_downstart}) {
    $self->_derive_coords('down');
  }

  return $self->{_downstart}
}

=head2 downend

  Arg        : none
  Example    : $self->downend;
  Description: Returns the end coordinate of the region 
               in the first intron of the transcript.  This 
               coordinate is always furthest from the first 
               exon (irregardless of strand).
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub downend {
  my $self = shift;

  if (@_) {
    $self->{_downend} = shift @_;
    return
  }

  if (! defined $self->{_downend}) {
    $self->_derive_coords('down');
  }

  return $self->{_downend}
}

=head2 transform

  Arg        : 
  Example    : 
  Description: Not yet implemented
  Returntype : 
  Exceptions : 
  Caller     : 
  Status     : At Risk

=cut


# Over-riding inherited class.  As yet unimplemented.

sub transform {
  my $self = shift;

  throw("No transform method implemented for " . $self);
}

=head2 derive_upstream_coords

  Arg        : none
  Example    : my ($upstart, $upend) 
                       = $self->derive_upstream_coords;
  Description: Derives upstream coordinates (for 
               compatability with older scripts).
  Returntype : arrayref
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub derive_upstream_coords {
  my $self = shift;

  return [$self->upstart, $self->upend]
}

=head2 derive_downstream_coords

  Arg        : none
  Example    : my ($downstart, $downend) 
                       = $self->derive_downstream_coords;
  Description: Derives downstream coordinates (for 
               compatability with older scripts).
  Returntype : arrayref
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub derive_downstream_coords {
  my $self = shift;

  return [$self->downstart, $self->downend]
}

=head2 _derive_coords

  Arg        : string $direction (either 'up' or 'down').
  Example    : $self->_derive_coords('up');
  Description: Determines the coordinates of either upstream 
               or downstream region.
  Returntype : none
  Exceptions : Throws if argument is not either 'up' or 'down'
  Caller     : $self->upstart, $self->upend, $self->downstart, 
               $self->downend
  Status     : Stable

=cut

sub _derive_coords {
  my ($self, $direction) = @_;

  # Check direction
  throw("Must specify either \'up\' of \'down\'-stream direction to derive coords.")
    unless (($direction eq 'up')||($direction eq 'down'));

  # Put things in easily accessible places.
  my $core_db_slice_adaptor = $self->transcript->slice->adaptor;
  my $region_length = $self->length;

  # Whatever coord system the gene is currently is, transform to the toplevel.
  my $transcript = $self->transcript->transform('toplevel');

  # Use our transformed transcript to determine the upstream region coords.
  # End should always be just before the coding start (like ATG), including 3' UTR.
  # Start is the outer limit of the region upstream (furthest from ATG).

  my $region_start;
  my $region_end;

  if ($transcript->strand == 1){
    if ($direction eq 'up'){
      $region_end = $transcript->coding_region_start - 1;
      $region_start = $region_end - $region_length;
    } elsif ($direction eq 'down'){
      $region_end = $self->_first_coding_Exon->end + 1;
      $region_start = $region_end + $region_length;
    }
  } elsif ($transcript->strand == -1) {
    if ($direction eq 'up'){
      $region_end = $transcript->coding_region_end + 1;
      $region_start = $region_end + $region_length;

    } elsif ($direction eq 'down'){
      $region_end = $self->_first_coding_Exon->start - 1;
      $region_start = $region_end - $region_length;
    }
  }

  # Trim the upstream/downstream region to remove extraneous coding sequences 
  # from other genes and/or transcripts.
    
  my ($slice_low_coord, $slice_high_coord) = sort {$a <=> $b} ($region_start, $region_end);

  my $region_slice 
      = $core_db_slice_adaptor->fetch_by_region($transcript->slice->coord_system->name,
						$transcript->slice->seq_region_name, 
						$slice_low_coord, 
						$slice_high_coord);

  if ($transcript->strand == 1) {
    if ($direction eq 'up') {
      $region_start += $self->_bases_to_trim('left_end', $region_slice);
    } elsif ($direction eq 'down') {
      $region_start -= $self->_bases_to_trim('right_end', $region_slice);
    }
  } elsif ($transcript->strand == -1) {
    if ($direction eq 'up') {
      $region_start -= $self->_bases_to_trim('right_end', $region_slice);
    } elsif ($direction eq 'down') {
      $region_start += $self->_bases_to_trim('left_end', $region_slice);
    }
  }

  # Always return start < end

  ($region_start, $region_end) = sort {$a <=> $b} ($region_start, $region_end);

  if ($direction eq 'up') {
    $self->upstart($region_start);
    $self->upend($region_end);
  } elsif ($direction eq 'down') {
    $self->downstart($region_start);
    $self->downend($region_end);
  }
}

=head2 _bases_to_trim

  Arg        : string $end_to_trim (either 'right_end' or 
               'left_end').
  Arg        : Bio::EnsEMBL::Slice
  Example    : $self->_derive_coords('right_end', $slice);
  Description: Finds exons from other genes/transcripts that
               invade our upstream/downstream slice and 
               returns the number of bases that should be 
               truncated from the appropriate end of the 
               upstream/downstream region.
  Returntype : in
  Exceptions : Throws if argument is not either 'right_end' 
               or 'left_end'
  Caller     : $self->_derive_coords
  Status     : Stable

=cut

# Method to look for coding regions that invade the upstream region.  For
# now, this method returns the number of bases to trim.  I doesn't yet
# do anything special if an exon is completely swallowed (truncates at 
# the end of the overlapping exon and discards any non-coding sequence 
# further upstream) or overlaps the 'wrong' end of the region (cases where 
# two alternate exons share one end of sequence - does this happen?).

# The input argument 'end' defines the end of the slice that should be 
# truncated. 

sub _bases_to_trim {
    my ($self, $end_to_trim, $slice) = @_;

    throw "Slice end argument must be either left_end or right_end" 
	unless ($end_to_trim eq 'right_end' || $end_to_trim eq 'left_end'); 

    my @overlap_coords;
    my $slice_length = $slice->length;
    my $right_trim = 0;
    my $left_trim  = 0;

    foreach my $exon (@{$slice->get_all_Exons}){
      next if $exon->stable_id eq $self->_first_coding_Exon->stable_id;

      my $start = $exon->start;
      my $end   = $exon->end;

      # Choose from four possible exon arrangements

      #  -----|********************|----- Slice
      #  --|=========================|--- Exon arrangement 1
      #  ----------|======|-------------- Exon arrangement 2
      #  --|=======|--------------------- Exon arrangement 3
      #  -------------------|=========|-- Exon arrangement 4


      if ($start <=  0 && $end >= $slice_length) {     # exon arrangement 1
	$right_trim = $slice_length - 1;
	$left_trim  = $slice_length - 1;
	last;

      } elsif ($start >= 0 && $end <= $slice_length) { # exon arrangement 2
	my $this_right_trim = ($slice_length - $start) + 1;

	$right_trim = $this_right_trim 
	  if $this_right_trim > $right_trim;

	$left_trim  = $end 
	  if $end > $left_trim;

      } elsif ($start <= 0 && $end < $slice_length) {  # exon arrangement 3
	$right_trim = $slice_length; # a bit draconian
	$left_trim  = $end 
	  if $end > $left_trim;

      } elsif ($start > 0 && $end >= $slice_length) {  # exon arrangement 4
	my $this_right_trim = ($slice_length - $start) + 1;

	$right_trim = $this_right_trim 
	  if $this_right_trim > $right_trim;

	$left_trim = $slice_length; # also a bit draconian
      }

    }

    return $right_trim  if $end_to_trim eq 'right_end';
    return $left_trim   if $end_to_trim eq 'left_end';
}

=head2 _first_coding_Exon

  Arg        : none
  Example    : $self->_first_coding_Exon;
  Description: Finds the first exon of our transcript that 
               contains coding bases.
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : none
  Caller     : $self->_derive_coords, $self->_bases_to_trim
  Status     : Stable

=cut

sub _first_coding_Exon {
  my $self = shift;

  unless ($self->{_first_coding_exon}){

    my $exons = $self->transcript->get_all_translateable_Exons;

    $self->{_first_coding_exon} =  $exons->[0]  
      if $self->transcript->strand == 1;
    $self->{_first_coding_exon} =  $exons->[-1] 
      if $self->transcript->strand == -1;
  }

  return $self->{_first_coding_exon}
}


return 1;
