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

Bio::EnsEMBL::Map::DitagFeature

=head1 SYNOPSIS

my $feature = Bio::EnsEMBL::Map::DitagFeature->new(
  -slice         => $slice,
  -start         => $qstart,
  -end           => $qend,
  -strand        => $qstrand,
  -hit_start     => $tstart,
  -hit_end       => $tend,
  -hit_strand    => $tstrand,
  -ditag_id      => $ditag_id,
  -ditag_side    => $ditag_side,
  -ditag_pair_id => $ditag_pair_id,
  -cigar_line    => $cigar_line,
  -analysis      => $analysis,
);

=head1 DESCRIPTION

Represents a mapped ditag object in the EnsEMBL database.  These are
the original tags separated into start ("L") and end ("R") parts if
applicable, successfully aligned to the genome. Two DitagFeatures
usually relate to one parent Ditag.  Alternatively there are CAGE tags
e.g. which only have a 5\'tag ("F").

=head1 METHODS

=cut

package Bio::EnsEMBL::Map::DitagFeature;

use strict;
use vars qw(@ISA);


use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate);
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );

@ISA = qw(Bio::EnsEMBL::Feature);

=head2 new

  Arg [1]    : (optional) int dbID
  Arg [2]    : (optional) Bio::EnsEMBL::DitagFeatureAdaptor $adaptor
  Arg [3]    : int start
  Arg [4]    : int end
  Arg [5]    : int strand
  Arg [6]    : Bio::EnsEMBL::Slice $slice
  Arg [7]    : (optional) Bio::EnsEMBL::Analysis
  Arg [8]    : int hit_start
  Arg [9]    : int hit_end
  Arg [10]   : int hit_strand
  Arg [11]   : int ditag_id
  Arg [12]   : string ditag_side
  Arg [13]   : (optional) sring cigar_line
  Arg [14]   : (optional) int ditag_pair_id
  Arg [15]   : (optional) int tag_count, only used for imported mappings where
               identical positions where collapsed into into one feature.
               Default: 1
  Arg [16]   : (optional) ditag object

  Example    : $ditag = Bio::EnsEMBL::Map::DitagFeature->new
                            (-dbID => 123, -adaptor => $adaptor, ...);
  Description: Creates a new DitagFeature
  Returntype : Bio::EnsEMBL::Map::DitagFeature
  Caller     : general
  Status     : At Risk

=cut

sub new {
  my ($caller, @args) = @_;
  my ( $dbID, $adaptor, $start, $end, $strand, $slice, $analysis, $hit_start, $hit_end, 
       $hit_strand, $ditag_id, $ditag_side, $cigar_line, $ditag_pair_id, $tag_count, $ditag  ) =
	 rearrange( [ 'dbid', 'adaptor' ,'start', 'end', 'strand', 'slice', 'analysis', 'hit_start',
	'hit_end', 'hit_strand', 'ditag_id', 'ditag_side', 'cigar_line', 'ditag_pair_id' ,'tag_count', 'ditag'],
       @args );
  my $class = ref($caller) || $caller;

  if($analysis) {
    if(!ref($analysis) || !$analysis->isa('Bio::EnsEMBL::Analysis')) {
      throw('-ANALYSIS argument must be a Bio::EnsEMBL::Analysis not '.
            $analysis);
    }
  }
  if(defined($strand)) {
    if(!($strand =~ /^-?\d$/) or !($strand == 1) && !($strand == -1) && !($strand == 0)) {
      throw('-STRAND argument must be 1, -1, or 0');
    }
  }
  if(defined($hit_strand)) {
    if(!($hit_strand == 1) && !($hit_strand == -1) && !($hit_strand == 0)) {
      throw('-HIT_STRAND argument must be 1, -1, or 0 not '.$hit_strand);
    }
  }
  if(defined($start) && defined($end)) {
    if($end+1 < $start) {
      throw('Start must be less than or equal to end+1.');
    }
  }
  else{
    throw('Need start and end location.');
  }
  if(!(defined($hit_start) && defined($hit_end))) {
    throw('Need hit start and hit end location.');
  }
  if(!defined($tag_count) or (!$tag_count =~ /^[\d]+$/)){
    $tag_count = 1;
  }

  my $self = bless( {'dbID'          => $dbID,
                     'analysis'      => $analysis,
                     'slice'         => $slice,
                     'start'         => $start,
                     'end'           => $end,
		     'strand'        => $strand,
		     'hit_start'     => $hit_start,
		     'hit_end'       => $hit_end,
		     'hit_strand'    => $hit_strand,
                     'ditag_id'      => $ditag_id,
		     'ditag_pair_id' => $ditag_pair_id,
                     'ditag_side'    => $ditag_side,
                     'cigar_line'    => $cigar_line,
		     'tag_count'     => $tag_count,
		     'ditag'         => $ditag,
                    }, $class);

  $self->adaptor($adaptor);
  return $self;
}


=head2 ditag

  Arg [1]    : (optional) ditag object
  Description: Get/Set the ditag object of this DitagFeature
  Returntype : Bio::EnsEMBL::Map::Ditag
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub ditag {
  my $self = shift;

  if(@_) {
    $self->{'ditag'} = shift;
  } elsif(!$self->{'ditag'}) {
    if($self->{'ditag_id'}) {
      #lazy load the ditag
      my $ditag_adaptor = $self->analysis->adaptor->db->get_DitagAdaptor;
      $self->{'ditag'}  = $ditag_adaptor->fetch_by_dbID($self->ditag_id);
    }
    else{
      throw "Could not get Ditag for DitagFeature ".$self->dbID;
    }
  }

  return($self->{'ditag'});
}


=head2 get_ditag_location

  Arg [1]    : none
  Description: Get the start and end location (and strand ) of the start-end pair
               this DitagFeature belongs to.
               If it is not a paired ditag, these will be identical
               to DitagFeature->start() & DitagFeature->end().
               Please note that the returned start/end are min/max locations.
  Returntype : int (start, end, strand)
  Exceptions : throws if the 2 features of a pair are found on different strands
               or if the second one cannot be found.
  Caller     : general
  Status     : At Risk

=cut

sub get_ditag_location {
  my $self = shift;

  my ($start, $end, $strand);
  if($self->ditag_side eq "F"){
    $start = $self->start;
    $end   = $self->end;
  }
  else{
    my ($ditag_a, $ditag_b, $more);
    eval{
     ($ditag_a, $ditag_b, $more) = @{$self->adaptor->fetch_all_by_ditagID($self->ditag_id,
									  $self->ditag_pair_id,
									  $self->analysis->dbID)};
    };
    if($@ or !defined($ditag_a) or !defined($ditag_b)){
      throw("Cannot find 2nd tag of pair (".$self->dbID.", ".$self->ditag_id.", ".
	    $self->ditag_pair_id.", ".$self->analysis->dbID.")\n".$@);
    }
    else{
#      if(defined $more){
#	throw("More than two DitagFeatures were returned for ".$self->dbID.", ".$self->ditag_id
#	      .", ".$self->ditag_pair_id);
#      }

      ($ditag_a->start < $ditag_b->start) ? ($start = $ditag_a->start) : ($start = $ditag_b->start);
      ($ditag_a->end   > $ditag_b->end)   ? ($end   = $ditag_a->end)   : ($end   = $ditag_b->end);
      if($ditag_a->strand != $ditag_b->strand){
	throw('the strands of the two ditagFeatures are different! '.$ditag_a->strand.'/'.$ditag_b->strand);
      }
    }
  }

  return($start, $end, $self->strand);
}


=head2 ditag_id

  Arg [1]    : (optional) value
  Description: Getter/Setter for the ditag_id
               of this DitagFeature
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub ditag_id {
  my $self = shift;

  if(@_) {
    $self->{'ditag_id'} = shift;
  }

  return $self->{'ditag_id'};
}

=head2 slice

  Arg [1]    : (optional) value
  Description: Getter/Setter for the slice
               of this DitagFeature
  Returntype : slice object
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub slice {
  my $self = shift;

  if(@_) {
    $self->{'slice'} = shift;
  }

  return $self->{'slice'};
}

=head2 ditag_pair_id

  Arg [1]    : (optional) value
  Description: Getter/Setter for the ditag_pair_id
               of this DitagFeature
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut
sub ditag_pair_id {
  my $self = shift;

  if(@_) {
    $self->{'ditag_pair_id'} = shift;
  }

  return $self->{'ditag_pair_id'};
}

=head2 ditag_side

  Arg [1]    : (optional) value
  Description: Getter/Setter for the ditag_side
               of this DitagFeature
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub ditag_side {
  my $self = shift;

  if(@_) {
    $self->{'ditag_side'} = shift;
  }

  return $self->{'ditag_side'};
}

=head2 hit_start

  Arg [1]    : (optional) value
  Description: Getter/Setter for the hit_start
               of this DitagFeature
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub hit_start {
  my $self = shift;

  if(@_) {
    $self->{'hit_start'} = shift;
  }

  return $self->{'hit_start'};
}

=head2 hit_end

  Arg [1]    : (optional) value
  Description: Getter/Setter for the hit_end
               of this DitagFeature
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub hit_end {
  my $self = shift;

  if(@_) {
    $self->{'hit_end'} = shift;
  }

  return $self->{'hit_end'};
}

=head2 hit_strand

  Arg [1]    : (optional) value
  Description: Getter/Setter for the hit_strand
               of this DitagFeature
  Returntype : 1/-1/0
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub hit_strand {
  my $self = shift;

  if(@_) {
    $self->{'hit_strand'} = shift;
  }

  return $self->{'hit_strand'};
}

=head2 cigar_line

  Arg [1]    : (optional) value
  Description: Getter/Setter for the cigar_line
               of this DitagFeature
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub cigar_line {
  my $self = shift;

  if(@_) {
    $self->{'cigar_line'} = shift;
  }

  return $self->{'cigar_line'};
}

=head2 start

  Arg [1]    : (optional) value
  Description: Getter/Setter for the start
               of this DitagFeature
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub start {
  my $self = shift;

  if(@_) {
    $self->{'start'} = shift;
  }

  return $self->{'start'};
}

=head2 end

  Arg [1]    : (optional) value
  Description: Getter/Setter for the end
               of this DitagFeature
  Returntype : int or string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub end {
  my $self = shift;

  if(@_) {
    $self->{'end'} = shift;
  }

  return $self->{'end'};
}

=head2 strand

  Arg [1]    : (optional) value
  Description: Getter/Setter for the strand
               of this DitagFeature
  Returntype : 1/-1/0
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub strand {
  my $self = shift;

  if(@_) {
    $self->{'strand'} = shift;
  }

  return $self->{'strand'};
}

=head2 dbID

  Arg [1]    : (optional) value
  Description: Getter/Setter for the dbID
               of this DitagFeature
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub dbID {
  my $self = shift;

  if(@_) {
    $self->{'dbID'} = shift;
  }

  return $self->{'dbID'};
}

=head2 sequence

  Arg [1]    : (optional) value
  Description: Getter/Setter for the sequence
               of this DitagFeature
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub sequence {
  my $self = shift;

  $self->{'sequence'} = $self->adaptor->sequence($self->dbID());

  return $self->{'sequence'};
}


=head2 tag_count

  Arg [1]    : (optional) value
  Description: Getter/Setter for the tag_count
               of this DitagFeature
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub tag_count {
  my $self = shift;

  if(@_) {
    $self->{'tag_count'} = shift;
  }

  return $self->{'tag_count'};
}

1;
