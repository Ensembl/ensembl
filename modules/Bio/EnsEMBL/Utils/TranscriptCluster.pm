=head1 NAME

TranscriptCluster

=head1 SYNOPSIS


=head1 DESCRIPTION

This object holds one or more transcripts which have been clustered according to 
comparison criteria external to this class (for instance, in the 
method _compare_Transcripts of the class GeneComparison).
Each TranscriptCluster object holds the IDs of the transcripts clustered and the beginning and end coordinates
of each one (taken from the start and end coordinates of the first and last exon in the correspondig
get_all_Exons array)

It inherits from Bio::Root::RootI and Bio::RangeI. A TranscriptCluster is a range in the sense that it convers
a defined extent of genomic sequence. It is also possible to check whether two clusters overlap (in range),
is included into another cluster, etc...

=head1 CONTACT

eae@sanger.ac.uk

=cut

# Let the code begin ...

package Bio::EnsEMBL::Utils::TranscriptCluster;
use Bio::EnsEMBL::Transcript;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Gene;
use Bio::Root::RootI;
use Bio::RangeI;

@ISA = qw(Bio::Root::RootI Bio::RangeI);

=head1 METHODS

=cut

#########################################################################


=head2 new()

new() initializes the attributes:
_transcript_array
_transcriptID_array
_start
_end

=cut

sub new {
  my ($class,$whatever)=@_;

  if (ref($class)){
    $class = ref($class);
  }
  my $self = {};
  bless($self,$class);
  
  if ($whatever){
    $self->throw( "Can't pass an object to new() method. Use put_Genes() to include Bio::EnsEMBL::Gene in cluster");
  }

  return $self;
}
#########################################################################

=head1 Range-like methods

Methods start and end are typical for a range. We also implement the boolean
and geometrical methods for a range.

=head2 start()

  Title   : start
  Usage   : $start = $transcript_cluster->end();
  Function: get/set the start of the range covered by the cluster. This is re-calculated and set everytime
            a new transcript is added to the cluster
  Returns : a number
  Args    : optionaly allows the start to be set

=cut

sub start{
  my ($self,$start) = @_;
  if ($start){
    $self->throw( "$start is not an integer") unless $start =~/^[-+]?\d+$/;
    $self->{'_start'} = $start;
  }
  return $self->{'_start'};
}

############################################################

=head2 end()

  Title   : end
  Usage   : $end = $transcript_cluster->end();
  Function: get/set the end of the range covered by the cluster. This is re-calculated and set everytime
            a new transcript is added to the cluster
  Returns : the end of this range
  Args    : optionaly allows the end to be set
          : using $range->end($end
=cut

sub end{
  my ($self,$end) = @_;
  if ($end){
    $self->throw( "$end is not an integer") unless $end =~/^[-+]?\d+$/;
    $self->{'_end'} = $end;
  }
  return $self->{'_end'};
}

############################################################

=head2 length

  Title   : length
  Usage   : $length = $range->length();
  Function: get/set the length of this range
  Returns : the length of this range
  Args    : optionaly allows the length to be set
          : using $range->length($length)

=cut

sub length{
  my $self = shift @_;
  if (@_){
    $self->confess( ref($self)."->length() is read-only");
  }
  return ( $self->{'_end'} - $self->{'_start'} + 1 );
}

############################################################

=head2 strand

  Title   : strand
  Usage   : $strand = $transcript->strand();
  Function: get/set the strand of the transcripts in the cluster.
            The strand is set in put_Transcripts when the first transcript is added to the cluster, in that
            that method there is also a check for strand consistency everytime a new transcript is added
  Returns : the strandidness (-1, 0, +1)
  Args    : optionaly allows the strand to be set
        
=cut

sub strand{
  my ($self,$strand) = @_;
  if ($strand){
    $self->{'_strand'} = $strand;
  }
  return $self->{'_strand'};
}


############################################################

=head1 Boolean Methods

These methods return true or false. They throw an error if start and end are
not defined. They are implemented in Bio::RangeI.

 $cluster->overlaps($other_cluster) && print "Clusters overlap\n";

=head2 overlaps

  Title   : overlaps
  Usage   : if($cluster1->overlaps($cluster)) { do stuff }
  Function: tests if $cluster2 overlaps $cluster1 overlaps in the sense of genomic-range overlap,
            it does NOT test for exon overlap.
  Args    : arg #1 = a range to compare this one to (mandatory)
            arg #2 = strand option ('strong', 'weak', 'ignore') 
  Returns : true if the clusters overlap, false otherwise

=cut

=head2 contains

  Title   : contains
  Usage   : if($cluster1->contains($cluster2) { do stuff }
  Function: tests whether $cluster1 totally contains $cluster2 
  Args    : arg #1 = a range to compare this one to (mandatory)
	             alternatively, integer scalar to test
            arg #2 = strand option ('strong', 'weak', 'ignore') (optional)
  Returns : true if the argument is totaly contained within this range

=cut

=head2 equals

  Title   : equals
  Usage   : if($cluster1->equals($cluster2))
  Function: test whether the range covered by $cluster1 has the same start, end, length as the range 
            covered by $cluster2
  Args    : a range to test for equality
  Returns : true if they are describing the same range

=cut

=head1 Geometrical methods

These methods do things to the geometry of ranges, and return
Bio::RangeI compliant objects or triplets (start, stop, strand) from
which new ranges could be built. They are implemented  here and not in Bio::RangeI, since we
want them to return a new TranscriptCluster object.

=head2 overlap_extent

 Title   : overlap_extent
 Usage   : ($a_unique,$common,$b_unique) = $a->overlap_extent($b)
 Function: Provides actual amount of overlap between two different
           ranges. Implemented already in RangeI
 Example :
 Returns : array of values for 
           - the amount unique to a
           - the amount common to both
           - the amount unique to b
 Args    : 

=cut

#########################################################################

=head2 intersection

  Title   : intersection
  Usage   : $intersection_cluster = $cluster1->intersection($cluster2)
  Function: gives a cluster with the transcripts which fall entirely within the intersecting range of
            $cluster1 and $cluster2
  Args    : arg #1 = a cluster to compare this one to (mandatory)
            arg #2 = strand option ('strong', 'weak', 'ignore') (not implemented here)
  Returns : a TranscriptCluster object, which is empty if the intersection does not contain
            any transcript
=cut

sub intersection{
  my ($self, $cluster) = @_;

  # if either is empty, return an empty cluster
  if ( scalar( $self->get_Transcripts) == 0 ){
    $self->warn( "cluster $self is empty, returning an empty TranscriptCluster");
    my $empty_cluster = Bio::EnsEMBL::Utils::TranscriptCluster->new();
    return $empty_cluster;
  }

  if ( scalar( $cluster->get_Transcripts ) == 0 ){
    $self->warn( "cluster $cluster is empty, returning an empty TranscriptCluster");
    my $empty_cluster = Bio::EnsEMBL::Utils::TranscriptCluster->new();
    return $empty_cluster;
  }

  my @transcripts = $self->get_Transcripts;
  push( @transcripts, $cluster->get_Transcripts);

  # make an unique list of transcripts, in case they are repeated
  my %list;
  foreach my $transcript (@transcripts){
    $list{$transcript} = $transcript;
  }
  @transcripts = values( %list );

  my ($inter_start,$inter_end);
  my ($start1,$end1) = ($self->start   ,   $self->end);
  my ($start2,$end2) = ($cluster->start,$cluster->end);

  my $strand = $cluster->strand;

  # if clusters overlap, calculate the intersection
  if ( $self->overlaps( $cluster ) ){
    if ( $start2 >= $start1 && $end2 >= $end1 ){
      if ( $strand == 1){
        $inter_start = $start2;
        $inter_end   = $end1;
      }
      else{
        $inter_start = $start1;
        $inter_end   = $end2;
      }
    }
    if ( $start2 >= $start1 && $end2 < $end1){
      if ( $strand == 1){
        $inter_start = $start2;
        $inter_end   = $end2;
      }
      else{
        $inter_start = $start1;
        $inter_end   = $end1;
      }
    }
    if ( $start2 < $start1 && $end2 < $end1 ){
      if ( $strand == 1){
        $inter_start = $start1;
        $inter_end   = $end2;
      }
      else{
        $inter_start = $start2;
        $inter_end   = $end1;
      }
    }
    if ( $start2 < $start1 && $end2 >= $end1 ){
      if ( $strand == 1){
        $inter_start = $start1;
        $inter_end   = $end1;
      }
      else{
        $inter_start = $start2;
        $inter_end   = $end2;
      }
    }
  }
  else{
    $self->warn( "clusters $self and $cluster do not intersect range-wise, returning an empty TranscriptCluster");
    my $empty_cluster = Bio::EnsEMBL::Utils::TranscriptCluster->new();
    return $empty_cluster;
  }

  my $inter_cluster = Bio::EnsEMBL::Utils::TranscriptCluster->new();
  my @inter_transcripts;

  # see whether any transcript falls within this intersection
  foreach my $transcript ( @transcripts ){
    my ($start,$end) = ($self->_get_start($transcript), $self->_get_end($transcript));
    if ($strand == 1 &&  $start >= $inter_start && $end <= $inter_end ){
       $inter_cluster->put_Transcripts( $transcript );
    }
    elsif ( $strand == -1 && $start <= $inter_start && $end >= $inter_end ){
      $inter_cluster->put_Transcripts( $transcript );
    }
  }

  if ( scalar( $inter_cluster->get_Transcripts ) == 0 ){
     $self->warn( "cluster $inter_cluster is empty, returning an empty TranscriptCluster");
     return $inter_cluster;
  }
  else{
    return $inter_cluster;
  }
}

############################################################

=head2 union

  Title   : union
  Usage   : $union_cluster = $cluster1->union(@clusters);
  Function: returns the union of clusters 
  Args    : a TranscriptCluster or list of TranscriptClusters to find the union of
  Returns : the TranscriptCluster object containing all of the ranges

=cut

sub union{
my ($self,@clusters) = @_;

if ( ref($self) ){
 unshift @clusters, $self;
}

my $union_cluster = Bio::EnsEMBL::Utils::TranscriptCluster->new();
my $union_strand;

foreach my $cluster (@clusters){
 unless ($union_strand){
  $union_strand = $cluster->strand;
 }
 unless ( $cluster->strand == $union_strand){
  $self->warn("You're making the union of clusters in opposite strands");
 }
 $union_cluster->put_Transcripts($cluster->get_Transcripts);
}

return $union_cluster;


}

############################################################

=head2 put_Transcripts()

  function to include one or more transcripts in the cluster.
  Useful when creating a cluster. It takes as argument an array of transcripts, it returns nothing.

=cut

sub put_Transcripts {
  my ($self, @new_transcripts)= @_;
  
  if ( !$new_transcripts[0]->isa('Bio::EnsEMBL::Transcript') ){
    $self->throw( "Can't accept a [ $new_transcripts[0] ] instead of a Bio::EnsEMBL::Transcript");
  }
  
  my @starts = sort { $a <=> $b } map( { $self->_get_start($_) } @new_transcripts );
  my @ends   = sort { $a <=> $b } map( { $self->_get_end($_)   } @new_transcripts );

  # check strand consistency among transcripts
  foreach my $transcript (@new_transcripts){
    my @exons = $transcript->get_all_Exons;
    unless ( $self->strand ){
      $self->strand( $exons[0]->strand );
    }
    if ( $self->strand != $exons[0]->strand ){
      $self->warn( "You're trying to put $transcript in a cluster of opposite strand");
    }
  }

  # if start is not defined, set it
  unless ( $self->start ){
    if ( $self->strand == 1 ){
      $self->start( $starts[0] );
    }
    if ( $self->strand == -1 ){
      $self->start( $starts[$#starts] );
    }
  }

  # if end is not defined, set it
  unless ( $self->end ){
    if ( $self->strand == 1 ){
      $self->end( $ends[$#ends]);
    }
    if ( $self->strand == -1 ){
      $self->end( $ends[0] );
    }
  }
  
  # extend start and end if necessary as we include more transcripts
  if ( $self->strand == 1 && $starts[0] < $self->start ){
    $self->start( $starts[0] );
  }
  if ( $self->strand == -1 && $starts[$#starts] > $self->start ){
    $self->start( $starts[$#starts] );
  }
  if ( $self->strand == 1 && $ends[$#ends] > $self->end ){
    $self->end( $ends[$#ends] );
  }
  if ( $self->strand == -1 && $ends[0] < $self->end ){
    $self->end( $ends[0] );
  }

  push ( @{ $self->{'_transcript_array'} }, @new_transcripts );

}

#########################################################################

=head2 get_Transcripts()

  it returns the array of transcripts in the GeneCluster object

=cut

sub get_Transcripts {
  my $self = shift @_;
  my @transcripts = @{ $self->{'_transcript_array'} };
  return @transcripts;
}

#########################################################################

=head2 to_String()

  it returns a string containing the information about the transcripts in the TranscriptCluster object

=cut

sub to_String {
  my $self = shift @_;
  my $data='';
  foreach my $tran ( @{ $self->{'_transcript_array'} } ){
    my @exons = $tran->get_all_Exons;
     
    $data .= sprintf "Id: %-16s"             , $tran->stable_id;
    $data .= sprintf "Contig: %-21s"         , $exons[0]->contig->id;
    $data .= sprintf "Exons: %-3d"           , scalar(@exons);
    $data .= sprintf "Start: %-9d"           , $self->_get_start($tran);
    $data .= sprintf "End: %-9d"             , $self->_get_end  ($tran);
    $data .= sprintf "Strand: %-3d"          , $exons[0]->strand;
    $data .= sprintf "Exon-density: %3.2f\n", $self->exon_Density($tran);
  }
  return $data;
}

#########################################################################

sub exon_Density{
  my ($self, $transcript) = @_;  
  my $density;
  my $exon_span;
  my @exons = $transcript->get_all_Exons;
  @exons = sort { $a->start <=> $b->start } @exons;
  my $transcript_length = $exons[$#exons]->end - $exons[0]->start;
  foreach my $exon ( @exons ){
    $exon_span += $exon->length;
  }
  $density = $exon_span/$transcript_length;
  return $density;
}

#########################################################################

=head2 _get_start()

 function to get the start position of a transcript - it reads the Bio::EnsEMBL::Transcript 
 object and it returns the start position of the first exon

=cut

sub _get_start {
  my ($self,$transcript) = @_;
  my @exons = $transcript->get_all_Exons;
  my $st;
  
  if ($exons[0]->strand == 1) {
    @exons = sort {$a->start <=> $b->start} @exons;
    $st = $exons[0]->start;
  } else {
    @exons = sort {$b->start <=> $a->start} @exons;
    $st = $exons[0]->end;
  }

  return $st;
}

#########################################################################

=head2 _get_end()

 function to get the end position of a transcript - it reads the Bio::EnsEMBL::Transcript 
 object and it returns the end position of the last exon

=cut

sub _get_end {
  my ($self,$transcript) = @_;
  my @exons = $transcript->get_all_Exons;
  my $end;
  
  if ($exons[0]->strand == 1) {
    @exons = sort {$a->start <=> $b->start} @exons;
    $end = $exons[$#exons]->end;
  } else {
    @exons = sort {$b->start <=> $a->start} @exons;
    $end = $exons[$#exons]->start;
  }
  return $end;
}
#########################################################################

1;
