=head1 NAME

TranscriptCluster

=head1 SYNOPSIS


=head1 DESCRIPTION

This object holds one or more transcripts which have been clustered according to 
comparison criteria external to this class (for instance, in the 
method _compare_Transcripts of the class GeneComparison).
Each TranscriptCluster object holds the IDs of the transcripts clustered and the beginning and end coordinates
of each one (taken from the start and end coordinates of the first and last exon in the correspondig
each_Exon array)

=head1 CONTACT

eae@sanger.ac.uk

=cut

# Let the code begin ...

package TranscriptCluster;
use Bio::EnsEMBL::Transcript;
use strict;

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
  my ($class,@args)=@_;

  if (ref($class)){
    $class = ref($class);
  }
  my $self = {};
  bless($self,$class);
  $self->{'_transcript_array'}= \@args;
  $self->{'_transcriptID_array'}=();  # array that holds the IDs of the transcripts
  $self->{'_start'}={};               # hash that holds the start position of each transcript
  $self->{'_end'}={};                 # hash that holds the end position of each transcript
  
  foreach my $transcript (@args){
    my ( $start , $end)=( _get_start($transcript) , _get_end($transcript) );
    push ( @ {$self->{'_transcriptID_array'} }, $transcript->id );
    $self->{'_start'}->{$transcript->id}=$start;
    $self-> {'_end'} ->{$transcript->id}=$end;
  }
  return $self;
}

#########################################################################

=head2 put_Transcripts()

  function to include one or more transcripts in the cluster.
  Useful when creating a cluster. It takes as argument an array of transcripts, it returns nothing.

=cut

sub put_Transcripts {
  my ($self, @args)= @_;
  my @new_transcripts = @args;
  push ( @{ $self->{'_transcript_array'} }, @new_transcripts );
  foreach my $new_transcript (@new_transcripts){
    push ( @{ $self->{'_transcriptID_array'} }, $new_transcript->id );
    my ( $start , $end)=( _get_start($new_transcript) , _get_end($new_transcript) );
    $self->{'_start'}->{$new_transcript->id}=$start;
    $self-> {'_end'} ->{$new_transcript->id}=$end;
  }
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

=head2 string()

  it returns a string containing the information from a TranscriptCluster object,
  i.e. the transcriptID, the start position and the end position.

=cut

sub string {
  my $self = shift @_;
  my $data='';
  foreach my $transcript ( @{ $self->{'_transcript_array'} } ){
    my $id = $transcript->id;
    while (length($id)<16){
      $id .=' ';
    }
    $data .= $id."  "._get_start($transcript)."\t"._get_end($transcript)."\n";
  }
  return $data;
}

#########################################################################

=head2 _get_start()

 function to get the start position of a transcript - it reads the Bio::EnsEMBL::Transcript 
 object and it returns the start position of the first exon

=cut

sub _get_start {
  my $transcript = shift @_;
  my @exons = $transcript->each_Exon;
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
  my $transcript = shift @_;
  my @exons = $transcript->each_Exon;
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
