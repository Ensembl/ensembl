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

@ISA = qw(Bio::Root::RootI);

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

=head2 put_Transcripts()

  function to include one or more transcripts in the cluster.
  Useful when creating a cluster. It takes as argument an array of transcripts, it returns nothing.

=cut

sub put_Transcripts {
  my ($self, @new_transcripts)= @_;
  
  if ( !$new_transcripts[0]->isa('Bio::EnsEMBL::Transcript') ){
    $self->throw( "Can't accept a [ $new_transcripts[0] ] instead of a Bio::EnsEMBL::Transcript");
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
     
    $data .= sprintf "Id: %-16s"      , $tran->stable_id;
    $data .= sprintf "Contig: %-20s"  , $exons[0]->contig->id;
    $data .= sprintf "Exons: %-3d"    , scalar(@exons);
    $data .= sprintf "Start: %-9d"    , _get_start($tran);
    $data .= sprintf "End: %-9d"      , _get_end  ($tran);
    $data .= sprintf "Strand: %-2d\n" , $exons[0]->strand;
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
  my $transcript = shift @_;
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
