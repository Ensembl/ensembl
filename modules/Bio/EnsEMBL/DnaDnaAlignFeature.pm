package Bio::EnsEMBL::DnaDnaAlignFeature;

# EnsEMBL module for storing dna-dna pairwise alignments
#
# Cared for by Michele Clamp <michele@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself
#

=head1 NAME

  Bio::EnsEMBL::DnaDnaAlignFeature - Ensembl specific dna-dna pairwise alignment feature

= head1 SYNOPSIS

  my $feat = new Bio::EnsEMBL::DnaDnaAlignFeature(-seqname => 'myseq',
						  -start   => 100,
						  -end     => 120,
						  -strand  => 1,
						  -hstart  => 200,
						  -hend    => 220,
						  -source_tag => 'blastn',
						  -primary_tag => 'similarity',
						  -analysis    => $analysis,
						  -cigar_string => '100,200,3:110,210,11');

  # $analysis is a Bio::EnsEMBL::Analysis object
  
  # Bio::EnsEMBL::SeqFeature methods can be used
  # Bio::EnsEMBL::FeaturePair methods can be used

  # The cigar_string contains the ungapped pieces that make up the gapped alignment
  #
  # It's format is qstart,qend,length*strand.
  #
  # So in the above example the gapped alignment contains 2 ungapped pieces from
  #
  # 100-102 in the query and 200-202 in the hit and
  # 110-120 in the query and 210-220 in the hit.
  #
  # The length parts of the cigar string are positive as the strand is +ve.

  # There is a method to manipulate the cigar_string into ungapped features

  my @ungapped_features = $feat->ungapped_features;

  # This converts the cigar string into an array of Bio::EnsEMBL::FeaturePair

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::FeaturePair);

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);

    my ($cigar_string) = $self->_rearrange([qw(CIGAR_STRING)],
					   @args);

    if (!defined($cigar_string)) {
      $self->throw("A DnaDnaAlignFeature must have a cigar string defining the alignment");
    }

    $self->cigar_string($cigar_string);

    return $self;
}

=head2 cigar_string

    Arg      : Nothing (for just returning an existing cigar_string) or
               A cigar string which looks something like

               '100,200,3:110,210,11'

               The cigar_string contains the ungapped pieces that make up the gapped alignment.

               Its format is qstart,qend,length*strand.

               So in the above example the gapped alignment contains 2 ungapped pieces from

                  100-102 in the query and 200-202 in the hit and
                  110-120 in the query and 210-220 in the hit.
  
               The length parts of the cigar string are positive as the strand is +ve.

    Usage    : $dna_feat->cigar_string('100,200,3:110,210,11');
               my $cigar_string = $dna_feat->cigar_string

    Function : This module is for storing gapped dna-dna alignments.  To save 
               space and memory all the structure of the gapped alignment
               is stored in a string and the ->start and ->end values
               are the minimum and maximum coords of all the gapped pieces.
               This allows the gapped alignment to also function as
               a feature pair.  It is mostly used for storing Blast HSPs
               (see Bio::EnsEMBL::Pipeline::Runnable::Blast)

               If a cigar_string is defined it is returned otherwise an exception is thrown

    Exception: If no cigar_string is defined an exception is thrown

    Caller   : Bio::EnsEMBL::Pipeline::Runnable::Blast

=cut



sub cigar_string {
  my ($self,$arg) = @_;

  if (defined($arg)) {

    # Do some checks to see whether this is valid
    my $tmp = $arg;
    
    $self->{_cigar_string} = $arg;
  }

  if (!defined($self->{_cigar_string})) {
    $self->throw("No cigar string defined - can't return one");
  }

  return $self->{_cigar_string};
}


=head2 _parse_cigar

    Arg      : None.  

    Usage    : Internal method - not used.

    Function : Converts the cigar_string contained by the module into 
               an array of ungapped Bio::EnsEMBL::FeaturePair.

               See sub cigar_string for an explanation of what that is.

    Exception: If no cigar_string is returned from $self->cigar_string but
               this should never happen as $self->cigar_string should throw 
               it first.

               If the cigar string is the wrong format.

               If the length of an ungapped piece is 0

    Caller   : Called internally to the module by ungapped_features

=cut


sub _parse_cigar {
  my ($self) = @_;

  my $string = $self->cigar_string;

  if (!defined($string)) {
    $self->throw("No cigar string defined in object.  This should be caught by the cigar_string method and never happen");
  }

  my @pieces = split(/:/,$string);

  my @features;

  foreach my $piece (@pieces) {

    my $start1;
    my $start2;
    my $length;
   
    if ($piece =~ /(\S+)\,(\S+)\,(\S+)/) {
      $start1 = $1; 
      $start2 = $2;
      $length = $3;
    } else {
      $self->throw("Can't parse cigar string element [$piece].  Invalid format.  Should be start,end,length*strand");
    }

    my $strand = 1;

    if ($length == 0) {
      $self->throw("Length of piece is 0 - impossible");
    }
    if ($length < 0) {
      $strand = -1;
    }

    $length = abs($length);

    my $feature1 = new Bio::EnsEMBL::SeqFeature();

    $feature1->start($start1);
    $feature1->end  ($start1 + $length - 1);
    $feature1->strand($strand);
    $feature1->score($self->score);
    $feature1->source_tag($self->source_tag);
    $feature1->primary_tag($self->primary_tag);
    $feature1->seqname($self->seqname);
    $feature1->phase($self->phase);
    $feature1->p_value($self->p_value);
    $feature1->percent_id($self->percent_id);

    my $feature2 = new Bio::EnsEMBL::SeqFeature();

    $feature2->start($start2);
    $feature2->end  ($start2 + $strand*($length - 1));
    $feature2->strand($strand);
    $feature2->score($self->score);
    $feature2->source_tag($self->source_tag);
    $feature2->primary_tag($self->primary_tag);
    $feature2->seqname($self->hseqname);
    $feature2->phase($self->phase);
    $feature2->p_value($self->p_value);
    $feature2->percent_id($self->percent_id);
    
    my $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $feature1,
					   -feature2 => $feature2);

 
    $fp->analysis($self->analysis);

    push(@features,$fp);
  }
  return @features;
}

=head2 ungapped_features

    Arg      : None.  

    Usage    : my @feature_pairs = $dna_feat->ungapped_features;

    Function : Converts the cigar_string contained by the module into 
               an array of ungapped Bio::EnsEMBL::FeaturePair.

               All the work is done in the subroutine _parse_cigar_string

    Exception: If no cigar string is defined

    Caller   : No specific caller.

=cut

sub ungapped_features {
  my ($self) = @_;

  if (defined($self->cigar)) {
    my @features = $self->_parse_cigar_string;
    return @features;
  } else {
    $self->throw("No cigar_string defined.  Can't return ungapped features");
  }

}


1;
