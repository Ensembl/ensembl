package Bio::EnsEMBL::PepDnaAlignFeature;

# EnsEMBL module for storing dna-dna pairwise alignments
#
# Cared for by Michele Clamp <michele@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself
#

=head1 NAME

  Bio::EnsEMBL::DnaDnaAlignFeature - Ensembl specific dna-dna pairwise alignment feature

=head1 SYNOPSIS

  my $feat = new Bio::EnsEMBL::DnaDnaAlignFeature(-seqname => 'myseq',
						  -start   => 100,
						  -end     => 120,
						  -strand  => 1,
						  -hstart  => 200,
						  -hend    => 220,
						  -analysis    => $analysis,
						  -cigar_string => '100,200,3:110,210,11');

  # Alternatively if you have an array of ungapped features

      my $feat = new Bio::EnsEMBL::DnaDnaAlignFeature(-features => \@features);

  # Where @features is an array of Bio::EnsEMBL::FeaturePair

  # There is a method to manipulate the cigar_string into ungapped features

      my @ungapped_features = $feat->ungapped_features;

  # This converts the cigar string into an array of Bio::EnsEMBL::FeaturePair

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


  # To make things clearer this is how a blast HSP would be parsed

  #>AK014066
  #       Length = 146

  #  Minus Strand HSPs:

  #  Score = 76 (26.8 bits), Expect = 1.4, P = 0.74
  #  Identities = 20/71 (28%), Positives = 29/71 (40%), Frame = -1

  #Query:   479 GLQAPPPTPQGCRLIPPPPLGLQAPLPTLRAVGSSHHHP*GRQGSSLSSFRSSLASKASA 300
  #             G  APPP PQG R   P P G + P   L             + + ++  R  +A   +
  #Sbjct:     7 GALAPPPAPQG-RWAFPRPTG-KRPATPLHGTARQDRQVRRSEAAKVTGCRGRVAPHVAP 64

  #Query:   299 SSPHNPSPLPS 267
  #                H P+P P+
  #Sbjct:    65 PLTHTPTPTPT 75

  #The alignment goes from 267 to 479 in sequence 1 and 7 to 75 in sequence 2 and the
  #strand is -1.

  #The alignment is made up of the following ungapped pieces :

  #sequence 1 start 447 , sequence 2 start 7  , match length 33 , strand -1
  #sequence 1 start 417 , sequence 2 start 18 , match length 27 , strand -1
  #sequence 1 start 267 , sequence 2 start 27 , match length 137 , strand -1

  #These ungapped pieces are made up into the following string (called a cigar string)

  #447,7,-33:417,18,-27:267,27,-137

  #i.e. seqstart1,seqstart2,length: etc


=cut 


use Bio::EnsEMBL::BaseAlignFeature;

use vars qw(@ISA);
use strict;

@ISA = qw( Bio::EnsEMBL::BaseAlignFeature );



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
  return $self->_generic_parse_cigar( 1, 3 );
}

      
=head2 _parse_features

    Arg      : Array of Bio::EnsEMBL::FeaturePair

    Usage    : Internal method - not used.

    Function : Converts an array of FeaturePairs into a gapped feature with
               a cigar string describing the 

               See sub cigar_string for an explanation of what that is.

    Exception: If the argument passed is not an array reference

               All the features must have arisen from the same source
               i.e. a blast HSP or some other alignment.  Thus
               exceptions are thrown when the scores,percent ids,p_values
               seqnames , hseqnames and strands differ amongst the input 
               features.

               All the features must not overlap in order to provide a 
               sane gapped alignment.  An exception is thrown if they do.

               If any element of the array is not a Bio::EnsEMBL::FeaturePair

               If there are no elements in the array

               If the hit length is not equal to the query length

    Caller   : Called internally to the module by the constructor

=cut

sub _parse_features {
  my ($self,$features) = @_;

  $self->_generic_parse_features( $features, 1, 3 );
}



sub transform{
  my ($self, $slice) = @_;

  if( ! defined $slice ) {
    #Since slice arg is not defined -  we want raw contig coords
    if(( defined  $self->contig ) && 
       ( $self->contig->isa( "Bio::EnsEMBL::RawContig" )) ) {
      print STDERR "DnaPepAlignFeature::tranform, you are already apparently in rawcontig coords so why try to transform to them\n";
      #we are already in rawcontig coords, nothing needs to be done
      return $self;
    } else {
      #transform to raw_contig coords from Slice coords
      return $self->_transform_to_rawcontig();
    }
  }

  if( defined $self->contig ) {  
    if($self->contig->isa( "Bio::EnsEMBL::RawContig" ))  {
      #transform to slice coords from raw contig coords
      return $self->_transform_to_slice( $slice );
    } elsif($self->contig->isa( "Bio::EnsEMBL::Slice" )) {
      #transform to slice coords from other slice coords
      return $self->_transform_between_slices( $slice );
    } else {
      #Unknown contig type - throw an exception
      return $self->throw("Exon's 'contig' is of unknown type " 
		   . $self->contig() . " - cannot transform to Slice coords");
    }
  } else {
    #Can't convert to slice coords without a contig to work with
    return $self->throw("Exon's contig is not defined - cannot transform to " .
			"Slice coords");
  }
}

sub _transform_to_slice{
  my ($self, $slice) = @_;

  $self->throw("_transform_to_slice not yet implemented in ".$self." complain to arne\n");
  #return $self->_generic_transform_to_slice($slice, 1, 1);
}

sub _transform_to_rawcontig{
  my ($self, $rc) = @_;
  
  return $self->_generic_transform_to_rawcontig($rc, 3, 1);
}

sub _transform_between_slices{
  my ($self, $to_slice) = @_;
  
  $self->throw("_transform_between_slices not yet implemented in ".$self." complain to arne\n");
   #return $self->_generic_transform_to_slice($to_slice, 1, 1);

}

sub _create_new_feature{
  my ($self, $query_unit, $hit_unit) = @_;

  my $new_feature;
  my $f1 = new Bio::EnsEMBL::SeqFeature();
  my $f2 = new Bio::EnsEMBL::SeqFeature();
 
  $new_feature = Bio::EnsEMBL::PepDnaAlignFeature->new( -cigar_string => $self->cigar_string, 
							-feature1 => $f1, 
							-feature2 => $f2);
 

  return $new_feature;
}

1;
