
# EnsEMBL module for storing dna-protein pairwise alignments
#
# Cared for by Michele Clamp <michele@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself
#

=head1 NAME

  Bio::EnsEMBL::DnaPepAlignFeature - Ensembl specific dna-protein pairwise alignment feature

=head1 SYNOPSIS

  my $feat = new Bio::EnsEMBL::DnaPepAlignFeature(-seqname => 'myseq',
						  -start   => 100,
						  -end     => 120,
						  -strand  => 1,
						  -hstart  => 200,
						  -hend    => 220,
						  -analysis    => $analysis,
						  -cigar_string => '100,200,3:109,203,12');

  # Alternatively if you have an array of ungapped features

      my $feat = new Bio::EnsEMBL::DnaPepAlignFeature(-features => \@features);

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
  # 100-102 in the query and 200-200 in the hit and
  # 109-120 in the query and 203-206 in the hit.
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

package Bio::EnsEMBL::BaseAlignFeature;

use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::FeaturePair);

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
    
    my ($cigar_string,$features) = $self->_rearrange([qw(CIGAR_STRING FEATURES)],
						     @args);

    if (defined($cigar_string) && defined($features)) {
      $self->throw("Can't input cigar_string and an array of features.");
    } elsif (defined($features)) {
      $self->_parse_features($features);
    } elsif (!defined($cigar_string)) {
      $self->throw("Must have a cigar string defining the alignment");
    } else {
      $self->throw("Wrong arguments in constructor");
    }

    $self->cigar_string($cigar_string);

    return $self;
}

=head2 cigar_string

    Arg      : Nothing (for just returning an existing cigar_string) or
               A cigar string which looks something like

               '100,200,3:109,203,12'

               The cigar_string contains the ungapped pieces that make up the gapped alignment.

               Its format is qstart,qend,length*strand.

               So in the above example the gapped alignment contains 2 ungapped pieces from

                  100-102 in the query and 200-200 in the hit and
                  109-120 in the query and 203-206 in the hit.
  
               The length parts of the cigar string are positive as the strand is +ve.

    Usage    : $pep_feat->cigar_string('100,200,3:109,203,12');
               my $cigar_string = $pep_feat->cigar_string

    Function : This module is for storing gapped dna-protein alignments.  To save 
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



=head2 _generic_parse_cigar

  Arg 1     : int query_unit, 
                  either 1 (peptide) or 3 (dna)

  Arg 2     : int hit_unit, 
                  either 1 (peptide) or 3 (dna)
  
             a DNA,DNA alignment should use 1,1 for the query_unit and hit_unit        
  Function  : Converts the cigar_string contained by the module into 
              an array of ungapped Bio::EnsEMBL::FeaturePair.

              See sub cigar_string for an explanation of what that is.

  Returntype: list of Bio::EnsEMBL::FeaturePairs
  Exceptions: Mainly internal, eg, no cigar line, wrong arguments or
              miswritten query_unit or hit_unit

  Caller    : _parse_cigar in derived class

=cut


sub _generic_parse_cigar {
  my ( $self, $query_unit, $hit_unit ) = @_;

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
    $feature1->seqname($self->seqname);
    $feature1->phase($self->phase);
    $feature1->p_value($self->p_value);
    $feature1->percent_id($self->percent_id);

    my $feature2 = new Bio::EnsEMBL::SeqFeature();

    my $mapped_length;

    # explicit if statements to avoid rounding problems
    # and make sure we have sane coordinate systems
    if( $query_unit == 1 && $hit_unit == 3 ) {
      $mapped_length = $length*3;
    } elsif( $query_unit == 3 && $hit_unit == 1 ) {
      $mapped_length = $length / 3;
    } elsif ( $query_unit == 1 && $hit_unit == 1 ) {
      $mapped_length = $length;
    } else {
      $self->throw("Internal error $query_unit $hit_unit, currently only allowing 1 or 3 ");
    }

    if( int($mapped_length) != $mapped_length ) {
      $self->throw("Internal error with mismapped length of hit, query $query_unit, hit $hit_unit, length $length");
    }
    

    if ($strand == 1) {
      $feature2->start($start2);
      $feature2->end  ($start2 + $strand*($mapped_length - 1));
    } else {
      $feature2->end  ($start2);
      $feature2->start($start2 + $strand*($mapped_length - 1));
    }      
    $feature2->strand($strand);
    $feature2->score($self->score);
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

      
=head2 _generic_parse_features

    Arg 1   : listref Bio::EnsEMBL::FeaturePair $features

    Arg 2   : int query_unit, 
                  either 1 (peptide) or 3 (dna)

    Arg 3   : int hit_unit, 
                  either 1 (peptide) or 3 (dna)
  
             a DNA,DNA alignment should use 1,1 for the query_unit and hit_unit    

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

               If the hit length is not exactly 3 times the query length

    Caller   : Called internally to the module by the constructor

=cut

sub _generic_parse_features {
  my ($self,$features, $query_unit, $hit_unit ) = @_;

  if (ref($features) ne "ARRAY") {
    $self->throw("features must be an array reference not a [" . ref($features) . "]");
  }

  my @f = sort {$a->start <=> $b->start} @$features;

  if (scalar(@f) == 0) {
    $self->throw("No features in the array to parse");
  }
  my $hstrand     = $f[0]->hstrand;
  my $strand     = $f[0]->strand;
  my $name        = $f[0]->seqname;
  my $hname       = $f[0]->hseqname;
  my $score       = $f[0]->score;
  my $percent     = $f[0]->percent_id;
  my $pvalue      = $f[0]->p_value;
  my $analysis    = $f[0]->analysis;
  my $phase       = $f[0]->phase;

  print STDERR "Analysis is $analysis\n";
  print STDERR "Strand is $hstrand\n";

  my $prev1;
  my $prev2;

  my $string;
 
  my $f1start = $f[0]->start;
  my $f1end   = $f[$#f]->end;

  my $f2start;
  my $f2end;

  if ($hstrand == 1) {
    $f2start = $f[0]->hstart;
    $f2end   = $f[$#f]->hend;
  } else {
    $f2end   = $f[0]->hend;
    $f2start = $f[$#f]->hstart;
  }

  foreach my $f (@f) {
    if (!$f->isa("Bio::EnsEMBL::FeaturePair")) {
      $self->throw("Array element [$f] is not a Bio::EnsEMBL::FeaturePair");
    }

    print STDERR "Processing " . $f->gffstring . "\n";
    if( defined $hstrand ) {
      if ($f->hstrand != $hstrand) {
	$self->throw("Inconsistent hstrands in feature array");
      }
    }
    if( defined $strand ) {
      if ($f->strand != $strand) {
	$self->throw("Inconsistent strands in feature array");
      }
    }

    if ($name ne $f->seqname) {
      $self->throw("Inconsistent names in feature array [$name - " . $f->seqname . "]");
    }
    if ($hname ne $f->hseqname) {
      $self->throw("Inconsistent names in feature array [$hname - " . $f->hseqname . "]");
    }
    if ($score ne $f->score) {
      $self->throw("Inconsisent scores in feature array [$score - " . $f->score . "]");
    }
    if ($percent ne $f->percent_id) {
      $self->throw("Inconsistent pids in feature array [$percent - " . $f->percent_id . "]");
    }

    my $start1 = $f->start;
    my $start2;

    if ($hstrand == 1) {
       $start2 = $f->hstart;
    } else {
       $start2 = $f->hend;
    }

    if (defined($prev1)) {
      if ($start1 < $prev1) {
        $self->throw("Inconsistent coordinates in feature array $start1 < $prev1");
      }
    }

    if (defined($prev2)) {
      if ($hstrand == 1) {
        if ($f->hstart < $prev2) {
          $self->throw("Inconsistent coordinates in feature array " . $f->hstart . " < $prev2");
        }
      } else {
        if ($f->hstart > $prev2) {
          $self->throw("Inconsistent coordinates in feature array " . $f->hstart . " > $prev2");
        }
      }

    }
    my $length = ($f->end - $f->start + 1)*$hstrand;
    my $hlength = ($f->hend - $f->hstart + 1)*$hstrand;

    # using multiplication to avoid rounding errors, hence the
    # switch from query to hit for the ratios
    if( $length * $hit_unit != $hlength * $query_unit ) {
      $self->throw( "Feature lengths not comparable Lengths:".$length." ".$hlength." Ratios:".$query_unit." ". $hit_unit );
    }

    $string = $string . $start1 . "," . $start2 . "," . $length . ":";

    $prev1 = $f->end;

    if ($hstrand == 1) {
      $prev2 = $f->hend;
    } else {
      $prev2 = $f->hstart;
    }

  }

  $string =~ s/:$//;

  my $feature1 = new Bio::EnsEMBL::SeqFeature();
 
  $feature1->start($f1start);
  $feature1->end  ($f1end);
  $feature1->strand($strand);
  $feature1->score($score);
  $feature1->percent_id($percent);
  $feature1->p_value($pvalue);
  $feature1->seqname($name);
  $feature1->percent_id($percent);
  $feature1->p_value($pvalue);
  $feature1->phase($phase);
  $feature1->analysis($analysis);

  my $feature2 = new Bio::EnsEMBL::SeqFeature();
 
  $feature2->start($f2start);
  $feature2->end  ($f2end);
  $feature2->strand($hstrand);
  $feature2->score($score);
  $feature2->seqname($hname);
  $feature2->percent_id($percent);
  $feature2->p_value($pvalue);
  $feature2->phase($phase);
  $feature2->analysis($analysis);

  $self->feature1($feature1);
  $self->feature2($feature2);
  $self->cigar_string($string);

}

=head2 ungapped_features

    Arg      : None.  

    Usage    : my @feature_pairs = $pep_feat->ungapped_features;

    Function : Converts the cigar_string contained by the module into 
               an array of ungapped Bio::EnsEMBL::FeaturePair.

               All the work is done in the subroutine _parse_cigar_string

    Exception: If no cigar string is defined

    Caller   : No specific caller.

=cut

sub ungapped_features {
  my ($self) = @_;

  if (defined($self->cigar_string)) {
    my @features = $self->_parse_cigar;
    return @features;
  } else {
    $self->throw("No cigar_string defined.  Can't return ungapped features");
  }

}


1;
