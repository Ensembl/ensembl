=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::BaseAlignFeature - Baseclass providing a common abstract
implmentation for alignment features

=head1 SYNOPSIS

  my $feat = new Bio::EnsEMBL::DnaPepAlignFeature(
    -slice        => $slice,
    -start        => 100,
    -end          => 120,
    -strand       => 1,
    -hseqname     => 'SP:RF1231',
    -hstart       => 200,
    -hend         => 220,
    -analysis     => $analysis,
    -cigar_string => '10M3D5M2I'
  );

  Alternatively if you have an array of ungapped features

  my $feat =
    new Bio::EnsEMBL::DnaPepAlignFeature( -features => \@features );

  Where @features is an array of Bio::EnsEMBL::FeaturePair

  There is a method to manipulate the cigar_string into ungapped features

  my @ungapped_features = $feat->ungapped_features();

  This converts the cigar string into an array of Bio::EnsEMBL::FeaturePair

  $analysis is a Bio::EnsEMBL::Analysis object

  Bio::EnsEMBL::Feature methods can be used
  Bio::EnsEMBL::FeaturePair methods can be used

  The cigar_string contains the ungapped pieces that make up the gapped 
  alignment

  It looks like: n Matches [ x Deletes or Inserts m Matches ]*
  but a bit more condensed like "23M4I12M2D1M"
  and evenmore condensed as you can ommit 1s "23M4I12M2DM"


  To make things clearer this is how a blast HSP would be parsed

  >AK014066
         Length = 146

    Minus Strand HSPs:

    Score = 76 (26.8 bits), Expect = 1.4, P = 0.74
    Identities = 20/71 (28%), Positives = 29/71 (40%), Frame = -1

  Query:   479 GLQAPPPTPQGCRLIPPPPLGLQAPLPTLRAVGSSHHHP*GRQGSSLSSFRSSLASKASA 300
               G  APPP PQG R   P P G + P   L             + + ++  R  +A   +
  Sbjct:     7 GALAPPPAPQG-RWAFPRPTG-KRPATPLHGTARQDRQVRRSEAAKVTGCRGRVAPHVAP 64

  Query:   299 SSPHNPSPLPS 267
                  H P+P P+
  Sbjct:    65 PLTHTPTPTPT 75

  The alignment goes from 267 to 479 in sequence 1 and 7 to 75 in sequence 2 
  and the strand is -1.

  The alignment is made up of the following ungapped pieces :

  sequence 1 start 447 , sequence 2 start 7  , match length 33 , strand -1
  sequence 1 start 417 , sequence 2 start 18 , match length 27 , strand -1
  sequence 1 start 267 , sequence 2 start 27 , match length 137 , strand -1

  These ungapped pieces are made up into the following string (called a cigar 
  string) "33M3I27M3I137M" with start 267 end 479 strand -1 hstart 7 hend 75 
  hstrand 1 and feature type would be DnaPepAlignFeature

=cut


package Bio::EnsEMBL::BaseAlignFeature;

use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::FeaturePair);


=head2 new

  Arg [..]   : List of named arguments. (-cigar_string , -features) defined
               in this constructor, others defined in FeaturePair and 
               SeqFeature superclasses.  Either cigar_string or a list
               of ungapped features should be provided - not both.
  Example    : $baf = new BaseAlignFeatureSubclass(-cigar_string => '3M3I12M');
  Description: Creates a new BaseAlignFeature using either a cigarstring or
               a list of ungapped features.  BaseAlignFeature is an abstract
               baseclass and should not actually be instantiated - rather its
               subclasses should be.
  Returntype : Bio::EnsEMBL::BaseAlignFeature
  Exceptions : thrown if both feature and cigar string args are provided
               thrown if neither feature nor cigar string args are provided
  Caller     : general
  Status     : Stable

=cut

sub new {
  
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my ($cigar_string,$features) = rearrange([qw(CIGAR_STRING FEATURES)], @_);

  if (defined($cigar_string) && defined($features)) {
    throw("CIGAR_STRING or FEATURES argument is required - not both.");
  } elsif (defined($features)) {
    $self->_parse_features($features);
    
  } elsif (defined($cigar_string)) {
    $self->{'cigar_string'} = $cigar_string;
  } else {
    throw("CIGAR_STRING or FEATURES argument is required");
  }
  
  return $self;
}


=head2 new_fast

  Arg [1]    : hashref $hashref
               A hashref which will be blessed into a PepDnaAlignFeature.
  Example    : none
  Description: This allows for very fast object creation when a large number
               of PepDnaAlignFeatures needs to be created.  This is a bit of
               a hack but necessary when thousands of features need to be
               generated within a couple of seconds for web display. It is
               not recommended that this method be called unless you know what
               you are doing.  It requires knowledge of the internals of this
               class and its superclasses.
  Returntype : Bio::EnsEMBL::BaseAlignFeature
  Exceptions : none
  Caller     : none currently
  Status     : Stable

=cut

sub new_fast {
  my ($class, $hashref) = @_;

  return bless $hashref, $class;
}


=head2 cigar_string

  Arg [1]    : string $cigar_string
  Example    : ( "12MI3M" )
  Description: get/set for attribute cigar_string
               cigar_string describes the alignment. "xM" stands for 
               x matches (mismatches), "xI" for inserts into query sequence 
               (thats the ensembl sequence), "xD" for deletions 
               (inserts in the subject). an "x" that is 1 can be omitted.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub cigar_string {
  my $self = shift;
  $self->{'cigar_string'} = shift if(@_);
  return $self->{'cigar_string'};
}


=head2 alignment_length

  Arg [1]    : None
  Example    : 
  Description: return the alignment length (including indels) based on the cigar_string
  Returntype : int
  Exceptions : 
  Caller     : 
  Status     : Stable

=cut

sub alignment_length {
  my $self = shift;

  if (! defined $self->{'_alignment_length'} && defined $self->cigar_string) {
    
    my @pieces = ( $self->cigar_string =~ /(\d*[MDI])/g );
    unless (@pieces) {
      print STDERR "Error parsing cigar_string\n";
    }
    my $alignment_length = 0;
    foreach my $piece (@pieces) {
      my ($length) = ( $piece =~ /^(\d*)/ );
      if (! defined $length || $length eq "") {
        $length = 1;
      }
      $alignment_length += $length;
    }
    $self->{'_alignment_length'} = $alignment_length;
  }
  return $self->{'_alignment_length'};
}

=head2 ungapped_features

  Args       : none
  Example    : @ungapped_features = $align_feature->get_feature
  Description: converts the internal cigar_string into an array of
               ungapped feature pairs
  Returntype : list of Bio::EnsEMBL::FeaturePair
  Exceptions : cigar_string not set internally
  Caller     : general
  Status     : Stable

=cut

sub ungapped_features {
  my ($self) = @_;

  if (!defined($self->{'cigar_string'})) {
    throw("No cigar_string defined.  Can't return ungapped features");
  }

  return @{$self->_parse_cigar()};
}

=head2 strands_reversed
 
  Arg [1]    : int $strands_reversed
  Example    : none
  Description: get/set for attribute strands_reversed
               0 means that strand and hstrand are the original strands obtained
                 from the alignment program used
               1 means that strand and hstrand have been flipped as compared to
                 the original result provided by the alignment program used.
                 You may want to use the reverse_complement method to restore the
                 original strandness.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable
 
=cut

sub strands_reversed {
   my ($self, $arg) = @_;
 
   if ( defined $arg ) {
      $self->{'strands_reversed'} = $arg ;
   }

   $self->{'strands_reversed'} = 0 unless (defined $self->{'strands_reversed'});

   return $self->{'strands_reversed'};
}

=head2 reverse_complement

  Args       : none
  Example    : none
  Description: reverse complement the FeaturePair,
               modifing strand, hstrand and cigar_string in consequence
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub reverse_complement {
  my ($self) = @_;

  # reverse strand in both sequences
  $self->strand($self->strand * -1);
  $self->hstrand($self->hstrand * -1);

  # reverse cigar_string as consequence
  my $cigar_string = $self->cigar_string;
  $cigar_string =~ s/(D|I|M)/$1 /g;
  my @cigar_pieces = split / /,$cigar_string;
  $cigar_string = "";
  while (my $piece = pop @cigar_pieces) {
    $cigar_string .= $piece;
  }
  
  $self->{'strands_reversed'} = 0 unless (defined $self->{'strands_reversed'});

  if ($self->strands_reversed) {
    $self->strands_reversed(0)
  } else {
    $self->strands_reversed(1);
  }

  $self->cigar_string($cigar_string);
}



=head2 transform

  Arg  1     : String $coordinate_system_name
  Arg [2]    : String $coordinate_system_version
  Example    : $feature = $feature->transform('contig');
               $feature = $feature->transform('chromosome', 'NCBI33');
  Description: Moves this AlignFeature to the given coordinate system.
               If the feature cannot be transformed to the destination 
               coordinate system undef is returned instead.
  Returntype : Bio::EnsEMBL::BaseAlignFeature;
  Exceptions : wrong parameters
  Caller     : general
  Status     : Medium Risk
             : deprecation needs to be removed at some time

=cut

sub transform {
  my $self = shift;

  # catch for old style transform calls
  if( ref $_[0] eq 'HASH') {
    deprecate("Calling transform with a hashref is deprecate.\n" .
              'Use $feat->transfer($slice) or ' .
              '$feat->transform("coordsysname") instead.');
    my (undef, $new_feat) = each(%{$_[0]});
    return $self->transfer($new_feat->slice);
  }

  my $new_feature = $self->SUPER::transform( @_ );
  if( ! defined $new_feature or 
      $new_feature->length != $self->length) {
    my @segments = $self->project( @_ );

    return undef if( ! @segments );  

    my @ungapped;
    foreach my $f ($self->ungapped_features) {
      $f = $f->transform( @_ );
      if (defined $f) {
        push @ungapped, $f;
      } else {
        warning("Failed to transform alignment feature; " .
                "ungapped component could not be transformed");
        return undef;
      }
    }

    eval {     
      $new_feature = $self->new(-features => \@ungapped );
    };
    if ($@) {
      warning($@);
      return undef;
    }
  }

  return $new_feature;
}


=head2 _parse_cigar

  Args       : none
  Example    : none
  Description: PRIVATE (internal) method - creates ungapped features from 
               internally stored cigar line
  Returntype : list of Bio::EnsEMBL::FeaturePair
  Exceptions : none
  Caller     : ungapped_features
  Status     : Stable

=cut

sub _parse_cigar {
  my ( $self ) = @_;

  my $query_unit = $self->_query_unit();
  my $hit_unit = $self->_hit_unit();

  my $string = $self->{'cigar_string'};

  throw("No cigar string defined in object") if(!defined($string));

  my @pieces = ( $string =~ /(\d*[MDI])/g );
  #print "cigar: ",join ( ",", @pieces ),"\n";

  my @features;
  my $strand1 = $self->{'strand'} || 1;
  my $strand2 = $self->{'hstrand'}|| 1;

  my ( $start1, $start2 );

  if( $strand1 == 1 ) {
    $start1 = $self->{'start'};
  } else {
    $start1 = $self->{'end'};
  }

  if( $strand2 == 1 ) {
    $start2 = $self->{'hstart'};
  } else {
    $start2 = $self->{'hend'};
  }

  #
  # Construct ungapped blocks as FeaturePairs objects for each MATCH
  #
  foreach my $piece (@pieces) {

    my ($length) = ( $piece =~ /^(\d*)/ );
    if( $length eq "" ) { $length = 1 }
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
      throw("Internal error $query_unit $hit_unit, currently only " .
            "allowing 1 or 3 ");
    }

    if( int($mapped_length) != $mapped_length and
        ($piece =~ /M$/ or $piece =~ /D$/)) {
      throw("Internal error with mismapped length of hit, query " .
            "$query_unit, hit $hit_unit, length $length");
    }

    if( $piece =~ /M$/ ) {
      #
      # MATCH
      #
      my ( $qstart, $qend);
      if( $strand1 == 1 ) {
        $qstart = $start1;
        $qend = $start1 + $length - 1;
        $start1 = $qend + 1;
      } else {
        $qend = $start1;
        $qstart = $start1 - $length + 1;
        $start1 = $qstart - 1;
      }

      my ($hstart, $hend);
      if( $strand2 == 1 ) {
        $hstart = $start2;
        $hend = $start2 + $mapped_length - 1;
        $start2 = $hend + 1;
      } else {
        $hend = $start2;
        $hstart = $start2 - $mapped_length + 1;
        $start2 = $hstart - 1;
      }


      push @features, Bio::EnsEMBL::FeaturePair->new
        (-SLICE      => $self->{'slice'},
         -SEQNAME   => $self->{'seqname'},
         -START      => $qstart,
         -END        => $qend,
         -STRAND     => $strand1,
         -HSLICE      => $self->{'hslice'},
         -HSEQNAME   => $self->{'hseqname'},
         -HSTART     => $hstart,
         -HEND       => $hend,
         -HSTRAND    => $strand2,
         -SCORE      => $self->{'score'},
         -PERCENT_ID => $self->{'percent_id'},
         -ANALYSIS   => $self->{'analysis'},
         -P_VALUE    => $self->{'p_value'},
         -EXTERNAL_DB_ID => $self->{'external_db_id'}, 
         -HCOVERAGE   => $self->{'hcoverage'},
         -GROUP_ID    => $self->{'group_id'},
         -LEVEL_ID    => $self->{'level_id'});
      

      # end M cigar bits 
    } elsif( $piece =~ /I$/ ) {
      #
      # INSERT
      #
      if( $strand1 == 1 ) {
        $start1 += $length;
      } else {
        $start1 -= $length;
      }
    } elsif( $piece =~ /D$/ ) {
      #
      # DELETION
      #
      if( $strand2 == 1 ) {
        $start2 += $mapped_length;
      } else {
        $start2 -= $mapped_length;
      }
    } else {
      throw( "Illegal cigar line $string!" );
    }
  }

  return \@features;
}




=head2 _parse_features

  Arg  1     : listref Bio::EnsEMBL::FeaturePair $ungapped_features
  Example    : none
  Description: creates internal cigarstring and start,end hstart,hend
               entries.
  Returntype : none, fills in values of self
  Exceptions : argument list undergoes many sanity checks - throws under many
               invalid conditions
  Caller     : new
  Status     : Stable

=cut

my $message_only_once = 1;

sub _parse_features {
  my ($self,$features ) = @_;

  my $query_unit = $self->_query_unit();
  my $hit_unit = $self->_hit_unit();

  if (ref($features) ne "ARRAY") {
    throw("features must be an array reference not a [".ref($features)."]");
  }

  my $strand  = $features->[0]->strand;

  throw ('FeaturePair needs to have strand == 1 or strand == -1') if(!$strand);

  my @f;

  #
  # Sort the features on their start position
  # Ascending order on positive strand, descending on negative strand
  #
  if( $strand == 1 ) {
    @f = sort {$a->start <=> $b->start} @$features;
  } else {
    @f = sort { $b->start <=> $a->start} @$features;
  }

  my $hstrand     = $f[0]->hstrand;
  my $slice       = $f[0]->slice();
  my $hslice       = $f[0]->hslice();
  my $name        = $slice->name() if($slice);
  my $hname       = $f[0]->hseqname;
  my $score       = $f[0]->score;
  my $percent     = $f[0]->percent_id;
  my $analysis    = $f[0]->analysis;
  my $pvalue      = $f[0]->p_value();
  my $external_db_id = $f[0]->external_db_id;
  my $hcoverage   = $f[0]->hcoverage;
  my $group_id    = $f[0]->group_id;
  my $level_id    = $f[0]->level_id;

  my $seqname = $f[0]->seqname;
  # implicit strand 1 for peptide sequences
  $strand  ||= 1;
  $hstrand ||= 1;
  my $ori = $strand * $hstrand;

  throw("No features in the array to parse") if(scalar(@f) == 0);

  my $prev1; # where last feature q part ended
  my $prev2; # where last feature s part ended

  my $string;

  # Use strandedness info of query and hit to make sure both sets of 
  # start and end  coordinates are oriented the right way around.
  my $f1start;
  my $f1end;
  my $f2end;
  my $f2start;

  if ($strand == 1) {
    $f1start = $f[0]->start;
    $f1end = $f[-1]->end;
  } else {
    $f1start = $f[-1]->start;
    $f1end = $f[0]->end;
  }

  if ($hstrand == 1) {
    $f2start = $f[0]->hstart;
    $f2end = $f[-1]->hend;
  } else {
    $f2start = $f[-1]->hstart;
    $f2end = $f[0]->hend;
  }

  #
  # Loop through each portion of alignment and construct cigar string
  #

  foreach my $f (@f) {
    #
    # Sanity checks
    #

    if (!$f->isa("Bio::EnsEMBL::FeaturePair")) {
      throw("Array element [$f] is not a Bio::EnsEMBL::FeaturePair");
    }
    if( defined($f->hstrand()) && $f->hstrand() != $hstrand ) {
      throw("Inconsistent hstrands in feature array");
    }
    if( defined($f->strand()) && ($f->strand != $strand)) {
      throw("Inconsistent strands in feature array");
    }
    if ( defined($name) && $name ne $f->slice->name()) {
      throw("Inconsistent names in feature array [$name - ".
            $f->slice->name()."]");
    }
    if ( defined($hname) && $hname ne $f->hseqname) {
      throw("Inconsistent hit names in feature array [$hname - ".
            $f->hseqname . "]");
    }
    if ( defined($score) && $score ne $f->score) {
      throw("Inconsisent scores in feature array [$score - " .
            $f->score . "]");
    }
    if (defined($f->percent_id) && $percent ne $f->percent_id) {
      throw("Inconsistent pids in feature array [$percent - " .
            $f->percent_id . "]");
    }
    if(defined($pvalue) && $pvalue != $f->p_value()) {
      throw("Inconsistant p_values in feature arraw [$pvalue " .
            $f->p_value() . "]");
    }
    if($seqname && $seqname ne $f->seqname){
      throw("Inconsistent seqname in feature array [$seqname - ".
            $f->seqname . "]");
    }
    my $start1 = $f->start;      #source sequence alignment start
    my $start2 = $f->hstart();   #hit sequence alignment start

    #
    # More sanity checking
    #
    if (defined($prev1)) {
      if ( $strand == 1 ) {
        if ($f->start < $prev1) {
          throw("Inconsistent coords in feature array (forward strand).\n" .
		       "Start [".$f->start()."] in current feature should be greater " .
           "than previous feature end [$prev1].");
        }
      } else {
        if ($f->end > $prev1) {
          throw("Inconsistent coords in feature array (reverse strand).\n" .
                "End [".$f->end() ."] should be less than previous feature " .
                "start [$prev1].");
        }
      }
    }

    my $length = ($f->end - $f->start + 1);    #length of source seq alignment
    my $hlength = ($f->hend - $f->hstart + 1); #length of hit seq alignment

    # using multiplication to avoid rounding errors, hence the
    # switch from query to hit for the ratios

    #
    # Yet more sanity checking
    #
    if($query_unit > $hit_unit){
      # I am going to make the assumption here that this situation will 
      # only occur with DnaPepAlignFeatures, this may not be true
      my $query_p_length = sprintf "%.0f", ($length/$query_unit);
      my $hit_p_length = sprintf "%.0f", ($hlength * $hit_unit);
      if( $query_p_length != $hit_p_length) {
        throw( "Feature lengths not comparable Lengths:" .$length .
               " " . $hlength . " Ratios:" . $query_unit . " " .
               $hit_unit );
      }
    } else{
      my $query_d_length = sprintf "%.0f", ($length*$hit_unit);
      my $hit_d_length = sprintf "%.0f", ($hlength * $query_unit);
      if( $length * $hit_unit != $hlength * $query_unit ) {
        throw( "Feature lengths not comparable Lengths:" . $length .
               " " . $hlength . " Ratios:" . $query_unit . " " .
               $hit_unit );
      }
    }

    my $hlengthfactor = ($query_unit/$hit_unit);

    #
    # Check to see if there is an I type (insertion) gap:
    #   If there is a space between the end of the last source sequence 
    #   alignment and the start of this one, then this is an insertion
    #

    my $insertion_flag = 0;
    if( $strand == 1 ) {
      if( ( defined $prev1 ) && ( $f->start > $prev1 + 1  )) {

        #there is an insertion
        $insertion_flag = 1;
        my $gap = $f->start - $prev1 - 1;
        if( $gap == 1 ) {
          $gap = ""; # no need for a number if gap length is 1
        }
        $string .= "$gap"."I";

      }

      #shift our position in the source seq alignment
      $prev1 = $f->end();
    } else {

      if(( defined $prev1 ) && ($f->end + 1 < $prev1 )) {

        #there is an insertion
        $insertion_flag = 1;
        my $gap = $prev1 - $f->end() - 1;
        if( $gap == 1 ) {
          $gap = ""; # no need for a number if gap length is 1
        }
        $string .= "$gap"."I";
      }

      #shift our position in the source seq alignment
      $prev1 = $f->start();
    }

    #
    # Check to see if there is a D type (deletion) gap
    #   There is a deletion gap if there is a space between the end of the
    #   last portion of the hit sequence alignment and this one
    #
    if( $hstrand == 1 ) {
      if((  defined $prev2 ) && ( $f->hstart() > $prev2 + 1 )) {

        #there is a deletion
        my $gap = $f->hstart - $prev2 - 1;
        my $gap2 = int( $gap * $hlengthfactor + 0.5 );
	
        if( $gap2 == 1 ) {
          $gap2 = "";  # no need for a number if gap length is 1
        }
        $string .= "$gap2"."D";

        #sanity check,  Should not be an insertion and deletion
        if($insertion_flag) {
          if ($message_only_once) {
            warning("Should not be an deletion and insertion on the " .
                    "same alignment region. cigar_line=$string\n");
            $message_only_once = 0;
          }
        }
      }
      #shift our position in the hit seq alignment
      $prev2 = $f->hend();

     } else {
      if( ( defined $prev2 ) && ( $f->hend() + 1 < $prev2 )) {

        #there is a deletion
        my $gap = $prev2 - $f->hend - 1;
        my $gap2 = int( $gap * $hlengthfactor + 0.5 );
	
        if( $gap2 == 1 ) {
          $gap2 = "";  # no need for a number if gap length is 1
        }
        $string .= "$gap2"."D";

        #sanity check,  Should not be an insertion and deletion
        if($insertion_flag) {
          if ($message_only_once) {
            warning("Should not be an deletion and insertion on the " .
                    "same alignment region. prev2 = $prev2; f->hend() = " .
                    $f->hend() . "; cigar_line = $string;\n");
            $message_only_once = 0;
          }
        }
      }
      #shift our position in the hit seq alignment
      $prev2 = $f->hstart();
    }

    my $matchlength = $f->end() - $f->start() + 1;
    if( $matchlength == 1 ) {
      $matchlength = "";
    }
    $string .= $matchlength."M";
  }

  $self->{'start'}      = $f1start;
  $self->{'end'}        = $f1end;
  $self->{'seqname'}    = $seqname;
  $self->{'strand'}     = $strand;
  $self->{'score'}      = $score;
  $self->{'percent_id'} = $percent;
  $self->{'analysis'}   = $analysis;
  $self->{'slice'}      = $slice;
  $self->{'hslice'}     = $hslice;
  $self->{'hstart'}     = $f2start;
  $self->{'hend'}       = $f2end;
  $self->{'hstrand'}    = $hstrand;
  $self->{'hseqname'}   = $hname;
  $self->{'cigar_string'} = $string;
  $self->{'p_value'}      = $pvalue;
  $self->{'external_db_id'} = $external_db_id; 
  $self->{'hcoverage'}    = $hcoverage; 
  $self->{'group_id'}     = $group_id;
  $self->{'level_id'}     = $level_id;
}






=head2 _hit_unit

  Args       : none
  Example    : none
  Description: abstract method, overwrite with something that returns
               one or three
  Returntype : int 1,3
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub _hit_unit {
  my $self = shift;
  throw( "Abstract method call!" );
}



=head2 _query_unit

  Args       : none
  Example    : none
  Description: abstract method, overwrite with something that returns
               one or three
  Returntype : int 1,3
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub _query_unit {
  my $self = shift;
  throw( "Abstract method call!" );
}




1;
