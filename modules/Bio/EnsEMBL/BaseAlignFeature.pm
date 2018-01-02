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
    -cigar_string => '10M3D5M2I',
    -align_type   => 'ensembl'
  );

  where $analysis is a Bio::EnsEMBL::Analysis object.

  Alternatively if you have an array of ungapped features:

    my $feat =
      new Bio::EnsEMBL::DnaPepAlignFeature( -features => \@features );

  where @features is an array of Bio::EnsEMBL::FeaturePair objects.

  There is a method to (re)create ungapped features from the cigar_string:

    my @ungapped_features = $feat->ungapped_features();

  where @ungapped_features is an array of Bio::EnsEMBL::FeaturePair's.

  Bio::EnsEMBL::BaseAlignFeature inherits from:
    Bio::EnsEMBL::FeaturePair, which in turn inherits from:
      Bio::EnsEMBL::Feature,
  thus methods from both parent classes are available.


  The cigar_string is a condensed representation of the matches and gaps
  which make up the gapped alignment (where CIGAR stands for
  Concise Idiosyncratic Gapped Alignment Report).

  CIGAR format is: n <matches> [ x <deletes or inserts> m <matches> ]*
  where M = match, D = delete, I = insert; n, m are match lengths;
  x is delete or insert length.

  Spaces are omitted, thus: "23M4I12M2D1M"
  as are counts for any lengths of 1, thus 1M becomes M: "23M4I12M2DM"


  To make things clearer this is how a blast HSP would be parsed:

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

  The alignment goes from 479 down to 267 in the query sequence on the reverse
  strand, and from 7 to 75 in the subject sequence.

  The alignment is made up of the following ungapped pieces:

  query_seq start 447 , sbjct_seq hstart  7 , match length  33 , strand -1
  query_seq start 417 , sbjct_seq hstart 18 , match length  27 , strand -1
  query_seq start 267 , sbjct_seq hstart 27 , match length 147 , strand -1

  When assembled into a DnaPepAlignFeature where:
    (seqname, start, end, strand) refer to the query sequence,
    (hseqname, hstart, hend, hstrand) refer to the subject sequence,
  these ungapped pieces are represented by the cigar string:
    33M3I27M3I147M
  with start 267, end 479, strand -1, and hstart 7, hend 75, hstrand 1.

=head1 CAVEATS

  AlignFeature cigar strings have the opposite 'sense'
  ('D' and 'I' swapped) compared with Exonerate cigar strings.

  Exonerate modules in Bio::EnsEMBL::Analysis use this convention:

   The longer genomic sequence specified by:
      exonerate:    target
      AlignFeature: (sequence, start, end, strand)

   A shorter sequence (such as EST or protein) specified by:
      exonerate:    query
      AlignFeature: (hsequence, hstart, hend, hstrand)

  The resulting AlignFeature cigar strings have 'D' and 'I'
  swapped compared with the Exonerate output, i.e.:

    exonerate:    M 123 D 1 M 11 I 1 M 39
    AlignFeature: 123MI11MD39M

=head1 METHODS

=cut


package Bio::EnsEMBL::BaseAlignFeature;

use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::FeaturePair);


=head2 new

  Arg [..]   : List of named arguments. (-cigar_string , -features, -align_type) defined
               in this constructor, others defined in FeaturePair and 
               SeqFeature superclasses.  Either cigar_string or a list
               of ungapped features should be provided - not both.
  Example    : $baf = new BaseAlignFeatureSubclass(-cigar_string => '3M3I12M', -align_type => 'ensembl');
  Description: Creates a new BaseAlignFeature using either a cigar string or
               a list of ungapped features.  BaseAlignFeature is an abstract
               baseclass and should not actually be instantiated - rather its
               subclasses should be.
  Returntype : Bio::EnsEMBL::BaseAlignFeature
  Exceptions : thrown if both feature and cigar string args are provided
               thrown if neither feature nor cigar string args are provided
               warn if cigar string is provided without cigar type
  Caller     : general
  Status     : Stable

=cut

sub new {
  
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my ($cigar_string,$align_type,$features) = rearrange([qw(CIGAR_STRING ALIGN_TYPE FEATURES)], @_);

  if (defined($align_type)) {
    $self->{'align_type'} = $align_type;
  } else {
    warning("No align_type provided, using ensembl as default");
    $self->{'align_type'} = 'ensembl';
  }

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


=head2 cigar_string

  Arg [1]    : string $cigar_string
  Example    : $feature->cigar_string( "12MI3M" );
  Description: get/set for attribute cigar_string.
               cigar_string describes the alignment:
                 "xM" stands for x matches (or mismatches),
                 "xI" for x inserts into the query sequence,
                 "xD" for x deletions from the query sequence
                 where the query sequence is specified by (seqname, start, ...)
                 and the subject sequence by (hseqname, hstart, ...).
               An "x" that is 1 can be omitted.
               See the SYNOPSIS for an example.
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


=head2 align_type

  Arg [1]    : type $align_type
  Example    : $feature->align_type( "ensembl" );
  Description: get/set for attribute align_type.
               align_type specifies which cigar string 
               is used to describe the alignment:
               The default is 'ensembl'
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub align_type {
  my $self = shift;
  $self->{'align_type'} = shift if(@_);
  return $self->{'align_type'};
}


=head2 alignment_length

  Arg [1]    : None
  Description: return the alignment length (including indels) based on the cigar_string
  Returntype : int
  Exceptions : 
  Caller     : 
  Status     : Stable

=cut

sub alignment_length {
  my $self = shift;

  if ($self->{'align_type'} eq 'ensembl') {
    return $self->_ensembl_cigar_alignment_length();
  } else {
    throw("No alignment_length method available for " . $self->{'align_type'});
  }

}


=head2 _ensembl_cigar_alignment_length

  Arg [1]    : None
  Description: return the alignment length (including indels) based on the cigar_string
  Returntype : int
  Exceptions :
  Caller     :
  Status     : Stable

=cut

sub _ensembl_cigar_alignment_length {
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
  Description: reverse complement the FeaturePair based on the cigar type
               modifing strand, hstrand and cigar_string in consequence
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub reverse_complement {
  my ($self) = @_;

  if ($self->{'align_type'} eq 'ensembl') {
    return $self->_ensembl_reverse_complement();
  } else {
    throw("no reverse_complement method implemented for " . $self->{'align_type'});
  }
}

=head2 _ensembl_reverse_complement

  Args       : none
  Description: reverse complement the FeaturePair for ensembl cigar string,
               modifing strand, hstrand and cigar_string in consequence
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub _ensembl_reverse_complement {
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

  my $new_feature = $self->SUPER::transform(@_);
  if ( !defined($new_feature)
    || $new_feature->length() != $self->length() )
  {
    my @segments = @{ $self->project(@_) };

    if ( !@segments ) {
      return undef;
    }

    my @ungapped;
    foreach my $f ( $self->ungapped_features() ) {
      $f = $f->transform(@_);
      if ( defined($f) ) {
        push( @ungapped, $f );
      } else {
        warning( "Failed to transform alignment feature; "
            . "ungapped component could not be transformed" );
        return undef;
      }
    }

    eval { $new_feature = $self->new( -features => \@ungapped ); };

    if ($@) {
      warning($@);
      return undef;
    }
  } ## end if ( !defined($new_feature...))

  return $new_feature;
}


=head2 _parse_ensembl_cigar

  Args       : none
  Description: PRIVATE (internal) method - creates ungapped features from 
               internally stored cigar line in ensembl format
  Returntype : list of Bio::EnsEMBL::FeaturePair
  Exceptions : none
  Caller     : ungapped_features
  Status     : Stable

=cut

sub _parse_ensembl_cigar {
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


sub _parse_cigar {
  my $self = shift;
  if ($self->{'align_type'} eq 'ensembl') {
    return $self->_parse_ensembl_cigar();
  } else {
    throw("No parsing method implemented for " . $self->{'align_type'});
  }
}




=head2 _parse_features

  Arg  [1]   : listref Bio::EnsEMBL::FeaturePair $ungapped_features
  Description: creates internal cigar_string and start,end hstart,hend
               entries.
  Returntype : none, fills in values of self
  Exceptions : argument list undergoes many sanity checks - throws under many
               invalid conditions
  Caller     : new
  Status     : Stable

=cut

sub _parse_features {
  my ($self, $features) = @_;
  if ($self->{'align_type'} eq 'ensembl') {
    $self->_parse_ensembl_features($features);
  } else {
    throw("No _parse_features method implemented for " . $self->{'align_type'});
  }
}

my $message_only_once = 1;

sub _parse_ensembl_features {
  my ($self,$features ) = @_;

  if (ref($features) ne "ARRAY") {
    throw("features must be an array reference not a [".ref($features)."]");
  } elsif (scalar(@$features) == 0) {
    throw("features array must not be empty");
  }

  my $query_unit = $self->_query_unit();
  my $hit_unit = $self->_hit_unit();

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
  my $name        = $slice ? $slice->name() : undef;
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
