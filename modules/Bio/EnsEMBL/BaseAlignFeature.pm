# EnsEMBL module for storing dna-protein pairwise alignments
#
# Cared for by Michele Clamp <michele@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself
#
=pod

=head1 NAME

Bio::EnsEMBL::DnaPepAlignFeature - Ensembl specific dna-protein 
                                   pairwise alignment feature

=head1 SYNOPSIS

  my $feat = new Bio::EnsEMBL::DnaPepAlignFeature(-seqname => 'myseq',
                                                  -start   => 100,
                                                  -end     => 120,
                                                  -strand  => 1,
                                                  -hstart  => 200,
                                                  -hend    => 220,
                                                  -analysis    => $analysis,
                                                  -cigar_string => '');

  Alternatively if you have an array of ungapped features

  my $feat = new Bio::EnsEMBL::DnaPepAlignFeature(-features => \@features);

  Where @features is an array of Bio::EnsEMBL::FeaturePair

  There is a method to manipulate the cigar_string into ungapped features

  my @ungapped_features = $feat->ungapped_features;

  This converts the cigar string into an array of Bio::EnsEMBL::FeaturePair

  $analysis is a Bio::EnsEMBL::Analysis object
  
  Bio::EnsEMBL::SeqFeature methods can be used
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
use Bio::EnsEMBL::SeqFeature;

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::FeaturePair);


=head2 new

  Arg [..]   : List of named arguments. (-cigar_string , -features) defined
               in this constructor, others defined in FeaturePair and 
               SeqFeature superclasses.  
  Example    : 
  Description: Creates a new BaseAlignFeature using either a cigarstring or
               a list of ungapped features.
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
    #print "calling new with @args\n";
    my ($cigar_string,$features) 
      = $self->_rearrange([qw(CIGAR_STRING FEATURES)], @args);

    if (defined($cigar_string) && defined($features)) {
      $self->throw("Can't input cigar_string and an array of features.");
    } elsif (defined($features)) {
      $self->_parse_features($features);
    } elsif (!defined($cigar_string)) {
      $self->throw("Must have a cigar string defining the alignment");
    } 
    
    $self->cigar_string($cigar_string);
    
    return $self;
}


=head2 cigar_string

  Arg [1]    : string $cigar_string
  Example    : ( "12MI3M" )
  Description: get/set for attribute cigar_string
               cigar_string describes the alignment. "xM" stands for 
               x matches (mismatches), "xI" for inserts into query sequence 
               (thats the ensembl sequence), "xD" for deletions (inserts in the 
               subject). an "x" that is 1 can be omitted.
  Returntype : string
  Exceptions : throws on a get without previous set
  Caller     : general

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


=head2 ungapped_features

  Args       : none
  Example    : none
  Description: converts the internal cigar_string into an array of
               ungapped feature pairs
  Returntype : list of Bio::EnsEMBL::FeaturePair
  Exceptions : cigar_string not set internally
  Caller     : general

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



=head2 transform

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Example    : none
  Description: if argument is given, transforms this feature into the slice
               coordinate system, invalidating this one.
               if no argument is given, transforms this feature into raw contig
               coordinates, invalidating this one.
               The process can produce more than one feature so we return an array. 
  Returntype : list of Bio::EnsEMBL::BaseAlignFeature
  Exceptions : none
  Caller     : general

=cut


sub transform{
  my ($self, $slice) = @_;
  #print "transforming ".$self." to ".$slice." coords\n";
  if( ! defined $slice ) {
    #Since slice arg is not defined -  we want raw contig coords
    if(( defined  $self->contig ) && 
       ( $self->contig->isa( "Bio::EnsEMBL::RawContig" )) ) {
   #   print STDERR "BaseAlignFeature::transform, you are already apparently in rawcontig coords so why try to transform to them\n";
      #we are already in rawcontig coords, nothing needs to be done
      return $self;
    } else {
      #transform to raw_contig coords from Slice coords
      my @array =  $self->_transform_to_rawcontig();
    #  print "transform to  rawcontig has returned ".@array." features\n";
    #  foreach my $a(@array){
#	print "feature is ".$a."\n";
#      }
      return @array;
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

=head2 dbID

  Arg [1]    : int $dbID
  Example    : none
  Description: get/set for the database internal id
  Returntype : int
  Exceptions : none
  Caller     : general, set from adaptor on store

=cut


sub dbID{
  my ($self, $arg) = @_;

  if(defined $arg){
    $self->{_database_id} = $arg;
  }

  return $self->{_database_id}; 

}

=head2 adaptor

  Arg [1]    : Bio::EnsEMBL::DBSQL::BaseAlignFeatureAdaptor $adaptor
  Example    : none
  Description: get/set for this objects Adaptor
  Returntype : Bio::EnsEMBL::DBSQL::BaseAlignFeatureAdaptor
  Exceptions : none
  Caller     : general, set from adaptor on store

=cut

sub adaptor {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'adaptor'} = $value;
    }
    return $self->{'adaptor'};

}


=head2 _parse_cigar

  Args       : none
  Example    : none
  Description: creates ungapped features from internally stored cigar line
  Returntype : list of Bio::EnsEMBL::FeaturePair
  Exceptions : none
  Caller     : ungapped_features

=cut

sub _parse_cigar {
  my ( $self ) = @_;

  my $query_unit = $self->_query_unit();
  my $hit_unit = $self->_hit_unit();

  
  my $string = $self->cigar_string;
 
  if (!defined($string)) {
    $self->throw("No cigar string defined in object.  This should be caught by the cigar_string method and never happen");
  }

  
  my @pieces = ( $string =~ /(\d*[MDI])/g );
  #print "cigar: ",join ( ",", @pieces ),"\n";
  
  my @features;
  my $strand1 = $self->strand() || 1;
  my $strand2 = $self->hstrand() || 1;

  my ( $start1, $start2 );
  
  if( $strand1 == 1 ) {
    $start1 = $self->start();
  } else {
    $start1 = $self->end();
  }

  if( $strand2 == 1 ) {
    $start2 = $self->hstart();
  } else {
    $start2 = $self->hend();
  }


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
      $self->throw("Internal error $query_unit $hit_unit, currently only allowing 1 or 3 ");
    }
    
    if( int($mapped_length) != $mapped_length ) {
      $self->throw("Internal error with mismapped length of hit, query $query_unit, hit $hit_unit, length $length");
    }

    if( $piece =~ /M$/ ) {
      my $feature1 = new Bio::EnsEMBL::SeqFeature();
      
      my ( $a, $b );
      if( $strand1 == 1 ) {
        $a = $start1;
        $b = $start1 + $length - 1;
        $start1 = $b + 1;
      } else {
        $b = $start1;
        $a = $start1 - $length + 1;
        $start1 = $a - 1;
      }
      
      $feature1->start($a);
      $feature1->end  ($b);
      $feature1->strand($self->strand() );
      $feature1->score($self->score);
      $feature1->seqname($self->seqname);
      $feature1->phase($self->phase);
      $feature1->p_value($self->p_value);
      $feature1->percent_id($self->percent_id);

      my $feature2 = new Bio::EnsEMBL::SeqFeature();
      if( $strand2 == 1 ) {
        $a = $start2;
        $b = $start2 + $mapped_length - 1;
        $start2 = $b + 1;
      } else {
        $b = $start2;
        $a = $start2 - $mapped_length + 1;
        $start2 = $a - 1;
      }

      $feature2->start($a);
      $feature2->end($b);
      $feature2->strand($self->hstrand());
      $feature2->score($self->score);
      $feature2->seqname($self->hseqname);
      $feature2->phase($self->phase);
      $feature2->p_value($self->p_value);
      $feature2->percent_id($self->percent_id);
    
      my $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $feature1,
                                             -feature2 => $feature2);

 
      $fp->analysis($self->analysis);

      push(@features,$fp);
      # end M cigar bits 
    } elsif( $piece =~ /I$/ ) {
      if( $strand1 == 1 ) {
        $start1 += $length;
      } else {
        $start1 -= $length;
      }
    } elsif( $piece =~ /D$/ ) {
      if( $strand2 == 1 ) {
        $start2 += $mapped_length;
      } else {
        $start2 -= $mapped_length;
      }
    } else {
      $self->throw( "Illegal cigar line $string!" );
    } 
      
  }
  # should the features be sorted ?
  # 
  return @features;
}




=head2 _parse_features

  Arg  1     : listref Bio::EnsEMBL::FeaturePair $ungapped_features
  Example    : none
  Description: creates internal cigarstring and start,end hstart,hend
               entries.
  Returntype : none, fills in values of self
  Exceptions : argument list is sanity checked
  Caller     : new

=cut

sub _parse_features {
  my ($self,$features ) = @_;

  my $query_unit = $self->_query_unit();
  my $hit_unit = $self->_hit_unit();
  my $verbose = 0;
  if($features->[0]->hseqname eq 'CE25688'){
    $verbose = 1;
  }
  if (ref($features) ne "ARRAY") {
    $self->throw("features must be an array reference not a [" . 
		 ref($features) . "]");
  }

  my $strand     = $features->[0]->strand;
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
  #print STDERR $f[0]->gffstring."\n";  
  my $hstrand     = $f[0]->hstrand;
  my $name        = $f[0]->seqname;
  my $hname       = $f[0]->hseqname;
  my $score       = $f[0]->score;
  my $percent     = $f[0]->percent_id;
  my $pvalue      = $f[0]->p_value;
  my $analysis    = $f[0]->analysis;
  my $phase       = $f[0]->phase;
 
  # implicit strand 1 for peptide sequences
  ( defined $strand ) || ( $strand = 1 );
  ( defined $hstrand ) || ( $hstrand = 1 );
  my $ori = $strand * $hstrand;
    
    
  if (scalar(@f) == 0) {
    $self->throw("No features in the array to parse");
  }

  my $prev1; # where last feature q part ended
  my $prev2; # where last feature s part ended
    
  my $string;

  # Use strandedness info of query and hit to make sure both sets of start and end 
  # coordinates are oriented the right way around.
  my $f1start;
  my $f1end;
  my $f2end;
  my $f2start;

  if($strand == 1 && $hstrand == 1){
    $f1start = $f[0]->start;
    $f1end = $f[-1]->end;
    $f2start = $f[0]->start;
    $f2end = $f[-1]->end;
  }elsif($strand == -1 && $hstrand == -1){
    $f1start = $f[-1]->start;
    $f1end = $f[0]->end;
    $f2start = $f[-1]->start;
    $f2end = $f[0]->end;
  }elsif($strand == 1 && $hstrand == -1){
    $f1start = $f[0]->start;
    $f1end = $f[-1]->end;
    $f2start = $f[0]->start;
    $f2end = $f[-1]->end;
  }elsif($strand == -1 && $hstrand == 1){
    $f1start = $f[-1]->start;
    $f1end = $f[0]->end;
    $f2start = $f[-1]->start;
    $f2end = $f[0]->end;
  }

    
  #
  # Loop through each portion of alignment and construct cigar string
  #
  foreach my $f (@f) {
    #
    # Sanity checks
    #
    #print STDERR "supporting feature to convert ".$f->gffstring."\n" if($verbose);
    if (!$f->isa("Bio::EnsEMBL::FeaturePair")) {
      $self->throw("Array element [$f] is not a Bio::EnsEMBL::FeaturePair");
    }
    if( defined $f->hstrand() ) {
      if ($f->hstrand != $hstrand) {
	$self->throw("Inconsistent hstrands in feature array");
      }
    }
    if( defined $f->strand() ) {
      if ($f->strand != $strand) {
	$self->throw("Inconsistent strands in feature array");
      }
    }
    if ($name ne $f->seqname) {
      $self->throw("Inconsistent names in feature array [$name - " . 
		   $f->seqname . "]");
    }
    if ($hname ne $f->hseqname) {
      $self->throw("Inconsistent names in feature array [$hname - " . 
		   $f->hseqname . "]");
    }
    if ($score ne $f->score) {
      $self->throw("Inconsisent scores in feature array [$score - " . 
		   $f->score . "]");
    }
    if ($percent ne $f->percent_id) {
      $self->throw("Inconsistent pids in feature array [$percent - " . 
		   $f->percent_id . "]");
    }
    
    my $start1 = $f->start;      #source sequence alignment start
    my $start2 = $f->hstart();   #hit sequence alignment start
 #   print STDERR "start1 ".$start1." start2 ".$start2."\n"  if($verbose);
    #
    # More sanity checking
    #
    if (defined($prev1)) {
      if ( $strand == 1 ) {
	if ($f->start < $prev1) {
	  $self->throw("Inconsistent coordinates feature is forward strand " .
		       "hstart in current feature should be greater than " .
		       "hend in previous feature " . $f->start . " < " .
		       $prev1."\n");
	}
      } else {
	if ($f->end > $prev1) {
	  $self->throw("Inconsistent coordinates in feature array feature " .
		       "is reverse strand hend should be less than previous " .
		       "hstart " . $f->end . " > $prev1");
	}
      }
    }

    my $length = ($f->end - $f->start + 1); #length of source seq alignment
    my $hlength = ($f->hend - $f->hstart + 1); #length of hit seq alignment
  #  print STDERR "length ".$length." hlength ".$hlength."\n"  if($verbose);
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
	print STDERR $length."/".$query_unit." ".$hlength."*".$hit_unit."\n";
	$self->throw( "Feature lengths not comparable Lengths:" .$length . 
		      " " . $hlength . " Ratios:" . $query_unit . " " . 
		      $hit_unit );
      }
    } else{
      my $query_d_length = sprintf "%.0f", ($length*$hit_unit);
      my $hit_d_length = sprintf "%.0f", ($hlength * $query_unit);
      if( $length * $hit_unit != $hlength * $query_unit ) {
	print STDERR $length."*".$hit_unit." ".$hlength."*".$query_unit."\n";
	$self->throw( "Feature lengths not comparable Lengths:" . $length . 
		      " " . $hlength . " Ratios:" . $query_unit . " " . 
		      $hit_unit );
      }
    }

    my $hlengthfactor = ($query_unit/$hit_unit);
   # print STDERR "hlengthfactor ".$hlengthfactor."\n" if($verbose);
    # if( $query_unit == 1 && $hit_unit == 3 ) {
    #	$hlengthfactor = (1/3);
    #     }
    #     if( $query_unit == 3 && $hit_unit == 1 ) {
    #	$hlengthfactor = 3;
    #      }
      

    #
    # Check to see if there is an I type (insertion) gap:
    #   If there is a space between the end of the last source sequence 
    #   alignment and the start of this one, then this is an insertion
    #

    my $insertion_flag = 0;
    if( $strand == 1 ) {
      if( ( defined $prev1 ) && ( $f->start > $prev1 + 1  )) {
	#print STDERR "there is an insertion\n" if($verbose);
	#there is an insertion
	$insertion_flag = 1;
	my $gap = $f->start - $prev1 - 1;
	if( $gap == 1 ) {
	  $gap = ""; # no need for a number if gap length is 1
	}
	$string .= "$gap"."I";
	#print STDERR "cigar stands at ".$string."\n" if($verbose);
      }

      #shift our position in the source seq alignment
      $prev1 = $f->end();
    } else {
      #print STDERR "there is a insertion\n";
      if(( defined $prev1 ) && ($f->end + 1 < $prev1 )) {

	#there is an insertion
	$insertion_flag = 1;
	my $gap = $prev1 - $f->end() - 1;
	if( $gap == 1 ) {
	  $gap = ""; # no need for a number if gap length is 1
	}
	$string .= "$gap"."I";
	#print STDERR "cigar stands at ".$string."\n" if($verbose);
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
	#print STDERR "there is a deletion\n" if($verbose);
	#there is a deletion
	my $gap = $f->hstart - $prev2 - 1;
	my $gap2 = int( $gap * $hlengthfactor + 0.05 );
	
	if( $gap2 == 1 ) {
	  $gap2 = "";  # no need for a number if gap length is 1
	}
	$string .= "$gap2"."D";
	#print STDERR "cigar stands at ".$string."\n"  if($verbose);
	#sanity check,  Should not be an insertion and deletion
	if($insertion_flag) {
	  $self->warn("Should not be an deletion and insertion on the " .
		       "same alignment region. cigar_line=$string\n");
	} 
      } 
      #shift our position in the hit seq alignment
      $prev2 = $f->hend();

     } else {
      if( ( defined $prev2 ) && ( $f->hend() + 1 < $prev2 )) {
	#print STDERR "there is a deletion\n" if($verbose);
	#there is a deletion
	my $gap = $prev2 - $f->hend - 1;
	my $gap2 = int( $gap * $hlengthfactor + 0.05 );
	
	if( $gap2 == 1 ) {
	  $gap2 = "";  # no need for a number if gap length is 1
	}
	$string .= "$gap2"."D";
	#print STDERR "cigar stands at ".$string."\n" if($verbose);
	#sanity check,  Should not be an insertion and deletion
	if($insertion_flag) {
	  $self->throw("Should not be an deletion and insertion on the " .
		       "same alignment region. prev2 = $prev2; f->hend() = " .
		       $f->hend() . "; cigar_line = $string;\n");
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
  #  print STDERR "cigar stands at ".$string."\n" if($verbose);
    #print STDERR "finished with this feature\n\n";
  }
#  print STDERR "creating align feature start ".$f1start." end ".$f1end." strand ".$strand." score ".$score." percent id ".$percent." pvalue ".$pvalue." seqname ".$name." phase ".$phase." analysis ".$analysis." hstart ".$f2start," hend ".$f2end." hstrand ".$hstrand." hid ".$hname."\n"; 
  if(!$score){
    #$self->warn("score is not set assume its 1");
    $score = 1;
  } 
  my $feature1 = new Bio::EnsEMBL::SeqFeature();
  
  $feature1->start($f1start);
  $feature1->end  ($f1end);
  $feature1->strand($strand);
  $feature1->score($score);
  $feature1->percent_id($percent);
  $feature1->p_value($pvalue);
  $feature1->seqname($name);
  $feature1->phase($phase);
  $feature1->analysis($analysis);
  #print STDERR "checking feature1 ".$feature1->gffstring."\n";
  $feature1->validate;
  
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
  $feature2->validate;
  $self->feature1($feature1);
  $self->feature2($feature2);
  $self->cigar_string($string);
  #print STDERR "finished ".$self->gffstring."\t ".$self->cigar_string."\n\n" if($verbose);

  #print STDERR "\n\n";
}



sub _transform_to_slice{
  my ($self, $slice ) = @_;
  $self->throw( "implented soon :-)" );
}



sub _transform_to_rawcontig{
  my ($self) = @_;
  #print STDERR "transforming to raw contig coord\n\n";
  if(!$self->contig){
    $self->throw("can't transform coordinates of ".$self." without some sort of contig defined");
  }
  #my $mapper = $self->contig->adaptor->db->get_AssemblyMapperAdaptor->fetch_by_type( $self->contig()->assembly_type() );
 
  my $rcAdaptor = $self->adaptor()->db()->get_RawContigAdaptor();
  #my $global_start = $self->contig->chr_start();
  my @out;

  my @features = $self->_parse_cigar();
  my @mapped_features;
  my %rc_features;

  foreach my $f(@features){
    my @new_features = $self->_transform_feature_to_rawcontig($f);
    push(@mapped_features, @new_features);
  }

  foreach my $mf(@mapped_features){
    my $contig_id = $mf->seqname;
    if(!$rc_features{$contig_id}){
      $rc_features{$contig_id} = [];
      push(@{$rc_features{$contig_id}}, $mf);
    }else{
      push(@{$rc_features{$contig_id}}, $mf);
    }
  }

  foreach my $contig_id(keys(%rc_features)){
    my $outputf = $self->new( -features => \@{$rc_features{$contig_id}} );
    $outputf->analysis( $self->analysis() );
    $outputf->score( $self->score() );
    $outputf->percent_id( $self->percent_id() );
    $outputf->p_value( $self->p_value() );
    push(@out, $outputf);
  }
  return @out;

}



sub _transform_between_slices {
  my ( $self, $to_slice ) = @_;
}





=head2 _hit_unit

  Args       : none
  Example    : none
  Description: abstract method, overwrite with something that returns
               one or three
  Returntype : int 1,3
  Exceptions : none
  Caller     : internal

=cut

sub _hit_unit {
  my $self = shift;
  $self->throw( "Abstract method call!" );
}

=head2 _query_unit

  Args       : none
  Example    : none
  Description: abstract method, overwrite with something that returns
               one or three
  Returntype : int 1,3
  Exceptions : none
  Caller     : internal

=cut

sub _query_unit {
  my $self = shift;
  $self->throw( "Abstract method call!" );
}


sub _transform_feature_to_rawcontig{
  my($self, $feature) =  @_;

  my $verbose = 0;
  #$verbose = 1 if($feature->hseqname eq 'CE25688');

  #print STDERR "transforming ".$feature->gffstring."\n" if($verbose);

  if(!$self->contig){
    $self->throw("can't transform coordinates of ".$self." without some sort of contig defined");
  }
  my $mapper = $self->contig->adaptor->db->get_AssemblyMapperAdaptor->fetch_by_type( $self->contig()->assembly_type() );
  
  my $rcAdaptor = $self->adaptor()->db()->get_RawContigAdaptor();
  my $global_start = $self->contig->chr_start();
  my @out;
  my @mapped = $mapper->map_coordinates_to_rawcontig
    (
     $self->contig()->chr_name(),
     $feature->start()+$global_start-1,
     $feature->end()+$global_start-1,
     $feature->strand()*$self->contig()->strand()
    );

  if( ! @mapped ) {
    $self->throw( "couldn't map ".$self."\n" );
    return $self;
  }
  if( scalar( @mapped ) > 1 ) {
    my $hit_start = $feature->hstart;
    #print STDERR " feature is being mapped across multiple contigs ".$feature->gffstring."\n";
  SPLIT: for( my $i=0; $i <= $#mapped; $i++ ) {
      if($mapped[$i]->isa("Bio::EnsEMBL::Mapper::Gap")){
	$self->warn("piece of evidence lies on gap\n");
	next SPLIT;
      }
      #print STDERR "query coords ".$mapped[$i]->end." ".$mapped[$i]->start."\n";
      my $query_length = ($mapped[$i]->end - $mapped[$i]->start + 1);
      #print STDERR " query length = ".$query_length."\n";
      my $hit_length;
      if($self->_query_unit == $self->_hit_unit){
	$hit_length = $query_length;
      }elsif($self->_query_unit > $self->_hit_unit){
	my $tmp =  ($query_length/$self->_query_unit);
	#print STDERR "tmp = ".$tmp."\n";
	#$hit_length = int($tmp);
	$hit_length = sprintf "%.0f", $tmp;
      }elsif($self->_hit_unit > $self->_query_unit){
	my $tmp = ($query_length*$self->_hit_unit);
	#print STDERR "tmp = ".$tmp."\n";
	$hit_length = sprintf "%.0f", $tmp;
      }
      #print STDERR "hit length ".$hit_length."\n";
     
      my $hit_end = ($hit_start + $hit_length) - 1;
      #print "hit start ".$hit_start." hit end ".$hit_end."\n";
      my $rawContig = $rcAdaptor->fetch_by_dbID( $mapped[$i]->id() );
      my $f1 = new Bio::EnsEMBL::SeqFeature();
      my $f2 = new Bio::EnsEMBL::SeqFeature();
      my $new_feature = Bio::EnsEMBL::FeaturePair->new(-feature1=>$f1,
						       -feature2=>$f2);
     
      $new_feature->start($mapped[$i]->start);
      $new_feature->end($mapped[$i]->end);
      $new_feature->strand($mapped[$i]->strand);
      $new_feature->seqname($mapped[$i]->id);
      $new_feature->score($feature->score);
      $new_feature->percent_id($feature->percent_id);
      $new_feature->p_value($feature->p_value);
      $new_feature->hstart($hit_start);
      $new_feature->hend($hit_end);
      $new_feature->hstrand($feature->hstrand);
      $new_feature->hseqname($feature->hseqname);
      #$new_feature->hscore($feature->score);
      $new_feature->analysis($feature->analysis);
      $new_feature->attach_seq($rawContig);
      #print STDERR "split feature ".$new_feature->gffstring."\n";
      push(@out, $new_feature);
      $hit_start = ($hit_end + 1);
    }
  }else{
    if($mapped[0]->isa("Bio::EnsEMBL::Mapper::Gap")){
      $self->warn("piece of evidence lies on gap\n");
      return;
    }
    my $rawContig = $rcAdaptor->fetch_by_dbID( $mapped[0]->id() );
    my $f1 = new Bio::EnsEMBL::SeqFeature();
    my $f2 = new Bio::EnsEMBL::SeqFeature();
    my $new_feature = Bio::EnsEMBL::FeaturePair->new(-feature1=>$f1,
						     -feature2=>$f2);
    $new_feature->start($mapped[0]->start);
    $new_feature->end($mapped[0]->end);
    $new_feature->strand($mapped[0]->strand);
    $new_feature->seqname($mapped[0]->id);
    $new_feature->score($feature->score);
    $new_feature->percent_id($feature->percent_id);
    $new_feature->p_value($feature->p_value);
    $new_feature->hstart($feature->hstart);
    $new_feature->hend($feature->hend);
    $new_feature->hstrand($feature->hstrand);
    $new_feature->hseqname($feature->hseqname);
    #$new_feature->hscore($feature->score);
    $new_feature->analysis($feature->analysis);
    $new_feature->attach_seq($rawContig);

    push(@out, $new_feature);
  }

  foreach my $sf(@out){
    #print STDERR "gff ".$sf->gffstring."\n";
  }

  return @out;

}

1;
