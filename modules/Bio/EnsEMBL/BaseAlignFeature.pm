# EnsEMBL module for storing dna-protein pairwise alignments
#
# Cared for by Michele Clamp <michele@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself
#
=pod

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

  The cigar_string contains the ungapped pieces that make up the gapped alignment
  
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

  The alignment goes from 267 to 479 in sequence 1 and 7 to 75 in sequence 2 and the
  strand is -1.

  The alignment is made up of the following ungapped pieces :

  sequence 1 start 447 , sequence 2 start 7  , match length 33 , strand -1
  sequence 1 start 417 , sequence 2 start 18 , match length 27 , strand -1
  sequence 1 start 267 , sequence 2 start 27 , match length 137 , strand -1

  These ungapped pieces are made up into the following string (called a cigar string)
  "33M3I27M3I137M" with start 267 end 479 strand -1 hstart 7 hend 75 hstrand 1 and
  feature type would be DnaPepAlignFeature
  

=cut 


package Bio::EnsEMBL::BaseAlignFeature;

use Bio::EnsEMBL::FeaturePairI;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::FeaturePairI);

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
    
    my ($cigar_string,$features) = $self->_rearrange([qw(CIGAR_STRING FEATURES)],

                                                     @args);

    #print STDERR "@args\n";
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

  if( ! defined $slice ) {
    #Since slice arg is not defined -  we want raw contig coords
    if(( defined  $self->contig ) && 
       ( $self->contig->isa( "Bio::EnsEMBL::RawContig" )) ) {
      print STDERR "BaseAlignFeature::transform, you are already apparently in rawcontig coords so why try to transform to them\n";
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

  if($arg){
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


=head2 contig

  Arg [1]    : Bio::EnsEMBL::Slice or
               Bio::EnsEMBL::RawContig $contig
  Example    : none
  Description: The sequnence where the coordinates make sense. Channeled to
               attach_seq/entire_seq from BioPerl
  Returntype : Slice/RawContig
  Exceptions : none
  Caller     : general

=cut


sub contig {
  my $self = shift;
  if( @_ ) {
    my $value = shift;
    #print "setting exons contig to ".$value." \n";
    $self->attach_seq($value);
  }
  else {
    return $self->entire_seq();
  }
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

  if (ref($features) ne "ARRAY") {
    $self->throw("features must be an array reference not a [" . ref($features) . "]");
  }

  for my $f ( @$features ) {

    my $feats = scalar (@$features) - 1;
    my $strand     = $features->[0]->strand;
    my @f;

    if( $strand == 1 ) {
      
      @f = sort {$a->start <=> $b->start} @$features;
      
    } else {
      
      @f = sort { $b->start <=> $a->start} @$features;
    }
    
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
    
    my $f1start = $f[0]->start;
    my $f1end   = $f[$#f]->end;
    
    if ( $strand == 1 ) {
      $f1start = $f[0]->start;
      $f1end   = $f[$#f]->end;
    } else {
      $f1end   = $f[0]->end;
      $f1start = $f[$#f]->start;
    }
    
    my $f2start;
    my $f2end;
    
    if ( $hstrand == 1 ) {
      $f2start = $f[0]->hstart;
      $f2end   = $f[$#f]->hend;
    } else {
      $f2end = $f[0]->hend;
      $f2start = $f[$#f]->hstart;
    }
    
    foreach my $f (@f) {
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
      my $start2 = $f->hstart();
      
      if (defined($prev1)) {
	if ( $strand == 1 ) {
	  if ($f->start < $prev1) {
	    $self->throw("Inconsistent coordinates feature is forward strand hstart in current feature should be greater than hend in previous feature " . $f->start . " < ".$prev1."\n");
	  }
	} else {
	  if ($f->end > $prev1) {
	    $self->throw("Inconsistent coordinates in feature array feature is reverse strand hend should be less than previous hstart " . $f->end . " > $prev1");
	  }
	}
	
      }
      my $length = ($f->end - $f->start + 1);
      my $hlength = ($f->hend - $f->hstart + 1);

      # using multiplication to avoid rounding errors, hence the
      # switch from query to hit for the ratios

      if( $length * $hit_unit != $hlength * $query_unit ) {
	$self->throw( "Feature lengths not comparable Lengths:".$length." ".$hlength." Ratios:".$query_unit." ". $hit_unit );
      }

      my $hlengthfactor = 1;
      if( $query_unit == 1 && $hit_unit == 3 ) {
	$hlengthfactor = (1/3);
      }
      if( $query_unit == 3 && $hit_unit == 1 ) {
	$hlengthfactor = 3;
      }

      # find out the type of gap
      if( $strand == 1 ) {
	if( ( defined $prev1 ) && ( $f->start > $prev1 + 1  )) {
	  # I type gap
	  my $gap = $f->start - $prev1 - 1;
	  if( $gap == 1 ) {
	    $gap = "";
	  }
	  $string .= "$gap"."I";
	}
	$prev1 = $f->end();
      } else {
	if(( defined $prev1 ) && ($f->end + 1 < $prev1 )) {
	  # I type gap
	  my $gap = $prev1 - $f->end() - 1;
	  if( $gap == 1 ) {
	    $gap = "";
	  }
	  $string .= "$gap"."I";
	}
	$prev1 = $f->start();
      }
      
      if( $hstrand == 1 ) {
	if((  defined $prev2 ) && ( $f->hstart() > $prev2 + 1 )) {
	  # D type gap
	  my $gap = $f->hstart - $prev2 - 1;
	  my $gap2 = int( $gap * $hlengthfactor + 0.05 );
	  
	  if( $gap2 == 1 ) {
	    $gap2 = "";
	  }
	  $string .= "$gap2"."D";
	} 
	$prev2 = $f->hend();
      } else {
	if( ( defined $prev2 ) && ( $f->hend() + 1 < $prev2 )) {
	  # D type gap 
	  my $gap = $prev2 - $f->hend - 1;
	  my $gap2 = int( $gap * $hlengthfactor + 0.05 );
	  
	  if( $gap2 == 1 ) {
	    $gap2 = "";
	  }
	  $string .= "$gap2"."D";
	}
	$prev2 = $f->hstart();
      }
      
      my $matchlength = $f->end() - $f->start() + 1;
      if( $matchlength == 1 ) {
	$matchlength = "";
      }
      $string .= $matchlength."M";
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
    
  }
}



sub _transform_to_slice{
  my ($self, $slice ) = @_;
  $self->throw( "implented soon :-)" );
}



sub _transform_to_rawcontig {
  my ( $self ) = @_;
  
  if( !$self->entire_seq->isa('Bio::EnsEMBL::Slice') ) {
    $self->throw("Called transform on a feature without a slice as argument, and no slice in entire_seq slot");
  }

  my @unmapped = $self->ungapped_features;
  my $rcAdaptor = $self->adaptor()->db()->get_RawContigAdaptor();
  my $slice = $self->entire_seq();

  if( !defined $slice ) {
    my %mapped;
    
    # mapping from slice to raw contigs
    # Basic logic -
    #    - get ungapped features
    #    - map each one
    #    - throw if there is a gap
    #    - if not form an ungapped feature on the raw contig
    #    - after all ungapped features are mapped, making new featurepairs
    #      from the mapped ungapped features keyed by raw contig id

    $slice = $self->entire_seq();
    my $global_start = $slice->chr_start;
    my $global_end   = $slice->chr_end;
    
    my $mapper = $slice->adaptor->db->
      get_AssemblyMapperAdaptor->fetch_by_type($slice->assembly_type);

    foreach my $f ( @unmapped ) {
      my @mapped = $mapper->map_coordinates_to_rawcontig
	(
	 $slice->chr_name,
	 $f->start()+$global_start-1,
	 $f->end()+$global_start-1,
	 $f->strand()*$slice->strand()
	);
      
      if( ! @mapped ) {
	$self->throw( "feature could not map" );
      }

      my $hstart = $f->hstart;

      foreach my $co ( @mapped ) {
	if( $co->isa("Bio::EnsEMBL::Mapper::Gap") ) {
	  $self->throw("Feature transform mapped to gap - currently throwing an exception here");
	}
	      
	my $f1 = Bio::EnsEMBL::SeqFeature->new();
	$f1->start($co->start - $global_start +1 );
	$f1->end($co->end - $global_start +1);
	$f1->strand($co->strand * $slice->strand());
	$f1->seqname($co->id());
	$f1->attach_seq( $rcAdaptor->fetch_by_dbID( $co->id() ));

	my $f2 = Bio::EnsEMBL::SeqFeature->new();
	$f2->start($hstart );
	$f2->end($hstart+$co->end-$co->start+1);
	$f2->strand($f->hstrand);
	$f2->seqname($f->hseqname);
	
	$f1->score($f->score);
	$f2->score($f->score);

	my $sf = Bio::EnsEMBL::FeaturePair->new(-feature1 => $f1,
						-feature2 => $f2);
	$hstart = $f2->end+1;

	if( !defined $mapped{$co->id} ) {
	  $mapped{$co->id} = [];
	}

	push(@{$mapped{$co->id}},$sf);
      }

    }

    my @outputf;

    foreach my $rawc ( keys %mapped ) {

      my $outputf = $self->new( -features => \@{$mapped{$rawc}} );
      $outputf->analysis( $self->analysis() );
      $outputf->score( $self->score() );
      $outputf->percent_id( $self->percent_id() );
      $outputf->pvalue( $self->pvalue() );

      push(@outputf,$outputf);
    }
	
    return @outputf;
  }
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



1;
