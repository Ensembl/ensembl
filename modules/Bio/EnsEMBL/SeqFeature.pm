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

Bio::EnsEMBL::SeqFeature - Ensembl specific sequence feature.

=head1 DESCRIPTION

Do not use this module if you can avoid it. It has been replaced by
Bio::EnsEMBL::Feature.  This module has a long history of usage but has
become very bloated, and quite unweildy.  It was decided to replace
it completely with a smaller, light-weight feature class rather than
attempting to refactor this class, and maintain strict backwards
compatibility.

Part of the complexity of this class was in its extensive
inheritance. As an example the following is a simplified inheritance
heirarchy that was present for Bio::EnsEMBL::DnaAlignFeature:

  Bio::EnsEMBL::DnaAlignFeature
  Bio::EnsEMBL::BaseAlignFeature
  Bio::EnsEMBL::FeaturePair
  Bio::EnsEMBL::SeqFeature
  Bio::SeqFeatureI
  Bio::RangeI
  Bio::Root::RootI

The new Bio::EnsEMBL::Feature class is much shorter, and hopefully much
easier to understand and maintain.

=head1 METHODS

=cut


# Let the code begin...


package Bio::EnsEMBL::SeqFeature;

use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::Analysis;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);

sub new {
  my($caller,@args) = @_;

  my $self = {};

  if(ref $caller) {
    bless $self, ref $caller;
  } else {
    bless $self, $caller;
  }

  $self->{'_gsf_tag_hash'} = {};
  $self->{'_gsf_sub_array'} = [];
  $self->{'_parse_h'} = {};
  $self->{'_is_splittable'} = 0;

  my ($start,$end,$strand,$frame,$score,$analysis,$seqname, $source_tag,
      $primary_tag, $percent_id, $p_value, $phase, $end_phase) =

      &rearrange([qw(START
                            END
                            STRAND
                            FRAME
                            SCORE
                            ANALYSIS
                            SEQNAME
                            SOURCE_TAG
                            PRIMARY_TAG
                            PERCENT_ID
                            P_VALUE
                            PHASE
                            END_PHASE
                            )],@args);

  #  $gff_string && $self->_from_gff_string($gff_string);

  if ( defined $analysis  && $analysis ne "")   { $self->analysis($analysis)};
  if ( defined ($start) && $start ne "" )       { $self->start($start)};
  if ( defined ($end )  && $end   ne "" )       { $self->end($end)}
  if ( defined $strand  && $strand ne "")       { $self->strand($strand)}
  if ( defined $frame  && $frame ne "")         { $self->frame($frame)}
  if ( defined $score  && $score ne "")         { $self->score($score)}
  if ( defined $seqname && $seqname ne "")      { $self->seqname($seqname)};
  if ( defined $percent_id && $percent_id ne ""){ $self->percent_id($percent_id)};
  if ( defined $p_value && $p_value ne "")      { $self->p_value($p_value)};
  if ( defined $phase && $phase ne "")          { $self->phase($phase)};
  if ( defined $end_phase && $end_phase ne "")  { $self->end_phase($end_phase)};

  return $self; # success - we hope!

}






=head2 start

 Title   : start
 Usage   : $start = $feat->start
           $feat->start(20)
 Function: Get/set on the start coordinate of the feature
 Returns : integer
 Args    : none


=cut

sub start{
    my ($self,$value) = @_;

    if (defined($value)) {
        if ($value !~ /^\-?\d+/ ) {
        $self->throw("$value is not a valid start");
    }
	
    $self->{'_gsf_start'} = $value
   }

    return $self->{'_gsf_start'};

}

=head2 end

 Title   : end
 Usage   : $end = $feat->end
           $feat->end($end)
 Function: get/set on the end coordinate of the feature
 Returns : integer
 Args    : none


=cut

sub end{
    my ($self,$value) = @_;

    if (defined($value)) {
        if( $value !~ /^\-?\d+/ ) {
            $self->throw("[$value] is not a valid end");
        }
	
        $self->{'_gsf_end'} = $value;
    }

   return $self->{'_gsf_end'};
}

=head2 length

 Title   : length
 Usage   :
 Function:
 Example :
 Returns :
 Args    :


=cut

sub length{
   my ($self,@args) = @_;

   return $self->end - $self->start +1;
}


=head2 strand

 Title   : strand
 Usage   : $strand = $feat->strand()
           $feat->strand($strand)
 Function: get/set on strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : none


=cut

sub strand {
    my ($self,$value) = @_;

    if (defined($value)) {
        if( $value eq '+' ) { $value = 1; }
        if( $value eq '-' ) { $value = -1; }
        if( $value eq '.' ) { $value = 0; }

        if( $value != -1 && $value != 1 && $value != 0 ) {
            $self->throw("$value is not a valid strand info");
        }
        $self->{'_gsf_strand'} = $value;
    }

    return $self->{'_gsf_strand'};
}


=head2 move

  Arg [1]    : int $start
  Arg [2]    : int $end
  Arg [3]    : (optional) int $strand 
  Example    : $feature->move(100, 200, -1);
  Description: Moves a feature to a different location.  This is faster
               then calling 3 seperate accesors in a large loop.
  Returntype : none
  Exceptions : none
  Caller     : BaseFeatureAdaptor

=cut

sub move {
  my ($self, $start, $end, $strand)  = @_;

  $self->{'_gsf_start'} = $start;
  $self->{'_gsf_end'} = $end;
  if(defined $strand) {
    $self->{'_gsf_strand'} = $strand;
  }
}


=head2 score

 Title   : score
 Usage   : $score = $feat->score()
           $feat->score($score)
 Function: get/set on score information
 Returns : float
 Args    : none if get, the new value if set


=cut

sub score {
    my ($self,$value) = @_;

    if(defined ($value) ) {
      if( $value !~ /^[+-]?\d+\.?\d*(e-\d+)?/ ) {
          $self->throw("'$value' is not a valid score");
      }
      $self->{'_gsf_score'} = $value;
  }

  return $self->{'_gsf_score'};
}

=head2 frame

 Title   : frame
 Usage   : $frame = $feat->frame()
           $feat->frame($frame)
 Function: get/set on frame information
 Returns : 0,1,2
 Args    : none if get, the new value if set


=cut

sub frame {
    my ($self,$value) = @_;

    if (defined($value)) {
        if( $value != 1 && $value != 2 && $value != 3 ) {
            $self->throw("'$value' is not a valid frame");
       }
        $self->{'_gsf_frame'} = $value;
    }

    return $self->{'_gsf_frame'};
}

=head2 primary_tag

 Title   : primary_tag
 Usage   : $tag = $feat->primary_tag()
           $feat->primary_tag('exon')
 Function: get/set on the primary tag for a feature,
           eg 'exon'
 Returns : a string
 Args    : none


=cut

sub primary_tag{
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    # throw warnings about setting primary tag
    my ($p,$f,$l) = caller;
    $self->warn("$f:$l setting primary_tag now deprecated." .
		"Primary tag is delegated to analysis object");
  }

  unless($self->analysis) {
    return '';
  }
  
  return $self->analysis->gff_feature();
}

=head2 source_tag

 Title   : source_tag
 Usage   : $tag = $feat->source_tag()
           $feat->source_tag('genscan');
 Function: Returns the source tag for a feature,
           eg, 'genscan'
 Returns : a string
 Args    : none


=cut

sub source_tag{
    my ($self,$arg) = @_;

    if (defined($arg)) {
        # throw warnings about setting primary tag
        my ($p,$f,$l) = caller;
        $self->warn("$f:$l setting source_tag now deprecated. " .
		    "Source tag is delegated to analysis object");
    }

    unless($self->analysis) {
      return "";
    }

   return $self->analysis->gff_source();
}


=head2 analysis

 Title   : analysis
 Usage   : $sf->analysis();
 Function: Store details of the program/database
           and versions used to create this feature.

 Example :
 Returns :
 Args    :


=cut

sub analysis {
   my ($self,$value) = @_;

   if (defined($value)) {
     unless(ref($value) && $value->isa('Bio::EnsEMBL::Analysis')) {
       $self->throw("Analysis is not a Bio::EnsEMBL::Analysis object "
		    . "but a $value object");
     }

     $self->{_analysis} = $value;
   } else {
     #if _analysis is not defined, create a new analysis object
     unless(defined $self->{_analysis}) {
       $self->{_analysis} = new Bio::EnsEMBL::Analysis();
     }
   }

   return $self->{_analysis};
}

=head2 validate

 Title   : validate
 Usage   : $sf->validate;
 Function: Checks whether all the data is present in the
           object.
 Example :
 Returns :
 Args    :


=cut

sub validate {
    my ($self) = @_;

    $self->vthrow("Seqname not defined in feature")     unless defined($self->seqname);
    $self->vthrow("start not defined in feature")       unless defined($self->start);
    $self->vthrow("end not defined in feature")         unless defined($self->end);
    $self->vthrow("strand not defined in feature")      unless defined($self->strand);
    $self->vthrow("score not defined in feature")       unless defined($self->score);
    $self->vthrow("analysis not defined in feature")    unless defined($self->analysis);

    if ($self->end < $self->start) {
      $self->vthrow("End coordinate < start coordinate");
    }

}



sub vthrow {
    my ($self,$message) = @_;

    print(STDERR "Error validating feature [$message]\n");
    print(STDERR "   Seqname     : [" . $self->{_seqname} . "]\n");
    print(STDERR "   Start       : [" . $self->{_gsf_start} . "]\n");
    print(STDERR "   End         : [" . $self->{_gsf_end} . "]\n");
    print(STDERR "   Strand      : [" .
        ((defined ($self->{_gsf_strand})) ? $self->{_gsf_strand} : "undefined") . "]\n");

    print(STDERR "   Score       : [" . $self->{_gsf_score} . "]\n");

    print(STDERR "   Analysis    : [" . $self->{_analysis}->dbID . "]\n");

    $self->throw("Invalid feature - see dump on STDERR");
}


=head2 validate_prot_feature

 Title   : validate_prot_feature
 Usage   :
 Function:
 Example :
 Returns :
 Args    :


=cut

# Shouldn't this go as "validate" into Pro_SeqFeature?
sub validate_prot_feature{
    my ($self,$num) = @_;
    $self->throw("Seqname not defined in feature")     unless defined($self->seqname);
    $self->throw("start not defined in feature")       unless defined($self->start);
    $self->throw("end not defined in feature")         unless defined($self->end);
    if ($num == 1) {
        $self->throw("score not defined in feature")       unless defined($self->score);
        $self->throw("percent_id not defined in feature") unless defined($self->percent_id);
        $self->throw("evalue not defined in feature") unless defined($self->p_value);
    }
    $self->throw("analysis not defined in feature")    unless defined($self->analysis);
}



# These methods are specified in the SeqFeatureI interface but we don't want
# people to store data in them.  These are just here in order to keep
# existing code working


=head2 has_tag

 Title   : has_tag
 Usage   : $value = $self->has_tag('some_tag')
 Function: Returns the value of the tag (undef if
           none)
 Returns :
 Args    :


=cut

sub has_tag{
   my ($self,$tag) = (shift, shift);

   return exists $self->{'_gsf_tag_hash'}->{$tag};
}

=head2 add_tag_value

 Title   : add_tag_value
 Usage   : $self->add_tag_value('note',"this is a note");
 Returns : nothing
 Args    : tag (string) and value (any scalar)


=cut

sub add_tag_value{
   my ($self,$tag,$value) = @_;

   if( !defined $self->{'_gsf_tag_hash'}->{$tag} ) {
       $self->{'_gsf_tag_hash'}->{$tag} = [];
   }

   push(@{$self->{'_gsf_tag_hash'}->{$tag}},$value);
}

=head2 each_tag_value

 Title   : each_tag_value
 Usage   :
 Function:
 Example :
 Returns :
 Args    :


=cut

sub each_tag_value {
   my ($self,$tag) = @_;
   if( ! exists $self->{'_gsf_tag_hash'}->{$tag} ) {
       $self->throw("asking for tag value that does not exist $tag");
   }

   return @{$self->{'_gsf_tag_hash'}->{$tag}};
}


=head2 all_tags

 Title   : all_tags
 Usage   : @tags = $feat->all_tags()
 Function: gives all tags for this feature
 Returns : an array of strings
 Args    : none


=cut

sub all_tags{
   my ($self,@args) = @_;

   return keys %{$self->{'_gsf_tag_hash'}};
}



=head2 seqname

  Arg [1]    : string $seqname
  Example    : $seqname = $self->seqname();
  Description: Obtains the seqname of this features sequence.  This is set
               automatically when a sequence with a name is attached, or may
               be set manually.
  Returntype : string
  Exceptions : none
  Caller     : general, attach_seq

=cut

sub seqname{
   my ($self,$seqname) = @_;

   my $seq = $self->contig();

   if(defined $seqname) {
     $self->{_seqname} = $seqname;
   } else {
     if($seq && ref $seq && $seq->can('name')) {
       $self->{_seqname} = $seq->name();
     } 
   }

   return $self->{_seqname};
}


=head2 attach_seq

 Title   : attach_seq
 Usage   : $sf->attach_seq($seq)
 Function: Attaches a Bio::PrimarySeqI object to this feature. This
           Bio::PrimarySeqI object is for the *entire* sequence: ie
           from 1 to 10000
 Example :
 Returns :
 Args    :


=cut

sub attach_seq{
   my ($self, $seq) = @_;

   $self->contig($seq);
}

=head2 seq

 Example : $tseq = $sf->seq()
 Function: returns the sequence (if any ) for this feature truncated to the range spanning the feature
 Returns : a Bio::PrimarySeq object (I reckon)

=cut

sub seq{
   my ($self,$arg) = @_;

   if( defined $arg ) {
       $self->throw("Calling SeqFeature::Generic->seq with an argument. " .
		    "You probably want attach_seq");
   }

   if( ! exists $self->{'_gsf_seq'} ) {
       return undef;
   }

   # assumming our seq object is sensible, it should not have to yank
   # the entire sequence out here.

   my $seq = $self->{'_gsf_seq'}->trunc($self->start(),$self->end());


   if( $self->strand == -1 ) {

       # ok. this does not work well (?)
       #print STDERR "Before revcom", $seq->str, "\n";
       $seq = $seq->revcom;
       #print STDERR "After  revcom", $seq->str, "\n";
   }

   return $seq;
}

=head2 entire_seq

 Title   : entire_seq
 Usage   : $whole_seq = $sf->entire_seq()
 Function: gives the entire sequence that this seqfeature is attached to
 Example :
 Returns :
 Args    :


=cut

sub entire_seq{
   my ($self) = @_;

   return $self->contig;
}


=head2 sub_SeqFeature

 Title   : sub_SeqFeature
 Usage   : @feats = $feat->sub_SeqFeature();
 Function: Returns an array of sub Sequence Features
 Returns : An array
 Args    : none


=cut

sub sub_SeqFeature {
  my ($self) = @_;

  if ( $self->{'_gsf_sub_array'} ) {
    return @{ $self->{'_gsf_sub_array'} };
  } else {
    return ();
  }
}

=head2 add_sub_SeqFeature

 Title   : add_sub_SeqFeature
 Usage   : $feat->add_sub_SeqFeature($subfeat);
           $feat->add_sub_SeqFeature($subfeat,'EXPAND')
 Function: adds a SeqFeature into the subSeqFeature array.
           with no 'EXPAND' qualifer, subfeat will be tested
           as to whether it lies inside the parent, and throw
           an exception if not.

           If EXPAND is used, the parents start/end/strand will
           be adjusted so that it grows to accommodate the new
           subFeature
 Returns : nothing
 Args    : An object which has the SeqFeatureI interface


=cut

sub add_sub_SeqFeature{
   my ($self,$feat,$expand) = @_;

   if( !$feat->isa('Bio::SeqFeatureI') ) {
       $self->warn("$feat does not implement Bio::SeqFeatureI. Will add it anyway, but beware...");
   }

   if( $expand eq 'EXPAND' ) {
       # if this doesn't have start/end set - forget it!
       if( !defined $self->start && !defined $self->end ) {
           $self->start($feat->start());
           $self->end($feat->end());
           $self->strand($feat->strand);
       } else {
           my ($start,$end);
           if( $feat->start < $self->start ) {
               $start = $feat->start;
           }

           if( $feat->end > $self->end ) {
               $end = $feat->end;
           }

           $self->start($start);
           $self->end($end);

       }
   } else {
     if( !defined($feat->start()) || !defined($feat->end()) ||
         !defined($self->start())  || !defined($self->end())) {
       $self->throw( "This SeqFeature and the sub_SeqFeature must define".
                     " start and end.");
     }
     if($feat->start() > $feat->end() || $self->start() > $self->end()) {
       $self->throw("This SeqFeature and the sub_SeqFeature must have " .
                    "start that is less than or equal to end.");
     }
     if($feat->start() < $self->start() || $feat->end() > $self->end() ) {
       $self->throw("$feat is not contained within parent feature, " .
                    "and expansion is not valid");
     }
   }

   push(@{$self->{'_gsf_sub_array'}},$feat);

}

=head2 flush_sub_SeqFeature

 Title   : flush_sub_SeqFeature
 Usage   : $sf->flush_sub_SeqFeature
 Function: Removes all sub SeqFeature
           (if you want to remove only a subset, take
            an array of them all, flush them, and add
            back only the guys you want)
 Example :
 Returns : none
 Args    : none


=cut

sub flush_sub_SeqFeature {
   my ($self) = @_;

   $self->{'_gsf_sub_array'} = []; # zap the array implicitly.
}


sub id {
    my ($self,$value) = @_;

    if (defined($value)) {
        $self->{_id} = $value;
    }

    return $self->{_id};

}

=head2 percent_id

 Title   : percent_id
 Usage   : $pid = $feat->percent_id()
           $feat->percent_id($pid)
 Function: get/set on percentage identity information
 Returns : float
 Args    : none if get, the new value if set

=cut

sub percent_id {
    my ($self,$value) = @_;

    if (defined($value))
    {
            $self->{_percent_id} = $value;
    }

    return $self->{_percent_id};
}

=head2 p_value

 Title   : p_value
 Usage   : $p_val = $feat->p_value()
           $feat->p_value($p_val)
 Function: get/set on p value information
 Returns : float
 Args    : none if get, the new value if set

=cut

sub p_value {
    my ($self,$value) = @_;

    if (defined($value))
    {
            $self->{_p_value} = $value;
    }

    return $self->{_p_value};
}

=head2 phase

 Title   : phase
 Usage   : $phase = $feat->phase()
           $feat->phase($phase)
 Function: get/set on start phase of predicted exon feature
 Returns : [0,1,2]
 Args    : none if get, 0,1 or 2 if set.

=cut

sub phase {
    my ($self, $value) = @_;

    if (defined($value) )
    {
        $self->throw("Valid values for Phase are [0,1,2]") if ($value < 0 || $value > 2);
            $self->{_phase} = $value;
    }

    return $self->{_phase};
}

=head2 end_phase

 Title   : end_phase
 Usage   : $end_phase = $feat->end_phase()
           $feat->end_phase($end_phase)
 Function: returns end_phase based on phase and length of feature
 Returns : [0,1,2]
 Args    : none if get, 0,1 or 2 if set.

=cut

sub end_phase {
   my ($self, $value) = @_;

    if (defined($value))
    {
            $self->throw("Valid values for Phase are [0,1,2]") if ($value < 0 || $value > 2);
            $self->{_end_phase} = $value;
    }

    return $self->{_end_phase};
}


sub gffstring {
   my ($self) = @_;

   my $str;

   my $strand = "+";
   
   if ((defined $self->strand)&&($self->strand == -1)) {
     $strand = "-";
   }
   
   $str .= (defined $self->seqname)     ?   $self->seqname."\t"      :  "\t";
   $str .= (defined $self->source_tag)  ?   $self->source_tag."\t"   :  "\t";
   $str .= (defined $self->primary_tag) ?   $self->primary_tag."\t"  :  "\t";
   $str .= (defined $self->start)       ?   $self->start."\t"        :  "\t";
   $str .= (defined $self->end)         ?   $self->end."\t"          :  "\t";
   $str .= (defined $self->score)       ?   $self->score."\t"        :  "\t";
   $str .= (defined $self->strand)      ?   $strand."\t"             :  ".\t";
   $str .= (defined $self->phase)       ?   $self->phase."\t"        :  ".\t";
   eval{
     $str .= (defined $self->end_phase)   ?   $self->end_phase."\t"        :  ".\t";
   };

   return $str;
}


=head2 external_db

 Title   : external_db
 Usage   : $pid = $feat->external_db()
           $feat->external_db($dbid)
 Function: get/set for an external db accession number (e.g.: Interpro)
 Returns :
 Args    : none if get, the new value if set

=cut

sub external_db {
    my ($self,$value) = @_;

    if (defined($value))
    {
            $self->{'_external_db'} = $value;
    }

    return $self->{'_external_db'};
}



=head2 contig

  Arg [1]    : Bio::PrimarySeqI $seq
  Example    : $seq = $self->contig;
  Description: Accessor to attach/retrieve a sequence to/from a feature
  Returntype : Bio::PrimarySeqI
  Exceptions : none
  Caller     : general

=cut

sub contig {
  my ($self, $arg) = @_;

  if ($arg) {
    unless (defined $arg && ref $arg && $arg->isa("Bio::PrimarySeqI")) {
      $self->throw("Must attach Bio::PrimarySeqI objects to SeqFeatures");
    }
    
    $self->{'_gsf_seq'} = $arg;
    
    # attach to sub features if they want it
    
    foreach my $sf ($self->sub_SeqFeature) {
      if ($sf->can("attach_seq")) {
	$sf->attach_seq($arg);
      }
    }
  }
  #print STDERR "contig is ".$self->{'_gsf_seq'}." with name ".$self->{'_gsf_seq'}->name."\n" unless(!$self->{'_gsf_seq'});
#  my ($p, $f, $l) = caller;
#  print STDERR "Caller = ".$f." ".$l."\n";
  return $self->{'_gsf_seq'};
}



sub is_splittable {
   my ($self, $arg) = @_;

   if (defined $arg) {
       $self->{'_is_splittable'} = $arg;
   }
   return $self->{'_is_splittable'};
}


sub transform {
  my ($self, $slice) = @_;

  unless (defined $slice) {

    if ((defined $self->contig) &&
      ($self->contig->isa("Bio::EnsEMBL::RawContig"))) {

      # we are already in rawcontig coords, nothing needs to be done
      return $self;

    }
    else {
      # transform to raw_contig coords from Slice coords
      return $self->_transform_to_RawContig();
    }
  }

  if (defined $self->contig) {

    if ($self->contig->isa("Bio::EnsEMBL::RawContig"))  {

      # transform to slice coords from raw contig coords
      return $self->_transform_to_Slice($slice);
    }
    elsif ($self->contig->isa( "Bio::EnsEMBL::Slice" ) or $self->contig->isa( "Bio::EnsEMBL::LRGSlice" )) {

      # transform to slice coords from other slice coords
      return $self->_transform_between_Slices($slice);
    }
    else {

      # Unknown contig type
      $self->throw("Cannot transform unknown contig type @{[$self->contig]}");
    }
  }
  else {

    #Can't convert to slice coords without a contig to work with
    return $self->throw("Object's contig is not defined - cannot transform");
  }

}


sub _transform_to_Slice {
  my ($self, $slice) = @_;

  $self->throw("can't transform coordinates of $self without a contig defined")
    unless $self->contig;

  unless($self->contig->adaptor) {
    $self->throw("cannot transform coordinates of $self without adaptor " .
		 "attached to contig");
  }

  my $dbh = $self->contig->adaptor->db;

  my $mapper = 
    $dbh->get_AssemblyMapperAdaptor->fetch_by_type($slice->assembly_type);
  my $rca = $dbh->get_RawContigAdaptor;

  my @mapped = $mapper->map_coordinates_to_assembly(
    $self->contig->dbID,
    $self->start,
    $self->end,
    $self->strand
  );

  unless (@mapped) {
    $self->throw("couldn't map $self to Slice");
  }

  unless (@mapped == 1) {
    $self->throw("$self should only map to one chromosome - " .
		 "something bad has happened ...");
  }

  if ($mapped[0]->isa("Bio::EnsEMBL::Mapper::Gap")) {
    $self->warn("feature lies on gap\n");
    return;
  }

  if( ! defined $slice->chr_name() ) {
    my $slice_adaptor = $slice->adaptor();
    %$slice = %{$slice_adaptor->fetch_by_chr_name( $mapped[0]->id() )};
  }
   
  # mapped coords are on chromosome - need to convert to slice
  if($slice->strand == 1) {
    $self->start  ($mapped[0]->start - $slice->chr_start + 1);
    $self->end    ($mapped[0]->end   - $slice->chr_start + 1);
    $self->strand ($mapped[0]->strand);
  } else {
    $self->start  ($slice->chr_end - $mapped[0]->end   + 1);
    $self->end    ($slice->chr_end - $mapped[0]->start + 1);
    $self->strand ($mapped[0]->strand * -1);
  }

  $self->seqname($mapped[0]->id);

  #set the contig to the slice
  $self->contig($slice);

  return $self;
}


sub _transform_between_Slices {
  my ($self, $to_slice) = @_;

  my $from_slice = $self->contig;

  $self->throw("New contig [$to_slice] is not a Bio::EnsEMBL::Slice")
   unless ($to_slice->isa("Bio::EnsEMBL::Slice") or $to_slice->isa("Bio::EnsEMBL::LRGSlice") );

  if ((my $c1 = $from_slice->chr_name) ne (my $c2 = $to_slice->chr_name)) {
    $self->warn("Can't transform between chromosomes: $c1 and $c2");
    return;
  }

  my($start, $end, $strand);

  #first convert to assembly coords
  if($from_slice->strand == 1) {
    $start  = $from_slice->chr_start + $self->start - 1;
    $end    = $from_slice->chr_start + $self->end   - 1;
    $strand = $self->strand;
  } else {
    $start  = $from_slice->chr_end - $self->end   + 1;
    $end    = $from_slice->chr_end - $self->start + 1;
    $strand = $self->strand;
  }

  #now convert to the other slice's coords 
  if($to_slice->strand == 1) {
    $self->start ($start - $to_slice->chr_start + 1); 
    $self->end   ($end   - $to_slice->chr_start + 1); 
    $self->strand($strand);
  } else {
    $self->start ($to_slice->chr_end - $end   + 1);
    $self->end   ($to_slice->chr_end - $start + 1);
    $self->strand($strand * -1);
  }

  $self->contig($to_slice);

  return $self;
}


sub _transform_to_RawContig {
  my($self) = @_;

  #print STDERR "transforming ".$self." to raw contig coords\n";
  $self->throw("can't transform coordinates of $self without a contig defined")
   unless $self->contig;

  my $slice = $self->contig;

  unless($slice->adaptor) {
    $self->throw("can't transform coordinates of $self without an adaptor " .
		 "attached to the feature's slice");
  }

  my $dbh = $slice->adaptor->db;

  my $mapper = 
    $dbh->get_AssemblyMapperAdaptor->fetch_by_type($slice->assembly_type);
  my $rca = $dbh->get_RawContigAdaptor;

  #first convert the features coordinates to assembly coordinates
  my($start, $end, $strand);
  if($slice->strand == 1) {
    $start  = $slice->chr_start + $self->start - 1;
    $end    = $slice->chr_start + $self->end   - 1;
    $strand = $self->strand;
  } else {
    $start = $slice->chr_end - $self->end   + 1;
    $end   = $slice->chr_end - $self->start + 1;
    $strand = $self->strand * -1;
  }

  #convert the assembly coordinates to RawContig coordinates
  my @mapped = $mapper->map_coordinates_to_rawcontig(
    $slice->chr_name,
    $start,
    $end,
    $strand
  );

  unless (@mapped) {
    $self->throw("couldn't map $self");
    return $self;
  }

  if (@mapped == 1) {

    if ($mapped[0]->isa("Bio::EnsEMBL::Mapper::Gap")) {
      $self->warn("feature lies on gap\n");
      return;
    }

    my $rc = $rca->fetch_by_dbID($mapped[0]->id);

    $self->start     ($mapped[0]->start);
    $self->end       ($mapped[0]->end);
    $self->strand    ($mapped[0]->strand);
    $self->seqname   ($mapped[0]->id);
    #print STDERR "setting contig to be ".$mapped[0]->id."\n";
    $self->contig($rca->fetch_by_dbID($mapped[0]->id));

    return $self;
  }
  else {

    # more than one object returned from mapper
    # possibly more than one RawContig in region

    my (@gaps, @coords);

    foreach my $m (@mapped) {

	if ($m->isa("Bio::EnsEMBL::Mapper::Gap")) {
	    push @gaps, $m;
	}
	elsif ($m->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
	    push @coords, $m;
	}
    }

    # case where only one RawContig maps
    if (@coords == 1) {

        $self->start  ($coords[0]->start);
        $self->end    ($coords[0]->end);
        $self->strand ($coords[0]->strand);
        $self->seqname($coords[0]->id);
	#print STDERR "2 setting contig to be ".$coords[0]->id."\n";
        $self->contig ($rca->fetch_by_dbID($coords[0]->id));

	$self->warn("Feature [$self] truncated as lies partially on a gap");
        return $self;
    }

    unless ($self->is_splittable) {
        $self->warn("Feature spans >1 raw contig - can't split\n");
        return;
    }

    my @out;
    my $obj = ref $self;

    SPLIT: foreach my $map (@mapped) {

      if ($map->isa("Bio::EnsEMBL::Mapper::Gap")) {
	$self->warn("piece of evidence lies on gap\n");
	next SPLIT;
      }

      my $feat = $obj->new;

      $feat->start  ($map->start);
      $feat->end    ($map->end);
      $feat->strand ($map->strand);
      #print STDERR "3 setting contig to be ".$mapped[0]->id."\n";
      $feat->contig ($rca->fetch_by_dbID($map->id));
      $feat->adaptor($self->adaptor) if $self->adaptor();
      $feat->display_label($self->display_label) if($self->can('display_label'));
      $feat->analysis($self->analysis);
      push @out, $feat;
    }
    
    return @out;
  }
}


1;
