#
# BioPerl module for Bio::EnsEMBL::SeqFeature
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::SeqFeature - Ensembl specific sequence feature.

=head1 SYNOPSIS

    my $feat = new Bio::EnsEMBL::SeqFeature(-seqname => 'pog',
                                            -start   => 100,
                                            -end     => 220,
                                            -strand  => -1,
                                            -frame   => 1,
                                            -source_tag  => 'tblastn_vert',
                                            -primary_tag => 'similarity',
                                            -analysis => $analysis
                                            );

    # $analysis is a Bio::EnsEMBL::Analysis object

    # SeqFeatureI methods can be used
    my $start = $feat->start;
    my $end   = $feat->end;

    # Bio::EnsEMBL::SeqFeature specific methods can be used
    my $analysis = $feat->analysis;

    # Validate all the data in the object
    $feat->validate  || $feat->throw("Invalid data in $feat");

=head1 DESCRIPTION

This is an implementation of the ensembl Bio::EnsEMBL::SeqFeatureI interface.
Extra
methods are to store details of the analysis program/database/version used
to create this data and also a method to validate all data in the object is
present and of the right type.  This is useful before writing into
a relational database for example.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::SeqFeature;

use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::SeqFeatureI;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Root Bio::EnsEMBL::SeqFeatureI);


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

      $self->_rearrange([qw(START
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

=head2 slide

  Arg [1]    : int $slide 
               the new feature start position
  Example    : $feature->move(100);
  Description: Shifts a features location.  This is faster
               then calling seperate accesors in a large loop.
  Returntype : none
  Exceptions : none
  Caller     : BaseFeatureAdaptor

=cut

sub slide {
  my $self = shift;
  my $slide = shift;

  $self->{'_gsf_start'} += $slide;
  $self->{'_gsf_end'} += $slide;
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
    #print STDERR "SeqFeature  Validate ".$self->strand."\n";
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
    print(STDERR "   Seqname     : [" . $self->{_gsf_seqname} . "]\n");
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

   my $seq = $self->entire_seq();

   if(defined $seqname) {
     $self->{_seqname} = $seqname;
   } else {
     if(defined $self->{_seqname}) {
       return $self->{_seqname};
     } elsif($seq && ref $seq && $seq->can('name')) {
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

 Title   : seq
 Usage   : $tseq = $sf->seq()
 Function: returns the truncated sequence (if there) for this
 Example :
 Returns :
 Args    :


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

sub sub_SeqFeature{
   my ($self) = @_;

   if($self->{'_gsf_sub_array'}){
     return @{$self->{'_gsf_sub_array'}};
   }else{
     return;
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

           If EXPAND is used, the parent's start/end/strand will
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
       if( !$self->contains($feat) ) {
           $self->throw("$feat is not contained within parent feature, and expansion is not valid");
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
   
   if ($self->strand == -1) {
     $strand = "-";
   }
   
   $str .= (defined $self->seqname)     ?   $self->seqname."\t"      :  "\t";
   $str .= (defined $self->source_tag)  ?   $self->source_tag."\t"   :  "\t";
   $str .= (defined $self->primary_tag) ?   $self->primary_tag."\t"  :  "\t";
   $str .= (defined $self->start)       ?   $self->start."\t"        :  "\t";
   $str .= (defined $self->end)         ?   $self->end."\t"          :  "\t";
   $str .= (defined $self->strand)      ?   $strand."\t"             :  ".\t";
   $str .= (defined $self->score)       ?   $self->score."\t"        :  "\t";
   $str .= (defined $self->phase)       ?   $self->phase."\t"        :  ".\t";

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


=head2 contig_id

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use Bio::EnsEMBL::SeqFeature::attach_seq or
               Bio::EnsEMBL::SeqFeature::entire_seq instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub contig_id{
   my ($self,$arg) = @_;

   $self->warn("Bio::EnsEMBL::SeqFeature::contig_id is deprecated. " .
	    "Use attach_seq instead to associate a contig with a SeqFeature");

   if($arg) {
     my $contig = $self->db->get_RawContigAdaptor->fetch_by_dbID($arg);
     $self->attach_seq($contig);
   }

   return $self->entire_seq();
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
    return $self->{'_gsf_seq'};
}


=head2 raw_seqname

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use Bio::EnsEMBL::SeqFeature::attach_seq->name or
               Bio::EnsEMBL::SeqFeature::raw_seqname instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub raw_seqname{
   my ($self,$arg) = @_;

   $self->warn("Bio::EnsEMBL::SeqFeature::raw_seqname is deprecated. " .
	       "Use Bio::EnsEMBL::entire_seq()->name() instead");

   my $seq = $self->entire_seq();

   if($arg && $seq && ref $seq && $seq->can('name')) {
     $self->entire_seq()->name($arg);
   }

   if($seq && ref $seq && $seq->can('name')) {
     return $seq->name;
   }

   return undef;
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
    elsif ($self->contig->isa( "Bio::EnsEMBL::Slice" )) {

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

  my $dbh = $self->contig->adaptor->db;

  my $mapper = $dbh->get_AssemblyMapperAdaptor->
   fetch_by_type($slice->assembly_type);
  my $rca = $dbh->get_RawContigAdaptor;

  my @mapped = $mapper->map_coordinates_to_assembly(
    $self->contig->dbID,
    $self->start,
    $self->end,
    $self->strand * $slice->strand
  );

  unless (@mapped) {
    $self->throw("couldn't map $self to Slice");
  }

  unless (@mapped == 1) {
    $self->throw("$self should only map to one chromosome - something bad has happened ...");
  }

  if ($mapped[0]->isa("Bio::EnsEMBL::Mapper::Gap")) {
    $self->warn("feature lies on gap\n");
    return;
  }

  # mapped coords are on chromosome - need to convert to slice
  my $global_start = $slice->chr_start;

  $self->start  ($mapped[0]->start - $global_start + 1);
  $self->end    ($mapped[0]->end   - $global_start + 1);
  $self->strand ($mapped[0]->strand);
  $self->seqname($mapped[0]->id);

  return $self;
}


sub _transform_between_Slices {
  my ($self, $slice) = @_;

  $self->throw("New contig [$slice] is not a Bio::EnsEMBL::Slice")
   unless $slice->isa("Bio::EnsEMBL::Slice");

  if ((my $c1 = $self->contig->chr_name) ne (my $c2 = $slice->chr_name)) {
    $self->warn("Can't transform between chromosomes: $c1 and $c2");
    return;
  }

  my $shift = $slice->chr_start - $self->contig->chr_start;

  $self->start ($self->start  - $shift);
  $self->end   ($self->end    - $shift);
  $self->strand($self->strand * $slice->strand);
  $self->contig($slice);

  return $self;
}


sub _transform_to_RawContig {
  my($self) = @_;

  $self->throw("can't transform coordinates of $self without a contig defined")
   unless $self->contig;

  my $dbh = $self->contig->adaptor->db;

  my $mapper = $dbh->get_AssemblyMapperAdaptor->
   fetch_by_type($self->contig->assembly_type);
  my $rca = $dbh->get_RawContigAdaptor;

  # need to convert to chromosomal coords before mapping
  my $global_start = $self->contig->chr_start;

  my @mapped = $mapper->map_coordinates_to_rawcontig(
    $self->contig->chr_name,
    $self->start  + $global_start - 1,
    $self->end    + $global_start - 1,
    $self->strand * $self->contig->strand
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
    $self->attach_seq($rca->fetch_by_dbID($mapped[0]->id));

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

        $self->start     ($coords[0]->start);
        $self->end       ($coords[0]->end);
        $self->strand    ($coords[0]->strand);
        $self->seqname   ($coords[0]->id);
        $self->attach_seq($rca->fetch_by_dbID($coords[0]->id));

	$self->warn("Feature [$self] truncated as lies partially on a gap");
        return $self;
    }

    unless ($self->is_splittable) {
        print STDERR "Feature spans >1 raw contig - can't split\n";
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

      $feat->start     ($map->start);
      $feat->end       ($map->end);
      $feat->strand    ($map->strand);
      $feat->attach_seq($rca->fetch_by_dbID($map->id));

      push @out, $feat;
    }
    return @out;
  }
}


1;
