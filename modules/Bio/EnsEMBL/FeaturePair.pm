
#
# BioPerl module for Bio::EnsEMBL::FeaturePair
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::FeaturePair - Stores sequence features which are
                            themselves hits to other sequence features.

=head1 SYNOPSIS

    my $feat  = new Bio::EnsEMBL::FeaturePair(-feature1 => $f1,
					      -feature2 => $f2,
					      );

    # Bio::SeqFeatureI methods can be used
    my $start = $feat->start;
    my $end   = $feat->end;

    # Bio::EnsEMBL::SeqFeatureI methods can be used
    my $analysis = $feat->analysis;
    
    $feat->validate  || $feat->throw("Invalid data in $feat");

    # Bio::FeaturePair methods can be used
    my $hstart = $feat->hstart;
    my $hend   = $feat->hend;

=head1 DESCRIPTION

A sequence feature object where the feature is itself a feature on another 
sequence - e.g. a blast hit where residues 1-40 of a  protein sequence 
SW:HBA_HUMAN has hit to bases 100 - 220 on a genomic sequence HS120G22.  
The genomic sequence coordinates are used to create one sequence feature 
$f1 and the protein coordinates are used to create feature $f2.  
A FeaturePair object can then be made

    my $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $f1,   # genomic
					   -feature2 => $f2,   # protein
					   );

This object can be used as a standard Bio::SeqFeatureI in which case

    my $gstart = $fp->start  # returns start coord on feature1 - genomic seq.
    my $gend   = $fp->end    # returns end coord on feature1.

In general standard Bio::SeqFeatureI method calls return information
in feature1.

Data in the feature 2 object are generally obtained using the standard
methods prefixed by h (for hit!)

    my $pstart = $fp->hstart # returns start coord on feature2 = protein seq.
    my $pend   = $fp->hend   # returns end coord on feature2.


If you wish to swap feature1 and feature2 around :

    $feat->invert

    $feat->start # etc. returns data in $feature2 object


=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::FeaturePair;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::SeqFeature;

@ISA = qw(Bio::EnsEMBL::SeqFeature);


sub new {
  my($class,@args) = @_;
  my $self = {};

  if( ref( $class) ) {
    $class = ref( $class );
  }

  bless ($self,$class);

  my ($feature1,$feature2) = 
      $self->_rearrange([qw(FEATURE1
			    FEATURE2
			    )],@args);

  if($feature1) {
    $self = $self->SUPER::new(-ANALYSIS => $feature1->analysis(),
		      -SEQNAME  => $feature1->seqname(),
		      -START    => $feature1->start(),
		      -END      => $feature1->end(),
		      -STRAND   => $feature1->strand(),
		      -FRAME    => $feature1->frame(),
		      -SCORE    => $feature1->score(),
		      -PERCENT_ID => $feature1->percent_id(),
		      -P_VALUE => $feature1->p_value(),
		      -PHASE => $feature1->phase(),
		      -END_PHASE => $feature1->end_phase());

    if($feature1->contig) {
      $self->contig($feature1->contig);
    }
  } else {
    $self = $self->SUPER::new( @args );
  }

  $feature2 && $self->feature2($feature2);

  # set stuff in self from @args
  return $self; # success - we hope!
}


=head2 feature1

 Title   : feature1
 Usage   : $f = $featpair->feature1
           $featpair->feature1($feature)
 Function: Get/set for the query feature
 Returns : Bio::SeqFeatureI
 Args    : none

=cut


sub feature1 {
  my ($self,$arg) = @_;

  if($arg) {
    #print STDERR "have ".$arg." for feature1 and it has ".$arg->entire_seq."\n";
      $self->start($arg->start());
      $self->end($arg->end());
      $self->strand($arg->strand());
      $self->frame($arg->frame());
      $self->score($arg->score());
      $self->seqname($arg->seqname());
      $self->percent_id($arg->percent_id());
      $self->p_value($arg->p_value());
      $self->phase($arg->phase());
      $self->end_phase($arg->end_phase());
      $self->analysis($arg->analysis);
      if($arg->entire_seq){
	$self->attach_seq($arg->entire_seq);
      }
    }

  return $self;
}

=head2 feature2

 Title   : feature2
 Usage   : $f = $featpair->feature2
           $featpair->feature2($feature)
 Function: Get/set for the hit feature
 Returns : Bio::SeqFeatureI
 Args    : none


=cut

sub feature2 {
    my ($self,$arg) = @_;
    #print "passing ".$arg." into feature2\n";
    if (defined($arg)) {
      unless(ref($arg) ne "" && $arg->isa("Bio::SeqFeatureI")) {
	$self->throw("Argument [$arg] must be a Bio::SeqFeatureI");
      }

      $self->hstart($arg->start());
      $self->hend($arg->end());
      $self->hstrand($arg->strand());
      $self->hseqname($arg->seqname());
      $self->hphase($arg->phase());
      $self->hend_phase($arg->end_phase());
     
      return $arg;
    } 
    
    my $seq = new Bio::EnsEMBL::SeqFeature(
		    -SEQNAME    => $self->{_hseqname},
		    -START      => $self->{_hstart},
		    -END        => $self->{_hend},
                    -STRAND     => $self->{_hstrand},
		    -SCORE      => $self->score(),
		    -PERCENT_ID => $self->percent_id(),
		    -P_VALUE    => $self->p_value(),
		    -PHASE      => $self->hphase,
		    -END_PHASE  => $self->hend_phase,
		    -ANALYSIS   => $self->analysis);

   
    return $seq;
}


=head2 hseqname

 Title   : hseqname
 Usage   : $featpair->hseqname($newval)
 Function: Get/set method for the name of
           feature2.
 Returns : value of $feature2->seqname
 Args    : newvalue (optional)


=cut

sub hseqname {
    my ($self,$arg) = @_;
    
    if (defined($arg)) {
     # print STDERR "hseqname being set ".$arg."\n";
      $self->{_hseqname} = $arg;
    }

    return $self->{_hseqname};
}


=head2 hstart

 Title   : hstart
 Usage   : $start = $featpair->hstart
           $featpair->hstart(20)
 Function: Get/set on the start coordinate of feature2
 Returns : integer
 Args    : none

=cut

sub hstart {
  my ($self,$value) = @_;
  
  if (defined($value)) {
    #print STDERR "htart being set to ".$value,"\n";
    $self->{_hstart} = $value;
  }
  
  return $self->{_hstart};
}

=head2 hend

 Title   : hend
 Usage   : $end = $featpair->hend
           $featpair->hend($end)
 Function: get/set on the end coordinate of feature2
 Returns : integer
 Args    : none

=cut

sub hend{
    my ($self,$value) = @_;

    if (defined($value)) {
      #print STDERR "hend being set to ".$value."\n"; 
      $self->{_hend} = $value;
    }

    return $self->{_hend};
}

=head2 hstrand

 Title   : hstrand
 Usage   : $strand = $feat->strand()
           $feat->strand($strand)
 Function: get/set on strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : none

=cut

sub hstrand{
    my ($self,$arg) = @_;

    if (defined($arg)) {
      #print "hstrand = ".$arg."\n";
      $self->{_hstrand} = $arg;
    } 
    
    return $self->{_hstrand};
}



=head2 invert

 Title   : invert
 Usage   : $tag = $feat->invert
 Function: Swaps feature1 and feature2 around
 Returns : Nothing
 Args    : none


=cut

sub invert {
    my ($self) = @_;

    my $tmp = $self->feature2;


    $self->feature2($self->feature1);
    $self->feature1($tmp);
}

=head2 validate

 Title   : validate
 Usage   : $sf->validate
 Function: Checks whether all data fields are filled
           in in the object and whether it is of
           the correct type.
           Throws an exception if it finds problems
 Example : $sf->validate
 Returns : nothing
 Args    : none


=cut

sub validate {
  my ($self) = @_;
  
  $self->SUPER::validate();
  $self->feature2->validate;
  #$f2->validate();
  
  # Now the analysis object
  if (defined($self->analysis)) {
    unless($self->analysis->isa("Bio::EnsEMBL::Analysis")) {
      $self->throw("Wrong type of analysis object");
    }
  } else {
    $self->throw("No analysis object defined");
  }
}

=head2 validate_prot_feature

 Title   : validate_prot_feature
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub validate_prot_feature{
   my ($self) = @_;

   $self->SUPER::validate_prot_feature(1);
   $self->feature2->validate_prot_feature(2);

   if (defined($self->analysis)) {
     unless($self->analysis->isa("Bio::EnsEMBL::Analysis")) {
       $self->throw("Wrong type of analysis object");
     }
   } else {
     $self->throw("No analysis object defined");
   }
 }

=head2 set_featurepair_fields

 Title   : set_featurepair_fields
 Usage   : $fp->set_featurepair_fields($start, $end, $strand,
           $score, $seqname, $hstart, $hend, $hstrand,
	   $hseqname, $analysis);
 Returns : nothing
 Args    : listed above, followed by optional $e_value, $perc_id, 
           $phase, $end_phase

=cut

sub set_featurepair_fields {
   my ($self, $start, $end, $strand, $score, $seqname, $hstart, $hend,
        $hstrand, $hseqname, $analysis, $e_value, $perc_id, 
        $phase, $end_phase) = @_;
   
   $self->throw('interface fault') if (@_ < 12 or @_ > 15);

   $self->start($start);
   $self->end($end);
   $self->strand($strand);
   $self->score($score);
   $self->seqname($seqname);
   $self->hstart($hstart);
   $self->hend($hend);
   $self->hstrand($hstrand);
   $self->hseqname($hseqname);
   $self->analysis($analysis);
   $self->p_value    ($e_value)   if (defined $e_value);
   $self->percent_id ($perc_id)   if (defined $perc_id);
   $self->phase      ($phase)     if (defined $phase);
   $self->end_phase  ($end_phase) if (defined $end_phase);
}



sub gffstring {
    my ($self) = @_;

    my $str = $self->SUPER::gffstring();

    my $hstrand = "+";
   
    if (($self->hstrand)&&($self->hstrand == -1)) {
      $hstrand = "-";
    }

    #Append a few FeaturePair specific things
    $str .= (defined $self->hseqname)   ?   $self->hseqname."\t"    :  "\t";
    $str .= (defined $self->hstart)     ?   $self->hstart."\t"      :  "\t";
    $str .= (defined $self->hend)       ?   $self->hend."\t"        :  "\t";
    $str .= (defined $self->hstrand)    ?   $hstrand."\t"           :  "\t";
    $str .= (defined $self->hphase)     ?   $self->hphase."\t"      :  ".\t";
    
    return $str;
}




=head2 hphase

 Title   : hphase
 Usage   : $hphase = $fp->hphase()
           $fp->hphase($hphase)
 Function: get/set on start hphase of predicted feature2
 Returns : [0,1,2]
 Args    : none if get, 0,1 or 2 if set. 

=cut

sub hphase {
  my ($self, $value) = @_;
  
  if (defined($value)) {
    $self->{_hphase} = $value;
  }
  
  return $self->{_hphase};
}


=head2 hend_phase

 Title   : hend_phase
 Usage   : $hend_phase = $feat->hend_phase()
           $feat->hend_phase($hend_phase)
 Function: get/set on end phase of predicted feature2
 Returns : [0,1,2]
 Args    : none if get, 0,1 or 2 if set. 

=cut

sub hend_phase {
  my ($self, $value) = @_;
    
  if (defined($value)) {
    $self->{_hend_phase} = $value;
  }
  
  return $self->{_hend_phase};
}



=head2 species

 Arg [1]    : string $genus_species_name (optional)
              e.g. Homo_sapiens or Mus_musculus
 Example    : 
 Description: get/set on the species of feature1
 Returntype : string
 Execeptions: none
 Caller     : general

=cut

sub species{
    my ($self,$arg) = @_;

    if (defined($arg)) {
        return $self->{'_species'} = $arg;
    }
    return $self->{'_species'};
}

=head2 hspecies

 Arg [1]    : string $genus_species_name (optional)
              e.g. Homo_sapiens or Mus_musculus
 Example    : 
 Description: get/set on the species of feature2
 Returntype : string
 Execeptions: none
 Caller     : general

=cut

sub hspecies{
    my ($self,$arg) = @_;

    if (defined($arg)) {
        return $self->{'_hspecies'} = $arg;
    } 
    return $self->{'_hspecies'};
}

1;
