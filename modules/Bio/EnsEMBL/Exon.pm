#
# BioPerl module for Exon
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Exon - Confirmed Exon 

=head1 SYNOPSIS

    $ex = new Bio::EnsEMBL::Exon;

    $ex->start(10);
    $ex->end(100);

Examples of creating an exon

    $ex = new Bio::EnsEMBL::Exon(1208,1506,1) # start = 1208, end = 1506, forward strand
    
    Start and end coordinates are always stored with start < end.  If they are input
    in the reverse order they will be swapped over.  The value for the strand will
    be kept as its input value;

    Strand values:  + or  1 = forward strand
                    - or -1 = reverse strand
                    . or  0 = unknown strand

    $ex->dna_seq($dna);                       # $dna is a Bio::Seq
    $ex->phase(0);                            # Sets the phase of the exon
    $ex->end_phase();                         # Calculates the end_phase of the exon from the
                                              # Length of the dna and the starting phase

    Phase values  are 0,1,2

    $trans = $ex->translate();                # Translates the exon. Returns a Bio::Seq

Frameshifts - these haven''t been coded yet :-)

    Frameshifts in the exon are stored as an array of positions and an array of lengths

    my @pos =  $ex->frameshift_position()   # Array of start coordinates of frameshifts
    my @len =  $ex->frameshift_length()     # Array of lengths of frameshifts

    Setting frameshifts
    
    $ex->add_frameshift(1208,1)             # Adds a frameshift at position 1208 and it
                                            # is an insertion of 1 base
    $ex->add_frameshift(6509,-2)            # Adds a frameshift at position 6509 and it
                                            # is a deletion of 2 bases

    $ex->translate()  translates the dna sequence in the correct phase taking the frameshifts into account.

=head1 DESCRIPTION

Exon object.  

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Exon;
use vars qw(@ISA $AUTOLOAD);
use strict;

# Object preamble - inherits from Bio::SeqFeature::Generic

use Bio::SeqFeature::Generic;
use Bio::Seq; # exons have to have sequences...

@ISA = qw(Bio::SeqFeature::Generic Exporter);


sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  # Array to store supporting evidence for this exon
  $self->{_supporting_evidence} = [];

  # add in EnsEMBL tag as 1.

  $self->primary_tag('exon');
  $self->source_tag('EnsEMBL');

  # Parse the input paramters (start,end,strand)
  if ($#args == 2) {
      $self->_parse_args(@args);
  }
  # set exon rank to be 1 be default
  $self->sticky_rank(1);

  # set stuff in self from @args
  return $self; # success - we hope!
}

# Parse routine called from the constructor to
# set the basic variables start,end and strand

sub _parse_args {

  my ($self,$start,$end,$strand) = @_;

  # Swap start and end if they're in the wrong order
  # We assume that the strand is correct and keep the input value.

  if ($start > $end) {
    my $tmp = $end;
    $end    = $start;
    $start  = $tmp;
  }

  
  $self->start ($start);
  $self->end   ($end);
  $self->strand($strand);

}

=pod 

=head1 Methods unique to exon


=pod 

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Example : 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id{
   my ($self) = shift;
   if( @_ ) {
       $self->{'id'} = shift;
   }
   return $self->{'id'};
   
}

=head2 version

 Title   : version
 Usage   : $obj->version($newval)
 Function: 
 Returns : value of version
 Args    : newvalue (optional)


=cut

sub version{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'version'} = $value;
    }
    return $obj->{'version'};

}

=head2 contig_id

 Title   : contig_id
 Usage   : $obj->contig_id($newval)
 Function: 
 Returns : value of contig_id
 Args    : newvalue (optional)


=cut

sub contig_id{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'contigid'} = $value;
    }
    return $self->{'contigid'};

}

=head2 clone_id

 Title   : clone_id
 Usage   : $obj->clone_id($newval)
 Function: 
 Returns : value of clone_id
 Args    : newvalue (optional)


=cut

sub clone_id{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'clone_id'} = $value;
    }
    return $obj->{'clone_id'};

}

=head2 created

 Title   : created
 Usage   : $obj->created($newval)
 Function: 
 Returns : value of created
 Args    : newvalue (optional)


=cut

sub created{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'created'} = $value;
    }
    return $self->{'created'};

}

=head2 modified

 Title   : modified
 Usage   : $obj->modified($newval)
 Function: 
 Returns : value of modified
 Args    : newvalue (optional)


=cut

sub modified{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'modified'} = $value;
    }
    return $self->{'modified'};

}

=head2 _stored

 Title   : _stored
 Usage   : $obj->_stored($newval)
 Function: 
 Returns : value of stored
 Args    : newvalue (optional)


=cut

sub _stored{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'stored'} = $value;
    }
    return $self->{'stored'};

}

=head2 pep_seq

  Title   : pep_seq
  Usage   : @pep = $feat->pep_seq
  Function: Returns the 3 frame peptide translations of the exon
  Returns : @Bio::Seq
  Args    : none

=cut

sub pep_seq {
    my ($self) = @_;
    my $pep = $self->_translate();
    return @$pep;
}

=head2 sticky_rank

 Title   : sticky_rank
 Usage   : $obj->sticky_rank($newval)
 Function: 
 Returns : value of sticky_rank
 Args    : newvalue (optional)


=cut

sub sticky_rank{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'sticky_rank'} = $value;
    }
    return $obj->{'sticky_rank'};

}





=head2 end_phase

  Title   : end_phase
  Usage   : $end_phase = $feat->end_phase
  Function: Returns the end phase of the exon
  Returns : int
  Args    : none

=cut

sub end_phase {
    my ($self) = @_;

    defined($self->phase()) || $self->throw("Can't return end_phase if phase is not set");
    defined($self->start()) || $self->throw("Can't return end_phase if start coordinate is not set");
    defined($self->end())   || $self->throw("Can't return end_phase if end coordinate is not set");

    my $len   = $self->end() - $self->start() + 1;
    my $phase = $self->phase();
    my( $end_phase );
    if ($phase == -1) {
        $end_phase = -1;
    } else {
        $end_phase = ($len + $phase) % 3;
    }
    
    return $end_phase;
}

=pod

=head2 phase

  Title   : phase
  Usage   : $phase = $feat->phase
  Function: Returns the phase of the exon
  Returns : int
  Args    : none

=cut

sub phase {
  my ($self,$value) = @_;
  
  if (defined($value)) {
    # Value must be 0,1,2, or -1 for non-coding
    if ($value =~ /^(-1|0|1|2)$/) {
#	print STDERR "Setting phase for " . $self->id . " to $value\n";
      $self->{'phase'} = $value;
    } else {
      $self->throw("Bad value ($value) for exon phase. Should only be -1,0,1,2\n");
    }
  }
  return $self->{'phase'};
}

=pod

=head2 frame

  Title   : frame
  Usage   : $frame = $feat->frame
  Function: Returns the frame of the exon in the parent DNA
  Returns : int
  Args    : none

=cut

sub frame {
  my ($self,$value) = @_;
  
  if (defined($value)) {
    # Value must be 0,1,2,
    if ($value < 0 || $value > 2) {
      $self->throw("Bad value ($value) for exon frame. Should only be 0,1,2\n");
    } else {
      $self->{'frame'} = $value;
    }
  }
  return $self->{'frame'};
}

=pod

=head2 type

  Title   : type
  Usage   : $type = $feat->type
  Function: Returns the type of the exon (Init,Intr,Term)
  Returns : String
  Args    : none

=cut

sub type {
  my ($self,$value) = @_;
  
  if (defined($value)) {
    $self->{'type'} = $value;
  }
  return $self->{'type'};
}


=pod

=head2 _translate

  Title   : _translate
  Usage   : @pep = $feat->_translate()
  Function: Returns the 3 frame translations of the exon
  Returns : @Bio::Seq
  Args    : none

=cut

sub _translate {
    my($self) = @_;
    my $pep;
    my $i;
  
    # changed this to work with the new SeqFeature stuff. I am still not
    # 100% happy about this. EB.
  
    # Get the DNA sequence and create the sequence string
    $self->seq() || $self->throw("No DNA in object. Can't translate\n");

    my $dna = $self->seq()->seq();
  
    # Translate in all frames - have to chop
    # off bases from the beginning of the dna sequence 
    # for frames 1 and 2 as the translate() method
    # only translates in one frame. Pah!
  
    for ($i = 0; $i < 3; $i++) {
	my $tmp = new Bio::Seq(-seq => substr($dna,$i));
	$pep->[$i] = $tmp->translate();
    }
    return $pep;
}

=pod

=head2 translate

  Title   : translate
  Usage   : $pep = $feat->translate()
  Function: Returns the translation of the exon in the defined phase
  Returns : Bio::Seq
  Args    : none

=cut

sub translate {
    my($self) = @_;
    my $pep = $self->_translate() || throw ("Can't translate DNA\n");
    my $phase=$self->phase();
    if($phase){$phase=3-$phase;}
    return $pep->[$phase];
}

=head2 strand

 Title   : strand
 Usage   : $strand = $feat->strand()
 Function: Returns strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : + or 1, - or -1, . or 0

=cut

sub strand {
  my ($self,$value) = @_;

  
   if( defined $value)  {
      if ($value eq "+") {$value =  1;}
      if ($value eq "-") {$value = -1;}
      if ($value eq ".") {$value =  0;}

      $self->{'strand'} = $value;
    }
    return $self->{'strand'};
}

=head2 start_translation

 Title   : start_translation
 Usage   : $start_translation = $feat->start_translation
 Function: Returns coordinate taking into account phase
 Returns : number
 Args    : none

=cut

sub start_translation {
    my ($self,$value) = @_;
  
    if( defined $value){
	$self->throw("cannot set translation start!");
    }

    my $phase=$self->phase;
    $phase=3-$phase if $phase;
    if($self->strand==1){
	return ($self->start + $phase);
    }else{
	return ($self->end - $phase);
    }
}

=head2 end_translation

 Title   : end_translation
 Usage   : $end_translation = $feat->end_translation
 Function: Returns coordinate taking into account phase
 Returns : number
 Args    : none

=cut

sub end_translation {
    my ($self,$value) = @_;
  
    if( defined $value){
	$self->throw("cannot set translation end!");
    }

    my $phase=$self->end_phase;
    if($self->strand==1){
	return ($self->end - $phase);
    }else{
	return ($self->start + $phase);
    }
}

=head2 _rephase_exon_genscan

 Title   : _rephase_exon_genscan
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _rephase_exon_genscan{
   my ($self) = @_;

   my $dna = $self->seq();
   my $phase;
   my $pep = $self->_genscan_peptide();

   # trim top and bottem
   my $dnaseq = $dna->seq();

   LOOP : {

       $dna->setseq(substr($dnaseq,3,-3));
       
       my $tr1 = $dna->translate();
       my $tr1pep = $tr1->str();
       chop $tr1pep;

#        print STDERR "[$tr1pep] to\n[$pep]\n";
       if( $tr1pep !~ /\*/ && $pep =~ /$tr1pep/ ) {
	   $phase = 0;
	   last LOOP;
       }
       
       $dna->setseq(substr($dnaseq,4,-3));
       
       $tr1 = $dna->translate();
       $tr1pep = $tr1->str();
       chop $tr1pep;

#       print STDERR "[$tr1pep] to\n[$pep]\n";
       if( $tr1pep !~ /\*/ && $pep =~ /$tr1pep/ ) {
	   $phase = 1;
	   last LOOP;
       }
       
       $dna->setseq(substr($dnaseq,5,-3));
       
       $tr1 = $dna->translate();
       $tr1pep = $tr1->str();
       chop $tr1pep;

#       print STDERR "[$tr1pep] to\n[$pep]\n";
       if( $tr1pep !~ /\*/ && $pep =~ /$tr1pep/ ) {
	   $phase = 2;
	   last LOOP;
       }
       
   }

   

#   print STDERR "For exon ",$self->strand," ",$self->phase," ",$phase,"\n";
   $self->phase($phase);
   

}

=head2 _genscan_peptide

 Title   : _genscan_peptide
 Usage   : $obj->_genscan_peptide($newval)
 Function: 
 Returns : value of _genscan_peptide
 Args    : newvalue (optional)


=cut

sub _genscan_peptide{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_genscan_peptide'} = $value;
    }
    return $obj->{'_genscan_peptide'};

}


=head2 add_Supporting_Feature

 Title   : add_Supporting_Feature
 Usage   : $obj->add_Supporting_Feature($feature)
 Function: 
 Returns : Nothing
 Args    : Bio::EnsEMBL::SeqFeature


=cut


sub add_Supporting_Feature {
    my ($self,$feature) = @_;

    $self->throw("Supporting evidence [$feature] not Bio::EnsEMBL::SeqFeatureI") unless 
	defined($feature) &&  $feature->isa("Bio::EnsEMBL::SeqFeatureI");

    $self->{_supporting_evidence} = [] unless defined($self->{_supporting_evidence});

    push(@{$self->{_supporting_evidence}},$feature);
}



=head2 each_Supporting_Feature

 Title   : each_Supporting_Feature
 Usage   : my @f = $obj->each_Supporting_Feature
 Function: 
 Returns : @Bio::EnsEMBL::Feature
 Args    : none


=cut


sub each_Supporting_Feature {
    my ($self) = @_;

    $self->{_supporting_evidence} = [] unless defined($self->{_supporting_evidence});

    return @{$self->{_supporting_evidence}};

}

=head2 find_supporting_evidence

 Title   : find_supporting_evidence
 Usage   : $obj->find_supporting_evidence(\@features)
 Function: Looks through all the similarity features and
           stores as supporting evidence any feature
           that overlaps with an exon.  I know it is
           a little crude but it\'s a start/
 Example : 
 Returns : Nothing
 Args    : Bio::EnsEMBL::Exon


=cut


sub find_supporting_evidence {
    my ($self,$features,$sorted) = @_;

    FEAT : foreach my $f (@$features) {
	# return if we have a sorted feature array
	if ($sorted == 1 && $f->start > $self->end) {
	    return;
	}
	if ($f->sub_SeqFeature) {
	  my @subf = $f->sub_SeqFeature;

	  $self->find_supporting_evidence(\@subf);
	} else {
	  if ($f->seqname eq $self->contig_id) {
	    if (!($f->end < $self->start || $f->start > $self->end)) {
	      $self->add_Supporting_Feature($f);
	    }
	  }
	}
      }
}


# Inherited methods
# but you do have all the SeqFeature documentation: reproduced here
# for convenience...

=pod

=head1 Methods inherited from SeqFeature

=head2 start

 Title   : start
 Usage   : $start = $feat->start
 Function: Returns the start coordinate of the feature
 Returns : integer
 Args    : none

=head2 end

 Title   : end
 Usage   : $end = $feat->end
 Function: Returns the end coordinate of the feature
 Returns : integer
 Args    : none

=head2 length

 Title   : length
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

=head2 sub_SeqFeature

 Title   : sub_SeqFeature
 Usage   : @feats = $feat->sub_SeqFeature();
 Function: Returns an array of sub Sequence Features
 Returns : An array
 Args    : none

=head2 primary_tag

 Title   : primary_tag
 Usage   : $tag = $feat->primary_tag()
 Function: Returns the primary tag for a feature,
           eg 'exon'
 Returns : a string 
 Args    : none

=head2 source_tag

 Title   : source_tag
 Usage   : $tag = $feat->source_tag()
 Function: Returns the source tag for a feature,
           eg, 'genscan' 
 Returns : a string 
 Args    : none

=head2 has_tag

 Title   : has_tag
 Usage   : $value = $self->has_tag('some_tag')
 Function: Returns the value of the tag (undef if 
           none)
 Returns : 
 Args    :

=head2 all_tags

 Title   : all_tags
 Usage   : @tags = $feat->all_tags()
 Function: gives all tags for this feature
 Returns : an array of strings
 Args    : none

=head2 gff_string

 Title   : gff_string
 Usage   : $str = $feat->gff_string
 Function: provides the feature information in GFF
           version 2 format.
 Returns : A string
 Args    : None



=head1 RangeI methods

These methods are inherited from RangeI and can be used
directly from a SeqFeatureI interface. Remember that a 
SeqFeature is-a RangeI, and so wherever you see RangeI you
can use a feature ($r in the below documentation).

=head2 overlaps

  Title   : overlaps
  Usage   : if($feat->overlaps($r)) { do stuff }
            if($feat->overlaps(200)) { do stuff }
  Function: tests if $feat overlaps $r
  Args    : a RangeI to test for overlap with, or a point
  Returns : true if the Range overlaps with the feature, false otherwise


=head2 contains

  Title   : contains
  Usage   : if($feat->contains($r) { do stuff }
  Function: tests whether $feat totally contains $r
  Args    : a RangeI to test for being contained
  Returns : true if the argument is totaly contained within this range


=head2 equals

  Title   : equals
  Usage   : if($feat->equals($r))
  Function: test whether $feat has the same start, end, strand as $r
  Args    : a RangeI to test for equality
  Returns : true if they are describing the same range


=head1 Geometrical methods

These methods do things to the geometry of ranges, and return
triplets (start, stop, strand) from which new ranges could be built.

=head2

  Title   : intersection
  Usage   : ($start, $stop, $strand) = $feat->intersection($r)
  Function: gives the range that is contained by both ranges
  Args    : a RangeI to compare this one to
  Returns : nothing if they don''t overlap, or the range that they do overlap


=head2 union

  Title   : union
  Usage   : ($start, $stop, $strand) = $feat->union($r);
          : ($start, $stop, $strand) = Bio::RangeI->union(@ranges);
  Function: finds the minimal range that contains all of the ranges
  Args    : a range or list of ranges to find the union of
  Returns : the range containing all of the ranges

=cut

1;
