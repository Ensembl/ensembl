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

Frameshifts - these haven't been coded yet :-)

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

@ISA = qw(Bio::SeqFeature::Generic Exporter);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  # add in EnsEMBL tag as 1.

  $self->primary_tag('exon');
  $self->source_tag('EnsEMBL');

  # Parse the input paramters (start,end,strand)
  if ($#args == 2) {
    $self->_parse(@args);
  }

  # set stuff in self from @args
  return $make; # success - we hope!
}

# Parse routine called from the constructor to
# set the basic variables start,end and strand

sub _parse {

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

=head2 dna_seq

  Title   : dna_seq
  Usage   : $dna = $feat->dna_seq
  Function: Returns the dna sequence of the exon
  Returns : Bio::Seq
  Args    : none

=cut

sub dna_seq {
  my ($self,$seq) = @_;

  # Bit more fussy on the input here.  Don't know if
  # we have to be so conscientious.
  if (defined($seq) && $seq->isa("Bio::Seq")) {
    $self->{'dna_seq'} = $seq;
  } elsif (defined($seq)) {
    $self->throw("Input to dna_seq is not a Bio::Seq");
  } 
  return $self->{'dna_seq'};

}

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

=head2 contigid

 Title   : contigid
 Usage   : $obj->contigid($newval)
 Function: 
 Returns : value of contigid
 Args    : newvalue (optional)


=cut

sub contigid{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'contigid'} = $value;
    }
    return $self->{'contigid'};

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

=head2 pep_seq

  Title   : pep_seq
  Usage   : @pep = $feat->pep_seq
  Function: Returns the 3 frame peptide translations of the exon
  Returns : @Bio::Seq
  Args    : none

=cut

sub pep_seq {
  my ($self) = @_;

  my @pep = $self->_translate();

  return @pep;
}

=pod 

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

  my $end = ($len - $phase) % 3;        # Jiggery pokery to find end phase
  $self->{'end_phase'} = (3 - $end)%3;  # Make sure only 0,1,2
  
  return $self->{'end_phase'};
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
    # Value must be 0,1,2,
    if ($value < 0 || $value > 2) {
      $self->throw("Bad value ($value) for exon phase. Should only be 0,1,2\n");
    } else {
      $self->{'phase'} = $value;
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
    $self->{'typ'} = $value;
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
  my @pep;
  my $i;
  
  # Get the DNA sequence and create the sequence string
  $self->dna_seq() || $self->throw("No DNA in object. Can't translate\n");

  my $dna = $self->dna_seq()->seq();
  
  # Translate in all frames - have to chop
  # off bases from the beginning of the dna sequence 
  # for frames 1 and 2 as the translate() method
  # only translates in one frame. Pah!
  
  for ($i = 0; $i < 3; $i++) {
    my $tmp = new Bio::Seq(-seq => substr($dna,$i));
    $pep[$i] = $tmp->translate();
  }
  
  return @pep;
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

  my @pep = $self->_translate() || throw ("Can't translate DNA\n");
  return $pep[$self->phase()];
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

# Inherited methods
# but you do have all the SeqFeature documentation: reproduced here
# for convience...

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
  Returns : nothing if they don't overlap, or the range that they do overlap


=head2 union

  Title   : union
  Usage   : ($start, $stop, $strand) = $feat->union($r);
          : ($start, $stop, $strand) = Bio::RangeI->union(@ranges);
  Function: finds the minimal range that contains all of the ranges
  Args    : a range or list of ranges to find the union of
  Returns : the range containing all of the ranges

=cut

1;
