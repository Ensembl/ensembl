
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

    # $analysis is a Bio::EnsEMBL::Analysis::Analysis object
    
    # SeqFeatureI methods can be used
    my $start = $feat->start;
    my $end   = $feat->end;

    # Bio::EnsEMBL::SeqFeature specific methods can be used
    my $analysis = $feat->analysis;

    # Validate all the data in the object
    $feat->validate  || $feat->throw("Invalid data in $feat");

=head1 DESCRIPTION

This is an implementation of the ensembl Bio::EnsEMBL::SeqFeatureI interface.  Extra
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

# Object preamble - inherits from Bio::Root::Object

use Bio::SeqFeatureI;
use Bio::Root::Object;

@ISA = qw(Bio::Root::Object Bio::SeqFeatureI);

# new is inherited from Bio::Root::Object

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  my($start,$end,$strand,$frame,$score,$analysis,$source_tag,$primary_tag,$seqname) = 
      $self->_rearrange([qw(START
			    END
			    STRAND
			    FRAME
			    SCORE
			    ANALYSIS
			    SOURCE_TAG
			    PRIMARY_TAG
			    SEQNAME
			    )],@args);

#  $gff_string && $self->_from_gff_string($gff_string);

  $start        && $self->start($start);
  $end          && $self->end($end);
  $strand       && $self->strand($strand);
  $primary_tag  && $self->primary_tag($primary_tag);
  $source_tag   && $self->source_tag($source_tag);
  $frame        && $self->frame($frame);
  $score        && $self->score($score);
  $analysis     && $self->analysis($analysis);
  $seqname      && $self->seqname($seqname);

  return $make; # success - we hope!

}

=head2 seqname

 Title   : seqname
 Usage   : $obj->seqname($newval)
 Function: There are many cases when you make a feature that you
           do know the sequence name, but do not know its actual
           sequence. This is an attribute such that you can store 
           the seqname.

           This attribute should *not* be used in GFF dumping, as
           that should come from the collection in which the seq
           feature was found.
 Returns : value of seqname
 Args    : newvalue (optional)


=cut

sub seqname{
   my ($self,$arg) = @_;

   if( $arg) {
      $self->{'_gsf_seqname'} = $arg;
    }

    return $self->{'_gsf_seqname'};

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

    if(defined($value)) {
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
	    $self->throw("$value is not a valid end");
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

sub strand{
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
  
    if (defined($value)) {
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
	$self->{'_primary_tag'} = $arg;
   }
   return $self->{'_primary_tag'};
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
	$self->{'_source_tag'} = $arg;
    }

   return $self->{'_source_tag'};
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
       $self->throw("Analysis is not a Bio::EnsEMBL::Analysis::Analysis object") unless 
	   ref($value) eq"Bio::EnsEMBL::Analysis::Analysis";
       $self->{_analysis} = $value;
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
   my ($self,$value) = @_;

   $self->_abstractDeath;

}


1;
