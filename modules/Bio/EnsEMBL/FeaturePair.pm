
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
sequence - e.g. a blast hit where residues 1-40 of a  protein sequence SW:HBA_HUMAN  
has hit to bases 100 - 220 on a genomic sequence HS120G22.  The genomic sequence 
coordinates are used to create one sequence feature $f1 and the protein coordinates
are used to create feature $f2.  A FeaturePair object can then be made

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

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::FeaturePair;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::EnsEMBL::SeqFeatureI;
use Bio::SeqFeature::FeaturePair;

@ISA = qw(Bio::SeqFeature::FeaturePair 
	Bio::EnsEMBL::SeqFeatureI);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize(@args);

  # set stuff in self from @args
  return $make; # success - we hope!
}


=head2 analysis

 Title   : analysis
 Usage   : $sf->analysis();
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub analysis{
   my ($self,$value) = @_;

   if( defined $value ) {
       $self->throw("Trying to add a non analysis object!") unless ref($value) eq 'Bio::EnsEMBL::Analysis::Analysis';
       $self->{_analysis} = $value;
   }

   if (defined($self->feature1)) { $self->feature1->analysis($value);}
   if (defined($self->feature2)) { $self->feature2->analysis($value);}

   if (defined($self->{_analysis})) {
       return $self->{_analysis};
   } else {
       return $self->feature1->analysis;
   }
   return $self->{_analysis};
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

    # First the features;

    $self->throw("Empty or wrong type of feature1 object") unless defined($self->feature1)   && 
	                                             ref($self->feature1) ne "" && 
						     $self->feature1->isa("Bio::EnsEMBL::SeqFeatureI");
    $self->throw("Empty or wrong type of feature1 object ") unless defined($self->feature2) &&
	                                             ref($self->feature2) ne "" && 
						     $self->feature2->isa("Bio::EnsEMBL::SeqFeatureI");

    $self->feature1->validate();
    $self->feature2->validate();

    # Now the analysis object
    if (defined($self->analysis)) {
	$self->throw("Wrong type of analysis object") unless ref($self->analysis) eq "Bio::EnsEMBL::Analysis::Analysis";
    } else {
	$self->throw("No analysis object defined");
    }
}

=head2 to_FTHelper

 Title   : to_FTHelper
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub to_FTHelper{
   my ($self) = @_;

   # Make new FTHelper, and fill in the key
   my $fth = Bio::AnnSeqIO::FTHelper->new;
   $fth->key('similarity');
   
   # Add location line
   my $g_start = $self->start;
   my $g_end   = $self->end;
   my $loc = "$g_start..$g_end";
   if ($self->strand == -1) {
        $loc = "complement($loc)";
    }
   $fth->loc($loc);
   
   # Add note describing similarity
   my $type    = $self->hseqname;
   my $r_start = $self->hstart;
   my $r_end   = $self->hend;
   $fth->add_field('note', "$type: matches $r_start to $r_end");
   $fth->add_field('note', "score=".$self->score);
   
   
   return $fth;
}


1;
