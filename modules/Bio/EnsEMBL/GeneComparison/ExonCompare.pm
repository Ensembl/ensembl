#
# BioPerl module for Bio::EnsEMBL::GeneComparison::ExonCompare
#
# Cared for by Simon Kay <sjk@sanger.ac.uk>
#
# Copyright Simon Kay
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::ExonCompare - Comparison of an exon against an array of exons

=head1 SYNOPSIS

$self->{'_standardExon'} is an exon being compared to an array of exons stored in
@{$self->{'_predictorExons'}}.


=head1 DESCRIPTION

This module represents a class that is used for gene comparison. It is a specialisation 
of Bio::Root::RootI. The constructor (new method) takes an array of exons which is stored
in the @{$self->{'_predictorExons'}} property. These are the exons against which all 
comparisons are made. The standard exon is set by the setStandardExon method. Four methods
perform comparisons between the standard exon and the array of predictor exons:
hasOverlap, hasExactOverlap, getOverlaps and getBaseOverlaps.



=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

=cut



package Bio::EnsEMBL::GeneComparison::ExonCompare;

use strict;
use vars qw(@ISA);
use Bio::Root::RootI;
use Bio::EnsEMBL::Exon;

@ISA = qw(Bio::Root::RootI);



=head2 new

 Title   : new
 Usage   : ExonCompare->new()
 Function: Constructor
 Example : 
 Returns : Reference to an object
 Args    : @exons - Array of exon objects to be used as predictors.

=cut

sub new {
    my ($class, @exons) = @_;
    my $self = bless {}, $class;
    @{$self->{'_predictorExons'}} = @exons;
    return $self;
}



=head2 setStandardExon

 Title   : setStandardExon
 Usage   : $obj->setStandardExon($exon)
 Function: 
 Example : 
 Returns : 
 Args    : Reference to a Bio::EnsEMBL::Exon


=cut

sub setStandardExon {
    my ($self, $exon) = @_;
    
    $self->throw("$exon is not an Exon") unless ($exon->isa('Bio::EnsEMBL::Exon'));
    $self->{'_standardExon'} = $exon;
}



=head2 _getStandardExon

 Title   : _getStandardExon
 Usage   : $obj->_getStandardExon()
 Function: 
 Example : 
 Returns : The standard exon
 Args    : 


=cut

sub _getStandardExon {
    my ($self) = @_;
    
    return $self->{'_standardExon'}; 
}



=head2 _getPredictorExons

 Title   : _getPredictorExons
 Usage   : $obj->_getPredictorExons
 Function: 
 Example : 
 Returns : Array of predictor exons
 Args    : 


=cut

sub _getPredictorExons {
    my ($self) = @_;
    
    return @{$self->{'_predictorExons'}}; 
}



=head2 hasOverlap

 Title   : hasOverlap
 Usage   : $overlap = $obj->hasOverlap()
 Function: Checks if any of the predictor exons have an overlap with the standard exon.
 Example : 
 Returns : 0 if none of the of the predictor exons overlap the standard, otherwise 1
 Args    : 


=cut

sub hasOverlap {
    my ($self) = @_;

    my $exon1 = $self->_getStandardExon();
    
    foreach my $exon2 ($self->_getPredictorExons()) {  
              
         if (($exon2->end   > $exon1->start && $exon2->start < $exon1->end) ||
             ($exon2->start < $exon1->end   && $exon2->end   > $exon1->start)) {                
                 return 1;
        }
    }
    
    return 0;
}



=head2 hasExactOverlap

 Title   : hasExactOverlap
 Usage   : $obj->hasExactOverlap()
 Function: Checks if any of the predictor exons exactly overlaps the standard exon.
 Example : 
 Returns : 0 if none of the predictor exons exactly matches the standard exon, otherwise 1.
 Args    : 


=cut

sub hasExactOverlap {
    my ($self) = @_;
    
    my $exon1 = $self->_getStandardExon();
   
    foreach my $exon2 ($self->_getPredictorExons()) {     
        if (($exon2->start == $exon1->start) && ($exon2->end == $exon1->end)) {
            return 1;
        }
    }
  
    return 0;
}



=head2 getOverlaps

 Title   : getOverlaps
 Usage   : $obj->getOverlaps()
 Function: Find all the predictor exons that the standard exons overlaps 
            adds them to an array and returns it.
 Example : 
 Returns : An array of references to Bio::EnsEMBL::Exon.
 Args    : 


=cut

sub getOverlaps {
    my ($self) = @_;
    
    my @overlaps = ();
    my $exon1 = $self->_getStandardExon();
    
    foreach my $exon2 ($self->_getPredictorExons()) { 
         if (($exon2->end   > $exon1->start && $exon2->start < $exon1->end) ||
             ($exon2->start < $exon1->end   && $exon2->end   > $exon1->start)) {   
                 push @overlaps, $exon2;
        }
    }
  
    return @overlaps;
}
   


=head2 getBaseOverlaps

 Title   : getBaseOverlaps
 Usage   : $obj->getBaseOverlaps()
 Function: Calculates the number of bases involved in true positive, true negative
            and false positive overlaps between the standard exon and overlapping predictor exons. 
            True positive means that a standard base is found on a predictor exon;
            false negative means that a standard base is not found on a predictor exon;
            false positve means that a base is found on an overlapping predictor exon
            but not on the standard exon.
 Example : 
 Returns : List of 3 integers.
 Args    : 


=cut

sub getBaseOverlaps {
    my ($self) = @_;
    
    my $exon1 = $self->_getStandardExon();
    my @overlaps = ();
    my $truePositive = $exon1->length();
    my $falsePositive = 0;
    my $falseNegative = 0;            
    
    # Find all the predictor exons with overlaps 
    foreach my $exon2 ($self->_getPredictorExons()) { 
        if (($exon2->end   > $exon1->start && $exon2->start < $exon1->end) ||
            ($exon2->start < $exon1->end   && $exon2->end   > $exon1->start)) {
                push @overlaps, $exon2;        
        }
    }
            
    # Sort the predictor exons in order of lowest start       
    @overlaps = sort { $a->start <=> $b->start } @overlaps;
    # Get the predictor exon with the lowest start.
    my ($exon2) = @overlaps;
    
    # Correct the scores for the difference between the standard and predictor starts
    my $leftEnd = $exon1->start - $exon2->start;
    if ($leftEnd > 0) {
        $falsePositive += $leftEnd;
    } else {
        $truePositive += $leftEnd;
        $falseNegative -= $leftEnd;
    }    
    
    # Sort the predictor exons in order of highest end
    @overlaps = sort { $b->end <=> $a->end } @overlaps;
    # Get the predictor exon with the highest end.
    ($exon2) = @overlaps;    
 
    # Correct the scores for the difference between the predictor and standard ends.
    my $rightEnd = $exon2->end - $exon1->end;
    if ($rightEnd > 0) {
        $falsePositive += $rightEnd;
    } else {
        $falseNegative += -$rightEnd;
        $truePositive += $rightEnd;
    }   
    
    return ($truePositive, $falsePositive, $falseNegative);
}


1;
