#
# BioPerl module for Bio::EnsEMBL::GeneComparison::GeneCompare
#
# Cared for by Simon Kay <sjk@sanger.ac.uk>
#
# Copyright Simon Kay
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::GeneComparison::GeneCompare - Comparison of a gene against an array of genes

=head1 SYNOPSIS

$self->{'_standardExons'} is an array of exons found on the standard gene.
These are compared to the array of exons found on the genes stored in
@{$self->{'_predictorGenes'}}.


=head1 DESCRIPTION

This module represents a class that is used for gene comparison. It is a specialisation 
of Bio::Root::RootI. The constructor (new method) takes an array of genes which is stored
in the @{$self->{'_predictorGenes'}} property. These are the genes against which all 
comparisons are made. An array of exons on these genes is obtained by calling the
_getPredictorExons method. The standard gene is set by the setStandardGene method which
finds all the standard exons and stores them in the $self->{'_standardExons'} property. 

Various methods perform comparisons between the standard gene and the array of predictor genes
by in turn creating objects of the class Bio::EnsEMBL::GeneComparison::ExonCompare and calling
appropriate methods on these:

isMissed returns 1 if none of the exons on the standard gene are overlapped by the predictor genes.
isExactlyMatched returns 1 if all the exons on the standard gene are exactly identified 
by only one of the predictor genes.
getGeneOverlapCount returns the number of predictor genes that the standard gene overlaps. 
getExactOverlapRatio returns a list. The first item is the number of exons on the standard 
gene that are exactly overlapped; the second is the number that are not.
getMissed returns the number of standard exons which are not overlapped by predictor exons.
getExonOverlaps returns an array containing predictor exons which the standard gene overlaps.
getBaseOverlaps returns the number of bases on the standard gene that have true positive, true negative
and false positive overlaps with all the overlapping exons from the predictor genes. 


=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

=cut


package Bio::EnsEMBL::GeneComparison::GeneCompare;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::GeneComparison::ExonCompare;
use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI);



=head2 new

 Title   : new
 Usage   : GeneCompare->new()
 Function: Constructor
 Example : 
 Returns : Reference to an object
 Args    : 

=cut

sub new {
    my ($class, @predictedGenes) = @_;
    my $self = bless {}, $class;
    @{$self->{'_predictorGenes'}} = @predictedGenes;           
    return $self;
}



=head2 setStandardGene

 Title   : setStandardGene
 Usage   : $obj->setStandardGene($standardGene)
 Function: Gets all the exons from the standard gene parameter and 
            stores them in $self->{'_standardExons'}. 
 Example : 
 Returns : 
 Args    : $standardGene - reference to a Gene object


=cut

sub setStandardGene {
    my ($self, $standardGene) = @_;
    
    $self->{'_standardExons'} = [$standardGene->each_unique_Exon()];
}



=head2 _getStandardExons

 Title   : _getStandardExons
 Usage   : $obj->_getStandardExons
 Function: 
 Example : 
 Returns : Array of standard exons
 Args    : 


=cut

sub _getStandardExons {
    my ($self) = @_;
    
    return @{$self->{'_standardExons'}}; 
}



=head2 _getPredictorGenes

 Title   : _getPredictorGenes
 Usage   : $obj->_getPredictorGenes
 Function: 
 Example : 
 Returns : Array of predictor genes
 Args    : 


=cut

sub _getPredictorGenes {
    my ($self) = @_;
    
    return @{$self->{'_predictorGenes'}}; 
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
    
    my @predictorExons;
    # Note that if different genes have the same exon, @predictorGenes will
    # include multiple copies of this exon but this is precisely what we want!
    foreach my $gene ($self->_getPredictorGenes()) {
        foreach my $exon ($gene->each_unique_Exon() ) {
            push(@predictorExons, $exon);
        }
    } 
    
    return @predictorExons;  
}



=head2 isMissed

 Title   : isMissed
 Usage   : $missed = $obj->isMissed()
 Function: A gene is considered missed if none of its exons are overlapped by any predicted exons.
 Example : 
 Returns : 1 or 0
 Args    : None


=cut

sub isMissed {
    my ($self) = @_;

    my $comparer = new Bio::EnsEMBL::GeneComparison::ExonCompare($self->_getPredictorExons());

    foreach my $standardExon ($self->_getStandardExons()) { 
        $comparer->setStandardExon($standardExon);            
        if ($comparer->hasOverlap()) {   
            return 0;
        }
    } 
    return 1;
}



=head2 isExactlyMatched

 Title   : isExactlyMatched
 Usage   : $obj->isExactlyMatched()
 Function: A gene is considered exactly matched if all the exons 
            on the standard gene are exactly identified by only one predictor gene
 Example : 
 Returns : 1 or 0 
 Args    : None


=cut

sub isExactlyMatched {
    my ($self) = @_;
     
    foreach my $gene ($self->_getPredictorGenes()) {
    
        my @predictorExons = $gene->each_unique_Exon();           
        my $comparer = new Bio::EnsEMBL::GeneComparison::ExonCompare(@predictorExons);    
        my $overlap = 1;
        
        foreach my $standardExon ($self->_getStandardExons()) {        
            $comparer->setStandardExon($standardExon);       
            unless ($comparer->hasExactOverlap()) { 
                $overlap = 0;
                # Break from inner loop as this predictor gene doesn't match         
                last;    
            }       
        } 
        if ($overlap) { 
            return 1;
        }
    }
   
    return 0;
}



=head2 getGeneOverlapCount

 Title   : getGeneOverlapCount
 Usage   : $obj->getGeneOverlapCount()
 Function: Calculates the number of predictor genes that the standard gene overlaps. 
 Example : 
 Returns : integer
 Args    : None


=cut

sub getGeneOverlapCount {
    my ($self) = @_;
    
    my $overlaps = 0;
    foreach my $gene ($self->_getPredictorGenes()) { 
        my $comparer = new Bio::EnsEMBL::GeneComparison::ExonCompare($gene->each_unique_Exon());
        foreach my $standardExon ($self->_getStandardExons()) {
            $comparer->setStandardExon($standardExon);            
            if ($comparer->hasOverlap()) {
                $overlaps++;
                last;
            }
        }
    }
    
    return $overlaps;
}



=head2 getExactOverlapRatio

 Title   : getExactOverlaps
 Usage   : $obj->getExactOverlapRatio()
 Function: Calculates the number of standard exons that are exactly overlapped and the number that are not.
 Example : 
 Returns : (integer, integer) - $overlaps, $nonOverlaps
 Args    : None


=cut

sub getExactOverlapRatio {
    my ($self) = @_;

    my $comparer = new Bio::EnsEMBL::GeneComparison::ExonCompare($self->_getPredictorExons());
    
    my $overlaps = 0;
    my $nonOverlaps = 0;
    foreach my $standardExon ($self->_getStandardExons()) {
    
        $comparer->setStandardExon($standardExon);            
        if ($comparer->hasExactOverlap()) {
            $overlaps++;
        } else {
            $nonOverlaps++;
        }
    }
    
    return ($overlaps, $nonOverlaps);
}



=head2 getMissed

 Title   : getMissed
 Usage   : $missed = $obj->getMissed()
 Function: The number of standard exons which are not overlapped by predictor exons.
 Example : 
 Returns : Integer
 Args    : None


=cut

sub getMissed {
    my ($self) = @_;

    my $comparer = new Bio::EnsEMBL::GeneComparison::ExonCompare($self->_getPredictorExons());

    my $missed = 0;
    my $count = 0;
    
    foreach my $standardExon ($self->_getStandardExons()) {
        $comparer->setStandardExon($standardExon);            
        unless ($comparer->hasOverlap()) {
            $missed++;
        }
        $count++
    }
    
    return ($missed, $count);
}



=head2 getExonOverlaps

 Title   : getExonOverlaps
 Usage   : @exons = $obj->getExonOverlaps()
 Function: Finds the predictor exons which the standard gene overlaps.
 Example : 
 Returns : Array of references to Bio::EnsEMBL::Exon objects.
 Args    : None


=cut

sub getExonOverlaps {
    my ($self) = @_;

    my $comparer = new Bio::EnsEMBL::GeneComparison::ExonCompare($self->_getPredictorExons());    
    my @overlaps = ();
    
    foreach my $standardExon ($self->_getStandardExons()) {
        $comparer->setStandardExon($standardExon); 
        @overlaps = (@overlaps, $comparer->getOverlaps());        
    }  
    
    return (@overlaps);
}



=head2 getBaseOverlaps

 Title   : getBaseOverlaps
 Usage   : $obj->getBaseOverlaps()
 Function: Calculates the number of bases involved in true positive, true negative
            and false positive overlaps between the standard gene and overlapping predictor genes. 
            True positive means that a standard base is found on a predictor exon;
            false negative means that a standard base is not found on a predictor exon;
            false positve means that a base is found on an overlapping predictor gene
            but not on the standard gene.
 Example : 
 Returns : (integer, integer, integer) - $truePositive, $falsePositive, $falseNegative 
 Args    : 


=cut

sub getBaseOverlaps {
    my ($self, $predictorExons) = @_;
        
    my $comparer = new Bio::EnsEMBL::GeneComparison::ExonCompare(@$predictorExons);
    
    my $truePositive = 0;
    my $falsePositive = 0;
    my $falseNegative = 0;
    
    foreach my $standardExon ($self->_getStandardExons()) {

        $comparer->setStandardExon($standardExon); 
        
        if ($comparer->hasOverlap) {           
            my ($tP, $fP, $fN) = $comparer->getBaseOverlaps();  
            $truePositive += $tP;
            $falsePositive += $fP;
            $falseNegative += $fN;
        }
        # If no overlaps are found then the whole length of the standard exon must be a false negative.
        else {
            $falseNegative += $standardExon->length;
        }
    }  
                     
    return ($truePositive, $falsePositive, $falseNegative);         
}

1;
