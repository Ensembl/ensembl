#
# BioPerl module for PairAlign object
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

PairAlign - Dna pairwise alignment module

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Contains list of sub alignments making up a dna-dna alignment

Creation:
   
    my $genomic = new Bio::EnsEMBL::SeqFeature(-start  => $qstart,
					       -end    => $qend,
					       -strand => $qstrand);

    my $cdna     = new Bio::EnsEMBL::SeqFeature(-start => $hstart,
						-end   => $hend,
						-strand => $hstrand);

    my $pair     = new Bio::EnsEMBL::FeaturePair(-feature1 => $genomic,
						 -feature2 => $cdna,
						 );

    my $pairaln   = new Bio::EnsEMBL::Analysis::PairAlign;
       $pairaln->addFeaturePair($pair);

Any number of pair alignments can be added to the PairAlign object


Manipulation:

To convert between coordinates : 

    my $cdna_coord = $pair->genomic2cDNA($gen_coord);
    my $gen_coord  = $pair->cDNA2genomic($cdna_coord);

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Analysis::PairAlign;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Object;

@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
    my($self,@args) = @_;
    
    $self->{'_homol'} = [];
    my $make = $self->SUPER::_initialize;
    
    return $self; # success - we hope!
}

sub addFeaturePair {
    my ($self,$pair) = @_;

    $self->throw("Not a Bio::SeqFeature::FeaturePair object") unless ($pair->isa("Bio::SeqFeature::FeaturePair"));

    push(@{$self->{'_pairs'}},$pair);
    
}


=head2 eachFeaturePair

 Title   : eachFeaturePait
 Usage   : my @pairs = $pair->eachFeaturePair
 Function: 
 Example : 
 Returns : Array of Bio::SeqFeature::FeaturePair
 Args    : none


=cut

sub eachFeaturePair {
    my ($self) = @_;

    if (defined($self->{'_pairs'})) {
	return @{$self->{'_pairs'}};
    }
}

=head2 genomic2cDNA

 Title   : genomic2cDNA
 Usage   : my $cdna_coord = $pair->genomic2cDNA($gen_coord)
 Function: Converts a genomic coordinate to a cdna coordinate
 Example : 
 Returns : int
 Args    : int


=cut

sub genomic2cDNA {
    my ($self,$coord) = @_;
    my @pairs = $self->eachFeaturePair;

    @pairs = sort {$a->start <=> $b->start} @pairs;

  HOMOL: while (my $sf1 = shift(@pairs)) {
      next HOMOL unless ($coord >= $sf1->start && $coord <= $sf1->end);
      
      if ($sf1->strand == 1 && $sf1->hstrand == 1) {
	  return ($sf1->hstart + ($coord - $sf1->start));
      } elsif ($sf1->strand == 1 && $sf1->hstrand == -1) {
	  return ($sf1->hend   - ($coord - $sf1->start));
      } elsif ($sf1->strand == -1 && $sf1->hstrand == 1) {
	  return ($sf1->hstart + ($sf1->end - $coord));
      } elsif ($sf1->strand == -1 && $sf1->hstrand == -1) {
	  return ($sf1->hend   - ($sf1->end - $coord));
      } else {
	  $self->throw("ERROR: Wrong strand value in FeaturePair (" . $sf1->strand . "/" . $sf1->hstrand . "\n");
      }
  }    
}

=head2 cDNA2genomic

 Title   : cDNA2genomic
 Usage   : my $gen_coord = $pair->genomic2cDNA($cdna_coord)
 Function: Converts a cdna coordinate to a genomic coordinate
 Example : 
 Returns : int
 Args    : int


=cut

sub cDNA2genomic {
    my ($self,$coord) = @_;

    my @pairs = $self->eachFeaturePair;

  HOMOL: while (my $sf1 = shift(@pairs)) {

      next HOMOL unless ($coord >= $sf1->hstart && $coord <= $sf1->hend);

      if ($sf1->strand == 1 && $sf1->hstrand == 1) {
	  return ($sf1->start + ($coord - $sf1->hstart));
      } elsif ($sf1->strand == 1 && $sf1->hstrand == -1) {
	  return ($sf1->start + ($sf1->hend -$coord));
      } elsif ($sf1->strand == -1 && $sf1->hstrand == 1) {
	  return ($sf1->end   - ($coord - $sf1->hstart));
      } elsif ($sf1->strand == -1 && $sf1->hstrand == -1) {
	  return ($sf1->end   - ($sf1->hend - $coord));
      } else {
	  $self->throw("ERROR: Wrong strand value in homol (" . $sf1->strand . "/" . $sf1->hstrand . "\n");
      }
  }
    
}

1;








