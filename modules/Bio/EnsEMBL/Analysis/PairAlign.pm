#
# BioPerl module for Transcript
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
   
    my $genomic = new Bio::SeqFeature::Homol  (-start  => $qstart,
					       -end    => $qend,
					       -strand => $qstrand);

    my $cdna     = new Bio::SeqFeature::Generic(-start => $hstart,
						-end   => $hend,
						-strand => $hstrand);

       $genomic->homol_SeqFeature($cdna);

    my $pair   = new Bio::EnsEMBL::Analysis::PairAlign;
       $pair->addPair($genomic);

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

sub addHomol {
    my ($self,$homol) = @_;

    $self->throw("Not a Bio::SeqFeature::Homol object") unless ($homol->isa("Bio::SeqFeature::Homol"));

    push(@{$self->{'_homol'}},$homol);
    
}


=head2 eachHomol

 Title   : eachHomol
 Usage   : my @homols = $pair->eachHomol
 Function: 
 Example : 
 Returns : Array of Bio::SeqFeature::Homol
 Args    : none


=cut

sub eachHomol {
    my ($self) = @_;

    if (defined($self->{'_homol'})) {
	return @{$self->{'_homol'}};
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
    my @homols = $self->eachHomol;

    @homols = sort {$a->start <=> $b->start} @homols;

  HOMOL: while (my $sf1 = shift(@homols)) {
      next HOMOL unless ($coord >= $sf1->start && $coord <= $sf1->end);
      
      my $sf2 = $sf1->homol_SeqFeature();

      if ($sf1->strand == 1 && $sf2->strand == 1) {
	  return ($sf2->start + ($coord - $sf1->start));
      } elsif ($sf1->strand == 1 && $sf2->strand == -1) {
	  return ($sf2->end   - ($coord - $sf1->start));
      } elsif ($sf1->strand == -1 && $sf2->strand == 1) {
	  return ($sf2->start + ($sf1->end - $coord));
      } elsif ($sf1->strand == -1 && $sf2->strand == -1) {
	  return ($sf2->end   - ($sf1->end - $coord));
      } else {
	  $self->throw("ERROR: Wrong strand value in homol (" . $sf1->strand . "/" . $sf2->strand . "\n");
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

    my @homols = $self->eachHomol;

  HOMOL: while (my $sf1 = shift(@homols)) {

      my $sf2 = $sf1->homol_SeqFeature();
      next HOMOL unless ($coord >= $sf2->start && $coord <= $sf2->end);

      if ($sf1->strand == 1 && $sf2->strand == 1) {
	  return ($sf1->start + ($coord - $sf2->start));
      } elsif ($sf1->strand == 1 && $sf2->strand == -1) {
	  return ($sf1->start + ($sf2->end -$coord));
      } elsif ($sf1->strand == -1 && $sf2->strand == 1) {
	  return ($sf1->end   - ($coord - $sf2->start));
      } elsif ($sf1->strand == -1 && $sf2->strand == -1) {
	  return ($sf1->end   - ($sf2->end - $coord));
      } else {
	  $self->throw("ERROR: Wrong strand value in homol (" . $sf1->strand . "/" . $sf2->strand . "\n");
      }
  }
    
}

1;








