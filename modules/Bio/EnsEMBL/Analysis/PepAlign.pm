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

PepAlign - Dna/peptide pairwise alignment module

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Contains list of sub alignments making up a peptide-dna alignment

Creation:
   
    my $cdna     = new Bio::SeqFeature::Homol      (-start  => $qstart,
						    -end    => $qend,
						    -strand => $qstrand);

    my $pep      = new Bio::EnsEMBL::Analysis::pep_SeqFeature(-start => $hstart,
						    -end   => $hend,
						    -strand => $hstrand);
       $pep->start_frac(2);
       $pep->end_frac(3);

       $cdna->homol_SeqFeature($pep);

    my $pair   = new Bio::EnsEMBL::Analysis::PepAlign;
       $pair->addHomol($cdna);

Any number of pair alignments can be added to the PepAlign object


Manipulation:

To convert between coordinates : 

    my $cdna_coord         = $pair->pep2cDNA($cdna_coord);
    my ($pep_coord,$frac)  = $pair->cDNA2pep($cdna_coord);
  
To access the homols

    my @homols = $pair->eachHomol;

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Analysis::PepAlign;

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

=pod 

=head1 Methods unique to PepAlign


=pod 

=head2 addHomol

 Title   : addHomol
 Usage   : $pair->addHomol($homol)
 Function: Adds a homol feature to the alignment object
 Example : 
 Returns : nothing
 Args    : Bio::SeqFeature::Homol


=cut

sub addHomol {
    my ($self,$homol) = @_;

    $self->throw("Not a Bio::SeqFeature::Homol object") unless ($homol->isa( "Bio::SeqFeature::Homol"));

    push(@{$self->{'_homol'}},$homol);
    
}

=head2 eachHomol

 Title   : eachHomol
 Usage   : $pair->eachHomol
 Function: Returns an array of all the homols in the alignment object
 Example : 
 Returns : @Bio::SeqFeature::Homol
 Args    : none


=cut

sub eachHomol {
    my ($self) = @_;

    if (defined($self->{'_homol'})) {
	return @{$self->{'_homol'}};
    }
}


=head2 cDNA2pep

 Title   : cDNA2pep
 Usage   : my ($pep_coord,$frac) = $pair->cDNA2pep($pep_coord);
 Function: Converts a cdna coordinate into a peptide coord and fraction
 Example : 
 Returns : int,int
 Args    : int


=cut

sub cDNA2pep {
    my ($self,$coord) = @_;
    my @homols = $self->eachHomol;

       @homols = sort {$a->start <=> $b->start} @homols;


    # We loop through all the homols to see whether our cDNA coordinate
    # is contained therein

  HOMOL: while (my $sf1 = shift(@homols)) {

      next HOMOL unless ($coord >= $sf1->start && $coord <= $sf1->end);
      

      # We have found the homol our coordinate is in - we need
      # to now find the peptide coord corresponding to our cDNA

      # This is the peptide homol
      my $sf2 = $sf1->homol_SeqFeature();


      # We now have four different combinations of strands
      # to cope with

      if ($sf1->strand == 1 && $sf2->strand == 1) {
	  #  ---------------->  cDNA
	  #  -  -  -  -  -  ->  peptide

	  # codon_start and end are the cDNA start and end
	  # coordinates of the 1st/last FULL codon.

	  my $codon_start = $sf1->start + (4-$sf2->start_frac)%3;
	  my $codon_end   = $sf1->end   - (  $sf2->end_frac  )%3;


	  # pep_start is the peptide coordinate of the first
	  # FULL peptide residue
	  my $pep_start   = $sf2->start;

	  if ($sf2->start_frac != 1) {
	      $pep_start++;
	  }

	  # Similarly for pep_end
	  my $pep_end     = $sf2->end;
	  
	  if ($sf2->end_frac != 3) {
	      $pep_end--;
	  }


	  # We have to deal with cDNA coordinates that lie
	  # outside codon_start/codon_end differently

	  if ($coord < $codon_start) {
	      return ($sf2->start,$coord-$sf1->start+ $sf2->start_frac);
	  } elsif ($coord > $codon_end) {
	      return ($sf2->end, $coord-$codon_end);
	  }


	  my $chop = (4 - $sf2->start_frac)%3;


	  my $ncodons = int(($coord - $chop - $sf1->start)/3);
	  my $frac    =     ($coord - $chop - $sf1->start)%3;

	  $ncodons++ if ($sf2->start_frac != 1);
	  $frac++;

	  my $pep = $sf2->start + $ncodons;

	  return ($pep,$frac);

      } elsif ($sf1->strand == 1 && $sf2->strand == -1) {
	  #  ---------------->  cDNA
	  #  <  -  -  -  -  -   peptide

	  # codon_start and end are the cDNA start and end
	  # coordinates of the 1st/last FULL codon.

	  my $codon_start = $sf1->start + $sf2->end_frac;
	  my $codon_end   = $sf1->end   - (4 - $sf2->start_frac)%3;

	  # pep_start is the peptide coordinate of the first
	  # FULL peptide residue
	  my $pep_start   = $sf2->start;

	  if ($sf2->start_frac != 1) {
	      $pep_start++;
	  }

	  # Similarly for pep_end
	  my $pep_end     = $sf2->end;
	  
	  if ($sf2->end_frac != 3) {
	      $pep_end--;
	  }

	  # We have to deal with cDNA coordinates that lie
	  # outside codon_start/codon_end differently
	  
	  if ($coord < $codon_start) {
	      return ($sf2->end,$codon_start-$coord);
	  } elsif ($coord > $codon_end) {
	      return ($sf2->start, $sf2->start_frac + $sf1->end - $coord);
	  }

	  my $chop = (4 - $sf2->start_frac)%3;

	  my $ncodons = int((-$coord - $chop +  $sf1->end)/3);
	  my $frac    =     (-$coord - $chop + $sf1->end)%3;

	  $ncodons++ if ($sf2->start_frac != 1);
	  $frac++;


	  my $pep = $sf2->start  + $ncodons;
	  

	  return ($pep,$frac);

      } elsif ($sf1->strand == -1 && $sf2->strand == 1) {
	  $self->throw("not implemented yet\n");

      } elsif ($sf1->strand == -1 && $sf2->strand == -1) {
	  $self->throw("not implemented yet\n");

      } else {
	  $self->throw("ERROR: Wrong strand value in homol (" . $sf1->strand . "/" . $sf2->strand . "\n");
      }
  }

}

=head2 pep2cDNA

 Title   : pep2cDNA
 Usage   : my ($cdna_coord) = $pair->pep2cDNA($pep_coord,$frac);
 Function: Converts a peptide coordinate and fraction into a cdna coordinate
 Example : 
 Returns : int
 Args    : int,int


=cut

sub pep2cDNA {
    my ($self,$coord,$frac) = @_;

    my @homols = $self->eachHomol;

    $frac = 1 unless $frac;


  HOMOL: while (my $sf1 = shift(@homols)) {
      # This is the peptide homol
      my $sf2 = $sf1->homol_SeqFeature();
      
      next HOMOL unless ($coord >= $sf2->start && $coord <= $sf2->end);
      next HOMOL if     ($coord == $sf2->start && $frac  < $sf2->start_frac);
      next HOMOL if     ($coord == $sf2->end   && $frac  > $sf2->end_frac);

      # We have found the homol our coordinate is in - we need
      # to now find the cDNA coord corresponding to our peptide coord


      # We now have four different combinations of strands
      # to cope with

      if ($sf1->strand == 1 && $sf2->strand == 1) {
	  #  ---------------->  cDNA
	  #  -  -  -  -  -  ->  peptide

	  # codon_start and end are the cDNA start and end
	  # coordinates of the 1st/last FULL codon.

	  my $codon_start = $sf1->start + (4-$sf2->start_frac)%3;
	  my $codon_end   = $sf1->end   - (  $sf2->end_frac  )%3;


	  # pep_start is the peptide coordinate of the first
	  # FULL peptide residue
	  my $pep_start   = $sf2->start;

	  if ($sf2->start_frac != 1) {
	      $pep_start++;
	  }

	  # Similarly for pep_end
	  my $pep_end     = $sf2->end;
	  
	  if ($sf2->end_frac != 3) {
	      $pep_end--;
	  }


	  # We have to deal with cDNA coordinates that lie
	  # outside codon_start/codon_end differently

	  if ($coord < $pep_start) {
	      return ($sf1->start + ($frac - $sf2->start_frac));
	  } elsif ($coord > $pep_end) {
	      return ($sf1->end -  ($sf2->end_frac - $frac));
	  }

	  my $chop = 4 - $sf2->start_frac;
	  if ($sf2->start_frac == 1) { $chop = 0;}

	  my $nbases = ($coord - $pep_start)* 3 + $chop;


	  my $cdna = $sf1->start + $nbases;

	  if ($frac == 2) { $cdna += 1;}
	  if ($frac == 3) { $cdna += 2;}

#	  print("Found cdna $cdna\n");
	  return ($cdna);


      } elsif ($sf1->strand == 1 && $sf2->strand == -1) {
	  #  ---------------->  cDNA
	  #  <  -  -  -  -  -   peptide

	  # codon_start and end are the cDNA start and end
	  # coordinates of the 1st/last FULL codon.

	  my $codon_start = $sf1->start + (4-$sf2->start_frac)%3;
	  my $codon_end   = $sf1->end   - (  $sf2->end_frac  )%3;


	  # pep_start is the peptide coordinate of the first
	  # FULL peptide residue
	  my $pep_start   = $sf2->start;

	  if ($sf2->start_frac != 1) {
	      $pep_start++;
	  }

	  # Similarly for pep_end
	  my $pep_end     = $sf2->end;
	  
	  if ($sf2->end_frac != 3) {
	      $pep_end--;
	  }


	  # We have to deal with cDNA coordinates that lie
	  # outside codon_start/codon_end differently

	  if ($coord < $pep_start) {
	      return ($sf1->end  - ($frac - $sf2->start_frac));
	  } elsif ($coord > $pep_end) {
	      return ($sf1->start +  ($sf2->end_frac - $frac));
	  }

	  my $chop = 4 - $sf2->start_frac;
	  
	  if ($sf2->start_frac == 1) { $chop = 0;}

	  my $nbases = ($coord - $pep_start)* 3 + $chop;


	  my $cdna = $sf1->end - $nbases;

	  if ($frac == 2) { $cdna -= 1;}
	  if ($frac == 3) { $cdna -= 2;}

#	  print("Found cdna $cdna\n");
	  return ($cdna);


      } elsif ($sf1->strand == -1 && $sf2->strand == 1) {
	  $self->throw("Not implemented yet\n");

      } elsif ($sf1->strand == -1 && $sf2->strand == -1) {
	  $self->throw("Not implemented yet\n");

      } else {
	  $self->throw("ERROR: Wrong strand value in homol (" . $sf1->strand . "/" . $sf2->strand . "\n");
      }
  }
    
}

1;








