#
# bioperl module for a genscan predicted peptide and its blast hits
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::GenscanPeptide - stores blast hits for a genscan peptide

=head1 SYNOPSIS

    my $gs = new Bio::EnsEMBL::Analysis::GenscanPeptide

Inputting blast hits

Extracting data


=head1 DESCRIPTION

Genscan peptide object.  Stores blast hits to the genscan peptide and
converts the coordinates into the genomic dna frame. Two conversion
objects are used : 

   1. Bio::EnsEMBL::Analysis::PepAlign - converts the coordinates between the genomic and the genscan peptide

   2. Bio::EnsEMBL::Analysis::PepAlign/PairAlign - converts the coordinates between the genscan peptide and the blast hit.

East blast hit to the peptide is converted into an array of Bio::SeqFeature::Homol if the hit
spans across more than one exon.
               g1    g2
 ------------------------------------------------> genomic
      \   \ |        | /      /
       \   \|        |/      /
        --------------------->    genscan predicted gene with 3 exons
               |         |
               |         |        
               <----------        blast hit to the genscan peptide across 2 exons of the predicted gene
              h1     h2
                     

For the hit to the 2nd exon we need to know the start and
end points in the genomic sequence and in the blast hit sequence.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::GenscanPeptide;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object;

use Bio::Root::Object;

use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Analysis::PairAlign;
use Bio::EnsEMBL::Analysis::PepAlign;
use Bio::EnsEMBL::Pep_SeqFeature;
use Bio::EnsEMBL::FeaturePair;

# Inherits from the base bioperl object
@ISA = qw(Bio::Root::Object);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  
  my $make = $self->SUPER::_initialize;

  # Input variables
  # ---------------
  # Bio::EnsEMBL::Transcript object - these are created by the Bio::EnsEMBL::Analysis::GenscanMC module

  my $transcript  = shift(@args);
  
  $self->throw("Input argument is not a Bio::EnsEMBL::Transcript") unless $transcript->isa("Bio::EnsEMBL::Transcript");

  $self->{_transcript} = $transcript;

  # Make the peptide alignment object
  $self->_make_pepAlign;
  
  # Stored data
  # -----------
  # These are the blast hits to the peptide.

  $self->{_pepHits} = [];
  $self->{_dnaHits} = [];

  return $self; # success - we hope!
}

=head2 _make_pepAlign 

  Title   : _make_pepAlign
  Usage   : $self->_make_pepAlign
  Function: Converts the transcript to a Bio::EnsEMBL::Analysis::PepAlign object
  Returns : Nothing
  Args    : None

=cut

sub _make_pepAlign {
    my ($self) = @_;

    # This works.

    $self->throw("No transcript in GenscanPeptide object. Can't convert to pepAlign") unless $self->{_transcript};

    my $g = $self->{_transcript};

    my $prev_coord = 0;           # Peptide coordinate of the end of the previous exon
    my $prev_phase = 0;           # Peptide end phase of the previous exon
    
    my $pepaln = new Bio::EnsEMBL::Analysis::PepAlign;
    
    foreach my $ex ($g->each_Exon) {
	
	my $orient = $ex->strand;
	

	# Find the peptide coords corresponding to this exon
	my ($pepstart,$pep_startphase,$pepend,$pep_endphase) = $self->_find_coord($ex,$prev_coord,$prev_phase);


	# Create the dna homol 

	my $dnah   = new Bio::EnsEMBL::SeqFeature               (-start  => $ex->start,
								 -end    => $ex->end,
								 -strand => 1);

	# Create the peptide homol
	my $peph  = new Bio::EnsEMBL::Pep_SeqFeature(-start  => $pepstart,
						     -end    => $pepend,
						     -strand => $orient,
						     -start_frac => $pep_startphase,
						     -end_frac   => $pep_endphase);
	
	# Add the peptide homol to the dna homol
	my $pair = new Bio::EnsEMBL::FeaturePair(-feature1 => $dnah,
						 -feature2 => $peph
						 );

	# Finally add the dna homol to the align object
	$pepaln->addFeaturePair($pair);

	# A bit of printing
#	print("\n" . $ex->id . " coord\t" .$ex->start . "\t" . $ex->end . "\n");
#	print("Pep start\t$pepstart\t$pep_startphase\n");
#	print("Pep end  \t$pepend\t$pep_endphase\n");


	# These store the previous coords for the peptide
	$prev_phase = $pep_endphase;
	$prev_coord = $pepend;

    }

    $self->{_pepaln} = $pepaln;
}

=head2 _find_coord

  Title   : _find_coord
  Usage   : my ($pepstart,$pep_startphase,$pepend,$pep_endphase) = $self->_find_coord($ex,$prev,$prevphase)
  Function: Converts genomic exon start/end coords to peptide start/end coords
            The start and end phases are 1,2 or 3
  Returns : int,int,int,int
  Args    : None

=cut

sub _find_coord {
    my ($self,$ex,$prev,$prevphase) = @_;

    my $pepstart;
    my $pep_startphase;

    my $pepend;
    my $pep_endphase;
    
    if ($prev == 0) {
	$pepstart = 1;
	$pep_startphase = $ex->phase + 1;
    } else {
	
	my $inc = 0;

	$pep_startphase = ++$prevphase;
	
	if ($pep_startphase == 4) {
	    $pep_startphase = 1;
	    $inc = 1;
	}
	
	$pepstart = $prev + $inc;
    }
    
    my $exlen = abs($ex->end - $ex->start) + 1;
       $exlen -= (4-$pep_startphase)%3;
    
    my $nopep = int($exlen/3);
       $pep_endphase = $exlen%3 + 1;
    
    $pepend    = $nopep + $pepstart -1;
    
    return ($pepstart,$pep_startphase,$pepend,$pep_endphase);
}


=head2 _pepHit2homol

  Title   : _pepHit2homol
  Usage   : my @homols = $self->_pepHit2homol
  Function: Converts a protein hit to the peptide to an array of genomic homols
  Returns : Array of Bio::SeqFeature::Homol;
  Args    : None

=cut

sub _pepHit2homol {
    my ($self,$pephit) = @_;




    my $pepaln  = $self->{_pepaln};

    my $pairaln = new Bio::EnsEMBL::Analysis::PairAlign;   # Make a coordinate conversion object from
       $pairaln->addFeaturePair($pephit);                        # the peptide exon and the homol.

    my @pep     = $pepaln->eachFeaturePair;
    my @homols;


#    print "Debug: looking at peptide",$pephit->seqname,"\n";

    # This is dodgy

    foreach my $pep_homol ($pairaln->eachFeaturePair) {
	my $pep_homol2 = $pep_homol->feature2;


	# We now need to loop over each pepaln exon to see which
	# ones the homol overlaps

	# Find the first exon and process
	while (my $gen_exon = shift @pep) {
	    my $pep_exon = $gen_exon->feature2;
	    
	    if ($pep_homol->start >= $pep_exon->start && $pep_homol->start <= $pep_exon->end) {
		my $pstart = $pep_homol->start;
		my $pend   = $pep_exon ->end;
		
		my $start_frac = $pep_exon->start_frac;
		my $end_frac;
		
		if ($pep_homol->end <= $pep_exon->end) {
		    $pend     = $pep_homol->end;
		    if( $pep_homol->end == $pep_exon->end ) {
			$end_frac = $pep_exon->end_frac;
		    } else {
			$end_frac = 3;
		    }

		} else {
		    $end_frac = $pep_exon->end_frac;
		}
		
#		print("DEBUG:: $pstart $pend :" .  $pep_homol->start . " " . $pep_homol->end . "\n");

		# We know whereabouts on the peptide we want the homol to lie
		# We need to convert these coords into genomic coords and 
		# also homol coords.

		my $newh = $self->_make_homol($pairaln,$gen_exon,$pep_exon,$pep_homol,$pstart,$start_frac,$pend,$end_frac);

		push(@homols,$newh);

		last;
	    }
	}

	# Process the rest of the exons
	while (my $gen_exon = shift @pep) {
	    # Check for overlap
	    my $pep_exon     = $gen_exon->feature2;

	    my $pstart = $pep_exon->start;
	    my $pend   = $pep_exon->end;

	    last if ($pep_homol->end < $pep_exon->start);
	    
	    my $start_frac = $pep_exon->start_frac;
	    my $end_frac;
	    
	    if ($pep_homol->end <= $pep_exon->end) {
		$pend     = $pep_homol->end;
		if( $pep_homol->end == $pep_exon->end ) {
		    $end_frac = $pep_exon->end_frac;
		} else {
		    $end_frac = 3;
		}
	    } else {
		$end_frac = $pep_exon->end_frac;
	    }

		    
	    # We know whereabouts on the peptide we want the homol to lie
	    # We need to convert these coords into genomic coords and 
	    # also homol coords.

	    my $newh = $self->_make_homol($pairaln,$gen_exon,$pep_exon,$pep_homol,$pstart,$start_frac,$pend,$end_frac);
	    
	    push(@homols,$newh);
	    
	}
	
    }

    return @homols;
}

=head2 each_Homol

  Title   : each_Homol
  Usage   : my @homols = $self->each_Homol
  Function: Converts the pep and dna hits to genomic homols
  Returns : array of Bio::SeqFeature::Homol
  Args    : none

=cut

sub each_Homol {
    my ($self) = @_;
    my @homols;

    foreach my $pep (@{$self->{_pepHits}}) {
	my @newh = $self->_pepHit2homol($pep);
	push(@homols,@newh);
    }

    foreach my $dna (@{$self->{_dnaHits}}) {
	my @newh = $self->_dnaHit2homol($dna);
	push(@homols,@newh);
    }
    return @homols;
}

=head2 _make_homol

  Title   : _make_homol
  Usage   : my $homol = $self->_make_homol($pairaln,$gen_exon,$pep_exon$pep_homol,$pstart,$start_frac,$pend,$end_frac);
  Function: Converts a range of coordinates in a pairAlign object to a homol
  Returns : Bio::SeqFeature::Homol
  Args    : Bio::EnsEMBL::Analysis::pep_SeqFeature

=cut


sub _make_homol {
    my ($self,$pairaln,$gen_exon,$pep_exon,$pep_homol,$pstart,$start_frac,$pend,$end_frac) = @_;

    my $pepaln = $self->{_pepaln};

    my $pep_homol2 = $pep_homol->feature2;

    my $g1 = $pepaln->pep2cDNA($pstart,$start_frac);
    my $g2 = $pepaln->pep2cDNA($pend,  $end_frac);

    my $h1 = $pairaln->genomic2cDNA($pstart);
    my $h2 = $pairaln->genomic2cDNA($pend);


    if ($g1 > $g2) {
	my $tmp = $g1;
	$g1 = $g2;
	$g2 = $tmp;
    }
    
    if ($h1 > $h2) {
	my $tmp = $h1;
	$h1 = $h2;
	$h2 = $tmp;
    }
    
    my $newh   = new Bio::EnsEMBL::SeqFeature(-start  => $g1,
					      -end    => $g2,
					      -strand => $pep_exon->strand,
					      );

    $newh->primary_tag($pep_homol->primary_tag);
    $newh->source_tag ($pep_homol->source_tag);
    $newh->seqname    ($pep_homol->seqname);
    $newh->score      ($pep_homol->score);

    if ($pep_homol->isa("Bio::EnsEMBL::SeqFeatureI") && defined($pep_homol->analysis) && $newh->isa("Bio::EnsEMBL::SeqFeatureI")) {
	$newh->analysis($pep_homol->analysis);
    } else {
	$self->throw("Can't find analysis for pep_homol " . $pep_homol->seqname . " " . $pep_homol->source_tag . "\n");
    }
    # Create the peptide homol
    my $peph  = new Bio::EnsEMBL::Pep_SeqFeature(-start  => $h1,
						 -end    => $h2,
						 -strand => $pep_exon->strand);
    $peph->primary_tag($pep_homol2->primary_tag);
    $peph->source_tag ($pep_homol2->source_tag);
    $peph->seqname    ($pep_homol2->seqname);
    $peph->start_frac ($start_frac);
    $peph->end_frac   ($end_frac);
    $peph->score      ($pep_homol->score);

    if ($pep_homol2->isa("Bio::EnsEMBL::SeqFeatureI") && defined($pep_homol2->analysis) && $peph->isa("Bio::EnsEMBL::SeqFeatureI")) {
	$peph->analysis($pep_homol2->analysis);
    } else {
	$self->throw("Can't find analysis for pep_homol " . $pep_homol2->seqname . " " . $pep_homol2->source_tag . "\n");
    }
	    
    # Add the peptide homol to the dna homol
    my $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $newh,
					   -feature2 => $peph,
					   );
    $fp->analysis($pep_homol->analysis);

    return $fp;
}

=head2 _make_dna_homol

  Title   : _make_dna_homol
  Usage   : my $homol = $self->_make_dna_homol($pairaln,$pstart,$start_frac,$pend,$end_frac);
  Function: Converts a range of coordinates in a pairAlign object to a homol
  Returns : Bio::SeqFeature::Homol
  Args    : Bio::EnsEMBL::Analysis::pep_SeqFeature

=cut


sub _make_dna_homol {
    my ($self,$pairaln,$gen_exon,$pep_exon,$pep_homol,$pep_homol2,$pstart,$start_frac,$pend,$end_frac) = @_;

    my $pepaln = $self->{_pepaln};

    my $g1  = $pepaln->pep2cDNA($pstart,$start_frac);
    my $g2  = $pepaln->pep2cDNA($pend,  $end_frac);
    
    my $h1 = $pairaln->pep2cDNA($pstart,$start_frac);
    my $h2 = $pairaln->pep2cDNA($pend,  $end_frac);

    if ($g1 > $g2) {
	my $tmp = $g1;
	$g1 = $g2;
	$g2 = $tmp;
    }
    
    if ($h1 > $h2) {
	my $tmp = $h1;
	$h1 = $h2;
	$h2 = $tmp;
    }

    my $newh   = new Bio::EnsEMBL::SeqFeature(-start  => $g1,
					      -end    => $g2,
					      -strand => $pep_exon->strand * $pep_homol2->strand);

    $newh->primary_tag($pep_homol2->primary_tag);
    $newh->source_tag ($pep_homol2->source_tag);
    $newh->seqname    ($pep_homol2->seqname);
    $newh->score      ($pep_homol->score);                # This is wrong - this score is for the whole hit


    if (defined($pep_homol2->analysis)) {
	$newh->analysis($pep_homol2->analysis);
    } else {
	$self->throw("Can't find analysis for pep_homol " . $pep_homol2->seqname . " " . $pep_homol2->source_tag . "\n");
    }
    

    # Create the peptide homol
    my $peph  = new Bio::EnsEMBL::SeqFeature  (-start  => $h1,
					       -end    => $h2,
					       -strand => $pep_exon->strand * $pep_homol2->strand);

    $peph->primary_tag($pep_homol->primary_tag);
    $peph->source_tag ($pep_homol->source_tag);
    $peph->seqname    ($pep_homol->seqname);
    $peph->score      ($pep_homol->score);

    if (defined($pep_homol->analysis)) {
	$peph->analysis($pep_homol->analysis);
    } else {
	$self->throw("Can't find analysis for pep_homol " . $pep_homol->seqname . " " . $pep_homol->source_tag . "\n");
    }

    
    # Add the peptide homol to the dna homol
    my $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $newh,
					   -feature2 => $peph,
					   );
    
    $fp->analysis($pep_homol->analysis);
    return $fp;
}
   
=head2 _dnaHit2homol

  Title   : _dnaHit2homol
  Usage   : my @homols = $self->_dnaHit2homol
  Function: Converts a dnahit to the peptide to an array of genomic homols
  Returns : Array of Bio::SeqFeature::Homol;
  Args    : None

=cut

sub _dnaHit2homol {
    my ($self,$dnahit) = @_;
    my $debug = 0;

    my $pepaln  = $self->{_pepaln};
    my $pairaln = new Bio::EnsEMBL::Analysis::PepAlign;

    # As of 5/1/2000 the dna hit has to be on the forward
    # strand as I haven't written the conversion code 
    # for the other orientation so we swap here

    if ($dnahit->strand == -1) {
	$dnahit->strand(1);
	$dnahit->feature2->strand(-1);
    }

    # Also for the dna hit the dna homol must be the 'parent' homol
    # so we swap these as well
   

#    print("in dnahit2 Homol is " . $dnahit->homol_SeqFeature->start . " " . 
#	  $dnahit->homol_SeqFeature->end . " " . 
#	  $dnahit->homol_SeqFeature->seqname . "\n");
#    print("in dnahit2 Pep   is " . $dnahit->start . " " . $dnahit->end . " " . $dnahit->seqname ."\n");

    $dnahit->invert;
#    $dnahit = Bio::EnsEMBL::Analysis::MSPcrunch->swaphomols($dnahit);
    $pairaln->addFeaturePair($dnahit);

    my @pep     = $pepaln->eachFeaturePair;
    my @homols;
    
    foreach my $dna_homol ($pairaln->eachFeaturePair) {
	my $pep_homol = $dna_homol->feature2;

	# We now need to loop over each pepaln exon to see which
	# ones the homol overlaps

	# Find the first exon and process
	while (my $gen_exon = shift @pep) {
	    
	    my $pep_exon= $gen_exon->feature2;

	    if( $debug == 1 ) {
		print STDOUT "Debug: Before start: looking at exon in peptide coordinates of ",$pep_exon->start," to ",$pep_exon->end," dna ",$gen_exon->start," to ",$gen_exon->end,"\n";
	    }

	    if ($pep_homol->start >= $pep_exon->start && $pep_homol->start <= $pep_exon->end) {
		my $pstart     = $pep_homol->start;
		my $pend       = $pep_exon->end;
		#my $start_frac = $pep_exon->start_frac;
		my $start_frac = 1;
		my $end_frac;

		if ($pep_homol->end <= $pep_exon->end) {
		    $pend     = $pep_homol->end;
		    if( $pep_homol->end == $pep_exon->end ) {
			$end_frac = $pep_exon->end_frac;
		    } else {
			$end_frac = 3;
		    }
		} else {
		    $end_frac = $pep_exon->end_frac;
		}
    
		# We know whereabouts on the peptide we want the homol to lie
		# We need to convert these coords into genomic coords and 
		# also homol coords.
		# Hideous method call. yuck.

		if( $debug == 1 ) {
		    print STDOUT "Debug: Calling make_dna_homol with $pstart, $pend:$end_frac. Exon end is ",$pep_exon->end,":",$pep_exon->end_frac,"\n";
		}

		my $newh = $self->_make_dna_homol($pairaln,$gen_exon,$pep_exon,$dna_homol,$pep_homol,$pstart,$start_frac,$pend,$end_frac);

		if( $debug ) {
		    print STDOUT "Debug Start: Got homol with ",$newh->start,":",$newh->end," and ",$newh->feature2->start,":",$newh->feature2->end,"\n";
		}

		push(@homols,$newh);

		last;
	    }
	}

	# Process the rest of the exons
	while (my $gen_exon = shift @pep) {
	    # Check for overlap
	    my $pep_exon    = $gen_exon->feature2;
	    my $loop = 0;
	    my $pstart = $pep_exon->start;
	    my $pend   = $pep_exon->end;


	    
	    last if ($pep_homol->end < $pep_exon->start);
	    
	    my $start_frac = $pep_exon->start_frac;
	    my $end_frac;
	    
	    if ($pep_homol->end <= $pep_exon->end) {
		$pend     = $pep_homol->end;
		if( $pep_homol->end == $pep_exon->end ) {
		    $end_frac = $pep_exon->end_frac;
		} else {
		    $end_frac = 3;
		}

	    } else { 
		$end_frac = $pep_exon->end_frac;

	    }

	    if( $debug == 1 ) {
		print STDOUT "Debug: [$loop] After start. looking at exon in peptide coordinates of ",$pep_exon->start," to ",$pep_exon->end," dna ",$gen_exon->start," to ",$gen_exon->end,"\n";
	    }

		    
	    # We know whereabouts on the peptide we want the homol to lie
	    # We need to convert these coords into genomic coords and 
	    # also homol coords.
	    # /yuskety yuckety.

	    if( $debug == 1 ) {
		print STDOUT "Debug: Calling make_dna_homol with $pstart, $pend:$end_frac. Exon end is ",$pep_exon->end,":",$pep_exon->end_frac,"\n";
	    }
	    if( $end_frac != 1 && $end_frac != 2 && $end_frac != 3 ) {
		$self->warn("WE ARE IN TROUBLE. end frac is not real! [$end_frac]");
	    }

	    my $newh = $self->_make_dna_homol($pairaln,$gen_exon,$pep_exon,$dna_homol,$pep_homol,$pstart,$start_frac,$pend,$end_frac);
	    
	    if( $debug ) {
		print STDOUT "Debug [$loop]: Got homol with ",$newh->start,":",$newh->end," and ",$newh->feature2->start,":",$newh->feature2->end,"\n";
	    }

	    push(@homols,$newh);
	    $loop++;
	}
	
    }

    return @homols;
}

=head2 add_pepHit

  Title   : add_pepHit
  Usage   : $self->add_pepHit($homol);
  Function: Adds a protein hit to the predicted peptide
  Returns : Nothing
  Args    : Bio::SeqFeature::Homol

=cut

sub add_pepHit {
    my ($self,$pephit) = @_;

    $self->throw("Argument to GenscanPeptide->add_pepHit must be Bio::SeqFeature::FeaturePair") 
	unless $pephit->isa("Bio::SeqFeature::FeaturePair");

    push(@{$self->{_pepHits}},$pephit);
    
}

=head2 add_dnaHit

  Title   : add_dnaHit
  Usage   : $self->add_pepHit($pepaln);
  Function: Adds a dna hit to the predicted peptide
  Returns : Nothing
  Args    : Bio::EnsEMBL::Analysis::PepAlign

=cut

sub add_dnaHit {
    my ($self,$dnahit) = @_;

    $self->throw("Argument to GenscanPeptide->add_dnaHit must be Bio::SeqFeature::FeaturePair") 
	unless $dnahit->isa("Bio::SeqFeature::FeaturePair");

    push(@{$self->{_dnaHits}},$dnahit);
}

