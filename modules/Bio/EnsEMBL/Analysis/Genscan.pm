#
# BioPerl module for Genscan
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

Bio::EnsEMBL::Analysis::Genscan - Parses genscan gene prediction output

=head1 SYNOPSIS

    my $gs = new Bio::EnsEMBL::Analysis::Genscan($file);      # $file is filename
or
    my $gs = new Bio::EnsEMBL::Analysis::Genscan($file,$dna); # $dna is Bio::Seq

Extracting the data

    my @genes    = $gs->each_Gene;    # Returns an array of the predicted genes
    my @peptides = $gs->each_Peptide; # Returns an array of the genscan peptides
    my $dna      = $gs->dna_seq       # Returns a Bio::Seq containging the DNA.

    Note: The genscan predicted peptides do NOT necessarily agree with the
          peptides in the genes array.  They are sometimes offset by 1.

=head1 DESCRIPTION

Genscan object. Parses the output file from the genscan gene prediction program.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::Genscan;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object;

use Bio::Root::Object;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;

# Inherits from the base bioperl object
@ISA = qw(Bio::Root::Object);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  
  my $make = $self->SUPER::_initialize;

  # Input variables
  # ---------------
  # Genscan filename
  # Bio::Seq DNA sequence

  die ("No genscan input file") unless $#args > 0;
  my $file      = shift(@args);

  # The DNA Bio::Seq object - optional
  $self->{dna}  = shift(@args);
  
  # Stored data:
  # ------------
  # These are the predicted genes
  @{$self->{genes}}    = ();

  # These are the peptides *as reported by genscan*
  # The translations of the genes are stored in
  # the gene objects
  @{$self->{peptides}} = ();

  # Now try and get some genes out
  $self->_parse($file);

  return $self; # success - we hope!
}



# Parses the genscan file.
# Fills up the genes array and
# also calculates the phases for the exons
# if a dna object exists

sub _parse {
  my ($self,$file) = @_;

  my $seqname;

  open(IN,"<$file") || die "Couldn't open genscan file $file\n";
  
  while(<IN>) {
    if (/^Sequence +(\S+)/) {
      $seqname = $1;
      close(IN);
    }
  } 
  
  open(IN,"<$file");
  
  my $No_Gene_Flag;
  
   PARSE: while(<IN>) {
       # Header gives Genscan version
       my $version;
       if (/^GENSCAN\s*(\S+)/o) {
	 $version = $1;
       }
       
       # Last line before predictions contains nothing
       # but spaces and dashes
       if (/^\s*-[-\s]+$/) {
	 
	 while (<IN>) {
	   
	   next if /^$/;
	   
	   # If sequence is too short;
	   if (m|NO EXONS/GENES PREDICTED IN SEQUENCE|) {
	     $No_Gene_Flag = 1;
	     last PARSE;
	   };
	   
	   # End of genes section
	   last PARSE if /Predicted peptide/;
	   
	   # We have a line containing exon info	  
	   my @l = split;
	   
	   my ($n) = $l[0] =~ /^(\d+)\./;
	   $n--;
	   
	   # Get the right gene from the set
	   my $gene = $self->_gene( $n );
	   
	   # Is it an exon line?
	   if ( $l[1] =~ /^(Sngl|Init|Intr|Term)/ ) {
	     # Pass type,strand, start, stop, phase to exons()
	     $self->_exons($gene, @l[1,2,3,4,7] );
	   }
	   
	   # or a Promoter?
	   elsif ( $l[1] =~ /^Prom/ ) {
	     #	    $gene->_prom( $l[3] );
	   }
	   # or a Poly-A?
	   elsif ( $l[1] =~ /^PlyA/ ) {
	     #	    $gene->_poly( $l[3] );
	   }
	   # Unknown line type
	   else {
	     chomp;
	     die "Unknown line type ('$_')";
	   }
	 }
       }
     }
  
  # Now deal with the predicted peptides - we need this
  # for finding the frame 
  
  my $tmp = <IN>; $tmp = <IN>;  # Read two blank lines
  
  my $in = Bio::SeqIO->new(-fh => \*IN, -format => 'Fasta');
  my $count = 0;
  
  while (my $seq = $in->next_seq) {
    push(@{$self->{peptides}},$seq);
  }
  
  # Set the phases of the exons now we have the peptides
  # This is a real hackeroony but the genscan ouput doesn't
  # give the right phase for the exons if there is no
  # initial exon in the prediction

  my $count = 0;
  
  foreach my $gene ($self->each_Gene) {
    
    # This sequence is what genscan thinks the translation is
    my $pep = $self->{peptides}->[$count]->seq();

    # Sort the coordinates according to strand
    $gene->sort();
    
    # Catch any exceptions where the phase can't be set for the gene
    if (defined($self->{dna})) {
      $self->_set_exon_phases($gene,$pep);
    }
    $count++;
  }

  close IN;
}

# Takes a gene out from the array.

sub _remove_gene {
  my ($self,$gene) = @_;
  
  my $count = 0;
  
  if (defined($self->{genes})) {
    foreach my $g ($self->each_Gene) {
      if ($g == $gene) {
	splice(@{$self->{genes}},$count,1);
      }
      $count++;
    }
  }
}

# This is a nightmare.  As far as I can see Genscan only reports
# the correct translation phase for each exon _if_ it has predicted a 
# full gene. i.e. the gene starts with a promoter or an Initial exon 
# (Initial exons always start with phase 0).  If there is no promoter or
# Initial exon then the frame and phase information is meaningless.
# N.b. I could be wrong here but I can't see any obvious pattern
#
# We calculate the phase by *ahem* translating the DNA sequence in
# all three frames and comparing the string to the full peptide sequence.

sub _set_exon_phases {
  my ($self,$gene,$pep) = @_;

  foreach my $exon ($gene->each_Exon) {
    my $seq = $exon->dna_seq->seq();
    my @trans = $exon->pep_seq;

    my $i = 0;
    my $phase;

    # Loop over all frames 0,1,2
    for ($i=0; $i < 3; $i++) {

      # If we have a stop codon at the end of the translation
      # chop it off before comparing
      my $tmp = $trans[$i]->seq();
      
      if (substr($tmp,-1) eq "*") {
	$tmp = substr($tmp,0,-1);
      }

      # Compare strings to see if the exon peptide is contained in 
      # the full sequence
      if ((my $pos = index($pep,$tmp)) >= 0) {
	$phase = $i;
      }
    }

    # Seet phase if poss.  If no phase is found the input DNA is
    # probably wrong.
    if (defined($phase)) {
      $exon->phase($phase);
    } else {
      $self->warn("Can not find frame for exon. Sequences do not match\n");
    }

  }


  # Genscan DNA coordinates include a stop codon at the end of the terminal exon.  We 
  # need to remove this and change the coordinates

  my @exons = $gene->each_Exon;
  my $ex = $exons[$#exons];

  my @pep = $ex->_translate();
  my $pep = $pep[$ex->phase]->seq;

  if (substr($pep,-1) eq "*") {
    $ex->end($ex->end()-3);
    my $seq = $ex->dna_seq->seq();
    $ex->dna_seq(new Bio::Seq(-seq => substr($seq,0,-3)));
  }
}

# Returns the nth Transcript object from the
# genes array.  Returns a new Transcript object
# if the object doesn't exist

sub _gene { 
  my ($self,$n) = @_;

  if ($#{$self->{genes}} >= $n) {
    return $self->{genes}[$n];
  } else {
    my $i;

    for ($i = $#{$self->{genes}} +1; $i <= $n; $i++){
      $self->{genes}[$i] = Bio::EnsEMBL::Transcript->new();
    }

    return $self->{genes}[$n];

  }

}

# Constructs an exon hash and attaches it
# to the parent gene exon array ref.

sub _exons {
  my ($self,$gene,$type,$strand,$start,$stop,$phase) = @_;

  # Create the exon object
  my $exon = new Bio::EnsEMBL::Exon($start,$stop,$strand);

  # Create the sequence object for the exon
  if (defined($self->{dna})) {
    my $seq = $self->{dna}->seq($exon->start,$exon->end);
    my $newseq = Bio::Seq->new(-seq => $seq);

    # Reverse complement if necessary
    if ($strand eq '-') {
      $newseq = $newseq->revcom();
    }
    
    $exon->dna_seq($newseq);
  }

  # Set the other variables
  $exon->type     ($type);
  $exon->phase    ($phase);  # This will get overwritten if the dna seq. is input
  $exon->end_phase();		

  # Finally add the exon to the gene
  $gene->add_Exon ($exon);

}


=head2 each_Gene

  Title   : each_Gene
  Usage   : @genes = $gs->each_Gene
 Function: Returns an array of predicted genes
  Returns : Bio::SeqFeature
  Args    : none

=cut

sub each_Gene {
  my ($self) = @_;

  return (@{$self->{genes}});
  
}

=head2 each_Peptide

  Title   : each_Peptide
  Usage   : @peps = $gs->each_Peptide
  Function: Returns an array of predicted peptides
  Returns : Bio::SeqFeature
  Args    : none

=cut

sub each_Peptide {
  my ($self) = @_;

  return (@{$self->{peptides}});
  
}

=head2 dna_seq

  Title   : dna_seq
  Usage   : $dna = $gs->dna_seq
  Function: Returns the genomic dna
  Returns : Bio::Seq
  Args    : none

=cut

sub dna_seq {
  my ($self) = @_;

  if (defined($self->{dna})) {
    return ($self->{dna});
  }
}
    
1;
