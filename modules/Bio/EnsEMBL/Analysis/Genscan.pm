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

    my @genes    = $gs->each_Transcript;    # Returns an array of the predicted genes
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
use Bio::SeqIO;

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

  # The DNA Bio::Seq object 
  my $seq = shift(@args);
  $seq->isa("Bio::Seq") || $self->throw("No DNA sequence passed into GenScan analysis code.");

  
  $self->{_dna}  = $seq;
  
  
  # Stored data:
  # ------------
  # These are the predicted genes
  @{$self->{_transcripts}}    = ();

  # These are the peptides *as reported by genscan*
  # The translations of the genes are stored in
  # the gene objects
  @{$self->{_peptides}} = ();

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
      $self->id($seqname);
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
	     return;
	   };
	   
	   # End of genes section
	   last PARSE if /Predicted peptide/;
	   
	   # We have a line containing exon info	  
	   my @l = split;
	   
	   my ($n) = $l[0] =~ /^(\d+)\./;
	   $n--;
	   
	   # Get the right gene from the set
	   my $transcript = $self->_transcript( $n );
	   
	   # Is it an exon line?
	   if ( $l[1] =~ /^(Sngl|Init|Intr|Term)/ ) {
	     # Pass type,strand, start, stop, frame,phase to exons()
	     $self->_exons($transcript, @l[1,2,3,4,6,7] );
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

  #
  # EB
  # Attach all the exons to the sequence, using the attach_seq method
  # on Exons. Exons in EnsEMBL are "real" seqfeatures, and so can 
  # have DNA sequence attached to them
  #

  foreach my $t ( $self->each_Transcript() ) {
      foreach my $e ( $t->each_Exon() ) {
	  $e->attach_seq($self->{_dna});
      }
  }

  
  # Now deal with the predicted peptides - we need this
  # for finding the frame 
  
  my $tmp = <IN>; $tmp = <IN>;  # Read two blank lines
  
  my $in = Bio::SeqIO->new( '-fh' => \*IN, '-format' => 'Fasta');
  
  while (my $seq = $in->next_seq) {
    push(@{$self->{_peptides}},$seq);
  }
  
  # Set the phases of the exons now we have the peptides
  # This is a real hackeroony but the genscan ouput doesn't
  # give the right phase for the exons if there is no
  # initial exon in the prediction

  my $count = 0;
  
  foreach my $transcript ($self->each_Transcript) {
    
    # This sequence is what genscan thinks the translation is
    my $pep = $self->{_peptides}->[$count]->seq();

    # Sort the coordinates according to strand
    $transcript->sort();
    
    # Catch any exceptions where the phase can't be set for the gene

    if (defined($self->{_dna})) {
      $self->_set_exon_phases($transcript,$pep);
    }
    $count++;
  }

  close IN;
}


sub id {
  my ($self,$value)  = @_;

  if (defined($value)) {
    $self->{_id} = $value;
  }

  return $self->{_id};
}
# Takes a gene out from the array.

sub _remove_transcript {
  my ($self,$transcript) = @_;
  
  my $count = 0;
  
  if (defined($self->{_transcripts})) {
    foreach my $g ($self->each_Transcript) {
      if ($g == $transcript) {
	splice(@{$self->{_transcripts}},$count,1);
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
  my ($self,$tran,$pep) = @_;

  my $contig_dna = $self->{_dna};
  
  foreach my $exon ($tran->each_Exon) {
    my $seq   = $contig_dna->str($exon->start,$exon->end);

    if ($exon->strand == -1) {
      $seq =~ tr/ATCGatcg/TAGCtagc/;
      $seq = reverse($seq);
    }

    my @trans;
    $trans[0] = $exon->seq->translate();
    # this is because a phase one intron leaves us 2 base pairs, whereas a phase 2
    # intron leaves one base pair.
    $trans[1] = $exon->seq->translate('*','X',2);
    $trans[2] = $exon->seq->translate('*','X',1);

    #my @trans = $self->_translate($seq);

#    for (my $i = 0; $i < 3; $i++) {
#      print("PHASE $i : " . $trans[$i]->seq() . "\n");
#    }

#    print("Genscan phase is " . $exon->phase . "\n");
#    print("Genscan frame is " . $exon->frame . "\n");

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
      #print STDERR "new phase is $phase\n";
      $exon->phase($phase);
    } else {
      $self->warn("Can not find frame for exon. Sequences do not match\n");
    }

  }

  return;

  # EB - no longer need this. Cope somewhere else.

  # Genscan DNA coordinates include a stop codon at the end of the terminal exon.  We 
  # need to remove this and change the coordinates
  $tran->sort();

  my @exons = $tran->each_Exon;
  my $ex    = $exons[$#exons];

  my $exstr = $contig_dna->str($ex->start,$ex->end);

  if ($ex->strand == -1) {
    $exstr =~ tr/ATCGatcg/TAGCtagc/;
    $exstr = reverse($exstr);
  }
  
  my @pep = $self->_translate($exstr);
  my $pep_p = $pep[$ex->phase]->seq;

  if (substr($pep_p,-1) eq "*") {
    if ($ex->strand == -1) {
      $ex->start($ex->start + 3);
    } else {
      $ex->end($ex->end()-3);
    }
  }
}


sub _translate {
  my ($self,$seq) = @_;

  my @trans;

  for (my $i = 0; $i < 3; $i++) {
    my $tmp = $seq;
       $tmp = substr($seq,$i);

    my $bs = new Bio::Seq(-seq => $tmp);
    my $tr = $bs->translate();
    
    push(@trans,$tr);
  }
  return @trans;
}

# Returns the nth Transcript object from the
# genes array.  Returns a new Transcript object
# if the object doesn't exist

sub _transcript { 
  my ($self,$n) = @_;

  if ($#{$self->{_transcripts}} >= $n) {
    return $self->{_transcripts}[$n];
  } else {
    my $i;

    for ($i = $#{$self->{_transcripts}} +1; $i <= $n; $i++){
      $self->{_transcripts}[$i] = Bio::EnsEMBL::Transcript->new();
    }

    return $self->{_transcripts}[$n];

  }

}

# Constructs an exon hash and attaches it
# to the parent gene exon array ref.

sub _exons {
  my ($self,$tran,$type,$strand,$start,$stop,$frame,$phase) = @_;

  # Create the exon object
  my $exon = new Bio::EnsEMBL::Exon($start,$stop,$strand);

  # Set the other variables
  $exon->type     ($type);
  $exon->phase    ($phase);  # This will get overwritten if the dna seq. is input
  $exon->frame    ($frame);
  $exon->end_phase();		

  # Finally add the exon to the gene
  $tran->add_Exon ($exon);

}


=head2 each_Transcript

  Title   : each_Transcript
  Usage   : @genes = $gs->each_Transcript
 Function: Returns an array of predicted genes
  Returns : Bio::SeqFeature
  Args    : none

=cut

sub each_Transcript {
  my ($self) = @_;
  
  return (@{$self->{_transcripts}});
  
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


sub toSQLfeatureset {

  my ($self) = @_;

  my @sqllines;

  my $date   = `date '+%Y-%m-%d'`; chomp $date;;
  my $gcount = 0;
  
  foreach my $tran ($self->each_Transcript) {
    
    # Create gene sql
    # the geneid is initially of the form
    # dJ401P4.00741.GENSCAN.1.2
    # i.e. 3rd exon of the 2nd gene (counting from 0)
    
    my $tran_id = $tran->id;  #!!!!!
    
    $tran_id =~ s/(.*)\..*$/$1/;

    # Extract clone and contig info;
    my $clone      = $tran_id;
    my $contig     = $tran_id;
    my $featureset = $tran_id;

    $clone  =~ s/^(.*?)\..*/$1/;
    $contig =~ s/^(.*?\..*?)\..*/$1/;

    my $tmp;
    
    foreach my $gf ($tran->each_Exon()) {
      # Exon sql
      
      $tmp = "insert into feature(id,contig,start,end,strand,featureset) values(\'" . 
           #$gf->seq->id()   .   "\',\'$contig\'," .
           "\',\'$contig\'," .
           $gf->start     .   "," .
           $gf->end       .   ",\'" . 
           $gf->strand    .   "\',\'" .
           $featureset    .   "\');\n";

      push(@sqllines,$tmp);

      # featureset sql
      $tmp = "insert into featureset(feature,id) values(\'" . $contig . 
        "\',\'$featureset\');\n";
      push(@sqllines,$tmp);
    }
    $gcount++;
  }
  return @sqllines;
}

    
1;
