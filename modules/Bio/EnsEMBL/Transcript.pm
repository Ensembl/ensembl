
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

Transcript - gene transcript object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Contains details of coordinates of all exons that make
up a gene transcript.

Creation:
   
     my $tran = new Bio::EnsEMBL::Transcript();
     my $tran = new Bio::EnsEMBL::Transcript(@exons);

Manipulation:

     my @exons = $tran->each_Exon         # Returns an array of Exon objects
     my $pep   = $tran->translate()       # Returns the peptide translation of the exons as a Bio::Seq
     
     $tran->sort()                        # Sorts exons into order (forward for + strand, reverse fo - strand)

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Transcript;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Object;
use Bio::EnsEMBL::Exon;


@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  $self->{'_trans_exon_array'} = [];
  my $make = $self->SUPER::_initialize;

  # set stuff in self from @args
  foreach my $a (@args) {
    $self->add_Exon($a);
  }
  $self->is_partial(0);

  return $self; # success - we hope!
}

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'id'} = $value;
    }
    return $self->{'id'};

}

=head2 add_Exon

 Title   : add_Exon
 Usage   : $trans->add_Exon($exon)
 Returns : Nothing
 Args    :


=cut

sub add_Exon{
   my ($self,$exon) = @_;

   #yup - we are going to be picky here...

   if( ! $exon->isa("Bio::EnsEMBL::Exon") ) {
       $self->throw("$exon is not a Bio::EnsEMBL::Exon!");
   }

   # at the moment, use the SeqFeature sub hash. But in the future,
   # possibly do something better?

   push(@{$self->{'_trans_exon_array'}},$exon);
   
}

=head2 each_Exon

 Title   : each_Exon
 Usage   : foreach $exon ( $trans->each_Exon)
 Function: Returns an array of exons in the transcript
 Example : my @exons = $tr->each_Exon
 Returns : An array of exon objects
 Args    : none


=cut

sub each_Exon{
   my ($self) = @_;
   my @sub;
   my @ret;

   return @{$self->{'_trans_exon_array'}};
}


=head2 flush_Exon

 Title   : flush_Exon
 Usage   : Removes all Exons from the array.
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub flush_Exon{
   my ($self,@args) = @_;

   $self->{'_trans_exon_array'} = [];
}


=head2 split_Transcript_to_Partial

 Title   : split_Transcript_to_Partial
 Usage   : @trans = $trans->split_Transcript_to_Partial
 Function: splits a transcript with potential non spliceable
           exons into a set of partial transcripts
 Example :
 Returns : an array of Bio::EnsEMBL::Transcript objects
 Args    :


=cut

sub split_Transcript_to_Partial{
   my ($self,@args) = @_;

   my @exons = $self->each_Exon;

   # one exon genes - easy to handle.
   if( $#exons == 0 ) {
       return $self;
   }

   my $l = $#exons;
   my $prev = shift @exons;


   my $t;
   my @out;


   TRANSCRIPT :
   while ( $#exons >= 0 ) {


       # make a new transcript, add the old exon
       $t = $self->new();
       $t->id($self->id);
       
       $t->add_Exon($prev);
       $t->is_partial(1);
       push(@out,$t);

       while( my $exon = shift @exons ) {
	   if( $exon->phase == $prev->end_phase ) {
	       # add it
	       $t->add_Exon($exon);
	       $prev = $exon;
	   } else {
	       $prev = $exon;
	       if( $#exons < 0 ) {
		   # this was the last exon!
		   $t = $self->new();
		   $t->id($self->id);
		   
		   $t->add_Exon($prev);
		   $t->is_partial(1);
		   push(@out,$t);
		   last TRANSCRIPT;
	       } else {
		   next TRANSCRIPT;
	       }
	   }
       }
   }

   if( $#out == 0 ) {
       $t->is_partial(0);
   }

   return @out;
}


=head2 translate

 Title   : translate
 Usage   : $pep = $feat->translate()
 Function: returns the peptide translation of the gene - in the correct phase
 Returns : Bio::Seq
 Args    : none

=cut

sub translate {
  my ($self) = @_;

  my $debug;
  
  my @trans = $self-> split_Transcript_to_Partial();

  if( $#trans == -1 ) {
      $self->throw("Bad internal error - split a transcript to zero transcripts! Doh!");
  }


  my $seqstr;
  foreach my $ptrans ( @trans ) {
      my $tseq = $ptrans->_translate_coherent($debug);
      # to be consistent with our EMBL dumping, we need a double X here.
      if( defined $seqstr ) { $seqstr .= 'XX'; } 
      $seqstr .= $tseq->str;
  }
  
  $seqstr =~ s/\*$//g;

  my $trans_seq = Bio::Seq->new( -seq => $seqstr , -id => $self->id() ) ;


  return $trans_seq;
}

=head2 dna_seq

  Title   : dna_seq
  Usage   : $dna = $feat->dna_seq
  Function: Returns the dna sequence of the gene, ie the mRNA
            sequence
  Returns : Bio::Seq
  Args    : none

=cut

sub dna_seq {
  my ($self) = @_;

  my $mrna = "";
  my $strand = $self->{_trans_exon_array}[0]->strand;

  my $prev = undef;
  foreach my $exon ($self->each_Exon) {

    # the seq call automatically truncates to the correct 
    # coordinates (handily) in SeqFeature

    my $tmp = $exon->seq->str();

    # we now have to figure out if the phase is compatible. If it
    # is not, we need to add some stuff in...

    if( $prev ) {
	if( $prev->end_phase != $exon->phase ) {
	    if( $prev->end_phase == 0 ) {
		if( $exon->phase == 1 ) {
		    $mrna .= 2 x "N";
		}

		if( $exon->phase == 2 ) {
		    $mrna .= 1 x "N";
		}
	    } elsif ( $prev->end_phase == 1 ) {
		if( $exon->phase == 0 ) {
		    $mrna .= 2 x "N";
		}
		
		if( $exon->phase == 2 ) {
		    $mrna .= 1 x "N";
		}
	    } elsif ( $prev->end_phase == 2 ) {
		if( $exon->phase == 0 ) {
		    $mrna .= 1 x "N";
		}
		
		if( $exon->phase == 1 ) {
		    $mrna .= 2 x "N";
		}
	    } else {
		$self->warn("Impossible phases in calculating fixing stuff");
	    }
	}
    } # end of if previous is there

    
    $mrna  .= $tmp;
  }


  my $seq = new Bio::Seq(-seq => $mrna);

  return $seq;
}

=head2 contig_dna

  Title   : contig_dna
  Usage   : $tran->contig_dna($dna);
   Function: Sets the dna sequence of the contig
  Returns : Bio::Seq
  Args    : Bio::Seq

=cut

sub contig_dna {
  my ($self,$dna) = @_;

  if (defined($dna)) {
    $self->{_contig_dna} = $dna;
  }

  return $self->{_contig_dna};
}

=head2 sort

 Title   : sort
 Usage   : $feat->sort()
 Function: Sorts the exon features by start coordinate
           Sorts forward for forward strand and reverse for reverse strand
 Returns : none
 Args    : none

=cut

sub sort {
  my $self = shift;

  # Fetch all the features
  my @exons = $self->each_Exon();

  # Empty the feature table
  $self->flush_Exon();

  # Now sort the exons and put back in the feature table
  my $strand = $exons[0]->strand;

  if ($strand == 1) {
    @exons = sort { $a->start <=> $b->start } @exons;
  } elsif ($strand == -1) {
    @exons = sort { $b->start <=> $a->start } @exons;
  }

  foreach my $e (@exons) {
    $self->add_Exon($e);
  }
}

=head2 _translate_coherent

 Title   : _translate_coherent
 Usage   : <internal function> translates a coherent transcript.
           Uncoherent transcripts need to be broken up with
           split to partial first.
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _translate_coherent{
   my ($self,$debug) = @_;
   
   my $prev;
   my $tstr;
#   my $debug;

   $self->sort();
   my @exons = $self->each_Exon;
   my $exon_start = $exons[0];

   foreach my $exon ( @exons ) {
      # if( $exon->id eq 'HE000030314' ) {
#	   print STDERR "setting debug to 1\n";
	#   $debug = 1;
#       }

       # trim down start ends on the basis of phase.
       if( $prev && $prev->end_phase != $exon->phase ) {
	   $self->throw("Called coherent translate but exon phases don't match. Yuk!");
       }

       # warn about non DNA passed in. 

       if( $exon->entire_seq()->type ne 'Dna' ) {
	   #$self->warn("Error. Whoever implemented this databases did not set type to Dna. Setting now!");
	   $exon->entire_seq()->type('Dna');
       }

       my $seq = $exon->seq();
       my $str = $seq->str();
       

       if( length $str == 0 ) {
	   $self->throw("Bad internal error - got a 0 length rstring...");
       }

       $tstr .= $str;
   }

   if( $exon_start->phase == 1 ) {
       $tstr = substr $tstr, 2;
   } elsif ( $exon_start->phase == 2 ) {
       $tstr = substr $tstr, 1;
   } 

   if ( $debug ) {
       print STDERR "Tstr is $tstr\n";
   }

   # phase 0 - no need.


   my $temp_seq = Bio::Seq->new( -seq => $tstr , -id => 'temp', -type => 'Dna' );
   my $trans_seq = $temp_seq->translate();

   return $temp_seq->translate();
}


sub exon_dna {
  my ($self,$exon) = @_;

  my $tmpseq = $self->contig_dna->str($exon->start,$exon->end);
  
  if ($exon->strand == -1) {
    $tmpseq =~ tr/ATGCatgc/TACGtacg/;
    $tmpseq = reverse($tmpseq);
  }
  return new Bio::Seq(-seq => $tmpseq);
}

  
sub pep_coords {
  my $self = shift;

  # for mapping the peptide coords back onto the dna sequence
  # it would be handy to have a list of the peptide start end coords
  # for each exon
  
  my @starts;
  my @ends;
  
  my $fullpep = $self->translate()->seq;

  foreach my $ex ($self->each_Exon) {
    
    my @tmp = $self->translate_exon($ex);
    my $pep = $tmp[$ex->phase]->seq;
    
    my $start = index($fullpep,$pep) + 1;
    my $end = $start + length($pep) - 1;
    
    push(@starts,$start);
    push(@ends,$end);
    
  }

  return \@starts,\@ends;
}




sub find_coord {
  my ($self,$coord,$type) = @_;
 
  my $count = 0;
  my @exons = $self->each_Exon;
  my $end   = $#exons;
  my $dna;

  my ($starts,$ends) = $self->pep_coords;
  my $strand = $exons[0]->strand;

  # $starts and $ends are array refs containing the _peptide_ coordinates
  # of each exon. We may have 1 missing residue that spans an intron.
  # We ignore these.

  if ($strand == 1) {
    foreach my $ex ($self->each_Exon) {
      
      if ($coord >= $starts->[$count] && $coord <= $ends->[$count]) {
	my $dna   = $ex->start + $ex->phase;
	my $nopep = $coord - $starts->[$count];
	
	$dna += 3 * $nopep;

	if ($type eq "end") {
	  $dna += 2;
	}
	
	return $dna;
	
      } elsif ($count < $end) {
	my $endpep = $ends->[$count]+1;
	if ($endpep == $coord) {

	  my $dna;

	  if ($type eq "end") {
	    my $end_phase = $ex->end_phase;
	    $dna = $ex->end - 3 + $end_phase;
	  } else {
	    $dna = $exons[$count+1]->start + $exons[$count+1]->phase;
	  }
	  return $dna;
	}
      }
      $count++;
    }
  } else {

    foreach my $ex ($self->each_Exon) {
      
      if ($coord >= $starts->[$count] && $coord <= $ends->[$count]) {
	
	my $dna   = $ex->end - $ex->phase;
	my $nopep = $coord - $starts->[$count];

	$dna -= 3*$nopep;

	if ($type eq "end") {
	  $dna -= 2;
	}
	
	return $dna;
	
      } elsif ($count < $end) {
	my $endpep = $ends->[$count]+1;

	if ($endpep == $coord) {
	  my $dna;

	  if ($type eq "end") {
	    my $end_phase = $ex->end_phase;
	    $dna = $ex->start + 3 - $end_phase;
	  } else {
	    $dna = $exons[$count+1]->end - $exons[$count+1]->phase;
	  }
	  return $dna;
	}
      }
      $count++;
    } 
  }
}

=head2 is_partial

 Title   : is_partial
 Usage   : $obj->is_partial($newval)
 Function: 
 Example : 
 Returns : value of is_partial
 Args    : newvalue (optional)


=cut

sub is_partial{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'is_partial'} = $value;
    }
    return $obj->{'is_partial'};

}


1;





