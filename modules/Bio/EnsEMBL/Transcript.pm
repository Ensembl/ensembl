
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
use Bio::EnsEMBL::Translation;
use Bio::DBLinkContainerI;


@ISA = qw(Bio::Root::Object Bio::DBLinkContainerI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  $self->{'_trans_exon_array'} = [];
  $self->{'_db_link'} = [];
  my $make = $self->SUPER::_initialize;

  # set stuff in self from @args
  foreach my $a (@args) {
    $self->add_Exon($a);
  }
  $self->is_partial(0);

  return $self; # success - we hope!
}


=head2 each_DBLink

 Title   : each_DBLink
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_DBLink {
   my ($self,@args) = @_;

   return @{$self->{'_db_link'}}
}

=head2 add_DBLink

 Title   : add_DBLink
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_DBLink{
   my ($self,$value) = @_;

   if( !defined $value || !ref $value || ! $value->isa('Bio::Annotation::DBLink') ) {
       $self->throw("This [$value] is not a DBLink");
   }

   push(@{$self->{'_db_link'}},$value);
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

=head2 version

 Title   : version
 Usage   : $obj->version($newval)
 Function: 
 Returns : value of version
 Args    : newvalue (optional)


=cut

sub version{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'version'} = $value;
    }
    return $obj->{'version'};

}


=head2 translation

 Title   : translation
 Usage   : $obj->translation($newval)
 Function: 
 Returns : value of translation
 Args    : newvalue (optional)


=cut

sub translation {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      if( ! ref $value || !$value->isa('Bio::EnsEMBL::Translation') ) {
	  $obj->throw("This [$value] is not a translation");
      }
      $obj->{'translation'} = $value;
    }
    return $obj->{'translation'};

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
           in order, ie the first exon is the 5' most exon
           in the transcript (the one closest to the 5' UTR).
 Example : my @exons = $tr->each_Exon
 Returns : An array of exon objects
 Args    : none


=cut

sub each_Exon{
   my ($self) = @_;

   return @{$self->{'_trans_exon_array'}};
}

=head2 number

 Title   : number
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub number{
   my ($self,@args) = @_;

   return scalar(@{$self->{'_trans_exon_array'}});   
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

=head2 first_exon

 Title   : first_exon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub first_exon{
   my ($self,@args) = @_;
   my @temp = @{$self->{'_trans_exon_array'}};

   return shift @temp;
}

=head2 last_exon

 Title   : last_exon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub last_exon{
   my ($self) = @_;
   my @temp = @{$self->{'_trans_exon_array'}};

   return pop @temp;
}

=head2 translatable_exons

 Title   : translatable_exons
 Usage   : @exons = $transcript->translateable_exons
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub translateable_exons{
   my ($self) = @_;
   my (@out);

   if( ! defined $self->translation ) {
       $self->throw("Attempting to split Transcript on translation, but not there...");
   }
   
   my @exons = $self->each_Exon;

   # one exon genes - easy to handle.
   if( $#exons == 0 ) {
       my $exon = shift @exons;
       if( $self->translation->start_exon_id ne $exon->id || $self->translation->end_exon_id ne $exon->id ) {
	   $self->throw("Single Exon transcript, but with start or stop outside of that exon");
       }

       my $retexon = new Bio::EnsEMBL::Exon;

       $retexon->contig_id ($exon->contig_id);
       $retexon->clone_id  ($exon->clone_id);
       $retexon->strand    ($exon->strand);
       $retexon->phase     ($exon->phase);         # MC exon phase can be 0,1,2
       $retexon->attach_seq($exon->entire_seq);
       $retexon->id($exon->id());
       
       if( $exon->strand == 1 ) {
	   $retexon->start($self->translation->start());
	   $retexon->end  ($self->translation->end());
       } else {
	   $retexon->end  ($self->translation->start());
	   $retexon->start($self->translation->end());
       }
       
       return $retexon;
   }

   while( my $exon = shift @exons ) {
       if( $exon->id eq $self->translation->start_exon_id() ) {
	   #print STDERR "New start exon " . $exon->id . "\n";
	   my $stexon = new Bio::EnsEMBL::Exon;

	   $stexon->contig_id ($exon->contig_id);
	   $stexon->clone_id  ($exon->clone_id);
	   $stexon->strand    ($exon->strand);
	   $stexon->phase     ($exon->phase);             # MC exon phase can be 0,1,2

	   $stexon->attach_seq($exon->entire_seq());

#	   print (STDERR "Exon entire seq is " . ref($stexon->entire_seq) . "\n");

	   $stexon->id($exon->id());
	   
	   if( $exon->strand == 1 ) {
	       # paranoid...
	       if( $self->translation->start < $exon->start || $self->translation->start > $exon->end ) {
		   $self->throw("For start exon " . $exon->id . 
				" translation start not within exon bounds. Start " 
				. $self->translation->start . 
				"Exon "  . $exon->start . ":" . 
				           $exon->end."\n");
	       }

	       $stexon->start($self->translation->start());
	       $stexon->phase(0);             # MC translation is always phase 0.
	       $stexon->end($exon->end);

	       #print STDERR "Setting start end to " . $stexon->start   . "\t" . $stexon->end  ."\n";
	       #print (STDERR "Exon sequence is "    . $stexon->seq->seq . "\n");	       
	   } else {
	       if( $self->translation->start < $exon->start || $self->translation->start > $exon->end ) {
		   $self->throw("For start exon ".$exon->id." translation start not within exon bounds. Start ". 
                        $self->translation->start . " Exon " .$exon->start." : ".$exon->end."\n");
	       }
	       $stexon->end($self->translation->start());
	       $stexon->start($exon->start);
	       $stexon->phase(0);             # MC translation is always phase 0.
	   }
	   push(@out,$stexon);
	   last;
       }
   }

   my $exon;
   while( $exon = shift @exons ) {
       if( $exon->id eq $self->translation->end_exon_id()) {
	   
	   my $endexon = new Bio::EnsEMBL::Exon;
	   
	   $endexon->contig_id($exon->contig_id);
	    $endexon->clone_id($exon->clone_id);
	    $endexon->strand($exon->strand);
	    $endexon->phase($exon->phase);
	    $endexon->attach_seq($exon->entire_seq());
	    $endexon->id($exon->id);

	    if( $exon->strand == 1 ) {
		if( $self->translation->end() > $exon->end ) {
		    $self->throw("Bad news. Attempting to say that this translation is inside this exon, but outside".$exon->id." ".$exon->end()." ".$self->translation->end()."\n");
		}

		$endexon->start($exon->start());
		$endexon->end  ($self->translation->end());

	    } else {
		if( $self->translation->end() < $exon->start ) {
		    $self->throw("Bad news. Attempting to say that this translation is inside this exon (reversed), but outside".$exon->id." ".$exon->start()." ".$self->translation->end()."\n");
		}
		$endexon->start($self->translation->end());
		$endexon->end($exon->end);
	    }
	    push(@out,$endexon);
	    last;
	} else {
	    push(@out,$exon);
	}
    }

    if( !defined $exon || $exon->id ne $self->translation->end_exon_id()) {
	$self->throw("Unable to find end translation exon");
    }


    return @out;
 }


=head2 split_Transcript_to_Partial

  Title   : split_Transcript_to_Partial
  Usage   : @trans = $trans->split_Transcript_to_Partial
  Function: splits a transcript with potential non spliceable
	    exons into a set of partial transcripts, from start/end
	    of the Translation if the second argument is true
  Example :
  Returns : an array of Bio::EnsEMBL::Transcript objects
  Args    :


=cut

 sub split_Transcript_to_Partial{
    my ($self,$on_translate) = @_;


    if( $on_translate == 1 && ! defined $self->translation ) {
	$self->throw("Attempting to split Transcript on translation, but not there...");
    }
    my @exons;
    if( $on_translate == 1 ) {
	@exons = $self->translateable_exons();
    } else {
	@exons = $self->each_Exon;
    }

    # one exon genes - easy to handle.
    if( $#exons == 0 ) {
	return $self;
    }

    # find the start exon;

    my $prev;


    $prev = shift @exons;

    my $l = $#exons;

    my $t;
    my @out;


    TRANSCRIPT :
    while ( $#exons >= 0 ) {
	# MC - we are allowing split phases within a contig now.

	# make a new transcript, add the old exon
	$t = $self->new();
	$t->id($self->id);

	$t->add_Exon($prev);
	$t->is_partial(1);

	push(@out,$t);

	while( my $exon = shift @exons ) {
	    if( $exon->contig_id eq $prev->contig_id && $exon->phase == $prev->end_phase) {
		# add it
		$t->add_Exon($exon);
		$prev = $exon;
	    } else {
 #	       if( $exon->contig_id eq $prev->contig_id ) {
 #		   $self->warn("Got two exons, ".$prev->id. "[".$prev->contig_id."]".$prev->start."/".$prev->phase .":". $exon->id ."[".$exon->contig_id."]: same contigs but incompatible phases!");
 #	       }

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


  if ( ! defined $self->translation ) {
      $self->throw("You have to have a translation now to make a translation. Doh!");
  }

  my @trans = $self-> split_Transcript_to_Partial(1);

  if( $#trans == -1 ) {
      $self->throw("Bad internal error - split a transcript to zero transcripts! Doh!");
  }

  

  my $seqstr;
  my $prevtrans;
  my $seen_start =0;

  foreach my $ptrans ( @trans ) {
      
      my $tseq = $ptrans->_translate_coherent($debug);
     
      # to be consistent with our EMBL dumping, we need to make the actual join
      # with fill-in N's and then pass into translate so that ambiguous codons
      # which still have a translation happen! This has to be the weirdest piece
      # of manipulation I have done in a long time. What we do is take the last
      # exon of the old transcript and the first exon of the current transcript,
      # fill both sides in so that they make nice reading frames and then translate
      # the 6 bases which this makes. Phase 0 exons are filled in by 3 to make sure
      # that there is some filler.
      
      if( defined $prevtrans ) {
	  my $last_exon  = $prevtrans->last_exon();
	  my $first_exon = $ptrans   ->first_exon();
	  my $filler;

	  # last exon
	  if( $last_exon->end_phase != 0 ) {
	      $filler = substr($last_exon->seq->seq,$last_exon->seq->length-$last_exon->end_phase);
	  } 
	  $filler .= 'N' x (3 - $last_exon->end_phase);

	  # do first exon now.

	  if( $first_exon->phase != 0 ) {
	      $filler .= 'N' x $first_exon->phase;
	  } else {
	      $filler .= 'NNN';
	  }
	  $filler .= substr($first_exon->seq->seq,0,(3-$first_exon->phase)%3);

	  # translate it.
	  if( length($filler) != 6 ) {
	      my $lphase = $last_exon->end_phase;
	      my $fphase = $first_exon->phase;
	      $self->throw("Wrong length of filler seq. Error in coding [$filler] $lphase:$fphase\n");
	  }
	  my $fillerseq = Bio::Seq->new( -seq => $filler, -moltype => 'dna');
	  my $tfillerseq = $fillerseq->translate();
	  $seqstr .= $tfillerseq->seq;
      } 
      $seqstr .= $tseq->seq;
      $prevtrans = $ptrans;
  }
  
  $seqstr =~ s/\*$//g;

  my $trans_seq = Bio::Seq->new( -seq => $seqstr , '-id' => $self->translation->id() ) ;


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

    my $tmp = $exon->seq->seq();

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


  my $seq = new Bio::PrimarySeq('-seq' => $mrna,'-id' => $self->id);

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

   #$self->sort();
   my @exons = $self->each_Exon;
   my $exon_start = $exons[0];


   foreach my $exon ( @exons ) {


       # trim down start ends on the basis of phase.
       if( $prev && $prev->end_phase != $exon->phase ) {
	   $self->throw("Called coherent translate but exon phases don't match. Yuk!");
       }

       # warn about non DNA passed in. 

       if( $exon->entire_seq()->moltype ne 'dna' ) {
	   #$self->warn("Error. Whoever implemented this databases did not set type to Dna. Setting now!");
	   $exon->entire_seq()->moltype('dna');
       }
#       print STDERR "Exon phase " . $exon->phase . "\t" . $exon->start . "\t" . $exon->end . "\n";
#       print STDERR "Exon sequence is " . $exon->seq->seq . "\n";

       my $seq = $exon->seq();
       my $str = $seq->seq();
       

       if( length $str == 0 ) {
	   $self->throw("Bad internal error - got a 0 length rstring...");
       }

       $tstr .= $str;
   }

   $debug = 0;

   if ( $debug ) {
#       print STDERR "Bstr is $tstr\n";
#       print STDERR "Exon phase is " . $exon_start->phase . "\n";
       my @trans;
       my $exseq = new Bio::PrimarySeq(-seq => $tstr , '-id' => 'dummy' , -moltype => 'dna');
       	$trans[0] = $exseq->translate();

	# this is because a phase one intron leaves us 2 base pairs, whereas a phase 2
	# intron leaves one base pair.

	$trans[1] = $exseq->translate('*','X',2);
	$trans[2] = $exseq->translate('*','X',1);

#       print(STDERR "Exon start end " . $exon_start->start . " " . $exon_start->end . " " . $exon_start->phase . " " . $exon_start->strand ."\n");
#       print(STDERR "Translation 0 " . $trans[0]->seq . "\n");
#       print(STDERR "Translation 1 " . $trans[1]->seq . "\n");
#       print(STDERR "Translation 2 " . $trans[2]->seq . "\n");
   }


   if( $exon_start->phase == 1 ) {
       $tstr = substr $tstr, 2;
   } elsif ( $exon_start->phase == 2 ) {
       $tstr = substr $tstr, 1;
   } 

   if ( $debug ) {
       print STDERR "Exon start phase is " . $exon_start->phase . "\n";
       print STDERR "Tstr is $tstr\n";
   }

   # phase 0 - no need.


   my $temp_seq = Bio::Seq->new( -seq => $tstr , '-id' => 'temp', -moltype => 'dna' );
  #my $trans_seq = $temp_seq->translate();

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
    
	my $tex=$ex->translate;
	my $pep=$tex->seq;
    
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

=head2 start_exon

 Title   : start_exon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub start_exon{
   my ($self,@args) = @_;

   return ${$self->{'_trans_exon_array'}}[0];

}

=head2 end_exon

 Title   : end_exon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub end_exon{
   my ($self,@args) = @_;
   return ${$self->{'_trans_exon_array'}}[$#{$self->{'_trans_exon_array'}}];
}


=head2 created

 Title   : created
 Usage   : $obj->created($newval)
 Function: 
 Returns : value of created
 Args    : newvalue (optional)


=cut

sub created{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'created'} = $value;
    }
    return $obj->{'created'};

}


=head2 modified

 Title   : modified
 Usage   : $obj->modified($newval)
 Function: 
 Returns : value of modified
 Args    : newvalue (optional)


=cut

sub modified{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'modified'} = $value;
    }
    return $obj->{'modified'};

}



1;





