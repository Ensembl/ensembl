
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

use Bio::Root::RootI;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Translation;
use Bio::DBLinkContainerI;


@ISA = qw(Bio::Root::RootI Bio::DBLinkContainerI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub new {
  my($class,@args) = @_;

  if( ref $class ) { 
      $class = ref $class;
  }

  my $self = {};
  bless $self,$class;

  $self->{'_trans_exon_array'} = [];
  $self->{'_db_link'} = [];

  # set stuff in self from @args
  foreach my $a (@args) {
    $self->add_Exon($a);
  }

  #$self->is_partial(0);

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

sub first_exon {
    my ($self) = @_;
    {
        my($pkg, $file, $line) = caller(0);
        $self->warn("In file '$file', package '$pkg' line $line\n".
             "please switch your code to call 'start_exon'\n".
             "'first_exon' is now deprecated\n");
    }
    return $self->start_exon;
}

sub last_exon {
    my ($self) = @_;
    {
        my($pkg, $file, $line) = caller(0);
        $self->warn("In file '$file', package '$pkg' line $line\n".
             "please switch your code to call 'end_exon'\n".
             "'last_exon' is now deprecated\n");
    }
    return $self->end_exon;
}

=head2 translatable_exons

 Title   : translatable_exons
 Usage   : @exons = $transcript->translateable_exons
 Function: return list of exons that translate with the start- and
           end-exons truncated to the CDS regions.
 Example :
 Returns : 
 Args    :


=cut

sub translateable_exons {
    my( $self ) = @_;
    
    my $translation = $self->translation
        or $self->throw("No translation attached to transcript object");
    my $start_exon_id   = $translation->start_exon_id;
    my $end_exon_id     = $translation->end_exon_id;
    my $t_start         = $translation->start;
    my $t_end           = $translation->end;
    
    my( @translateable );
    foreach my $ex ($self->each_Exon) {
        my $ex_id   = $ex->id;
        
        if ($ex_id ne $start_exon_id and ! @translateable) {
            next;   # Not yet in translated region
        }
        
        my $start   = $ex->start;
        my $end     = $ex->end;
        my $length  = $ex->length;
        my $strand  = $ex->strand;
        
        my $trunc_start = $start;
        my $trunc_end   = $end;
        
        # Adjust the translation start if this is the start exon
        if ($ex_id eq $start_exon_id) {
            if ($t_start < 1 or $t_start > $length) {
                $self->throw("Translation start '$t_start' is outside exon '$ex_id' length=$length");
            }
            if ($strand == 1) {
                $trunc_start = $start + $t_start - 1;
            } else {
                $trunc_end   = $end   - $t_start + 1;
            }
        }
        
        # Adjust the translation end if this is the end exon
        if ($ex_id eq $end_exon_id) {
            if ($t_end < 1 or $t_end > $length) {
                $self->throw("Translation end '$t_end' is outside exon '$ex_id' length=$length");
            }
            if ($strand == 1) {
                $trunc_end   = $end   - $length + $t_end;
            } else {
                $trunc_start = $start + $length - $t_end;
            }
        }
        
        # Make a truncated exon if the translation start or
        # end causes the coordinates to be altered.
        if ($trunc_start != $start or $trunc_end != $end) {
            my $new_exon = Bio::EnsEMBL::Exon->new;
            $new_exon->id($ex_id);
            
            # Use the new start and end
            $new_exon->start($trunc_start);
            $new_exon->end  ($trunc_end);
            
            $new_exon->strand($strand);
            $new_exon->clone_id($ex->clone_id);
            $new_exon->contig_id($ex->contig_id);
            $new_exon->phase($ex->phase);
            $new_exon->attach_seq($ex->entire_seq);
            push(@translateable, $new_exon);
        } else {
            # It's just an ordinary internal exon
            push(@translateable, $ex);
        }
        
        # Exit the loop when we've found the last exon
        last if $ex_id eq $end_exon_id;
    }
    
    return @translateable;
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

sub split_Transcript_to_Partial {
   my ($self,$on_translate) = @_;


   if( $on_translate == 1 && ! defined $self->translation ) {
       $self->throw("Attempting to split Transcript on translation, but no translation...");
   }
   # should we check here if the Exons know their entire_seq?

   my @exons;
   if( $on_translate == 1 ) {
       @exons = $self->translateable_exons();
   } else {
       @exons = $self->each_Exon;
   }

   # one exon genes - easy to handle. (unless of course they have UTRs ...)
   if (@exons == 1) {
      # can't just return self - spliced UTR exons should not be returned or they'll be translated ...
     # return $self;
     
     my $t = $self->new();
     $t->id($self->id);
     
     $t->add_Exon($exons[0]);
     return $t;
   }

   # find the start exon;

   my $prev = shift @exons;

   my $t;
   my @out;


   TRANSCRIPT :
   while (@exons) {

       # make a new transcript, add the old exon
       $t = $self->new();
       $t->id($self->id);

       $t->add_Exon($prev);
       $t->is_partial(1);

       push(@out,$t);

       while( my $exon = shift @exons ) {
	   if( $exon->phase == $prev->end_phase) {
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

   if (@out == 1) {
       $t->is_partial(0);
   }

   return @out;
}                                       # split_Transcript_to_Partial

=head2 translate

 Title   : translate
 Usage   : $pep = $feat->translate()

 Function: returns the peptide translation of the gene - in the correct phase. 
           This method now requires that the Transcript has a valid Translation
           object.

 Returns : Bio::Seq
 Args    : none

=cut

sub translate {
  my ($self) = @_;

  my $debug;


  if ( ! defined $self->translation ) {
      $self->throw("You have to have a translation now to make a translation. Doh!");
  }

  my @trans = $self->split_Transcript_to_Partial(1);

  unless (@trans) {
      $self->throw("Bad internal error - split a transcript to zero transcripts! Doh!");
  }
    #if (@trans > 1) {
    #    warn "Made ", scalar(@trans), " partial transcripts";
    #}
  

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
	  my $last_exon  = $prevtrans->end_exon();
	  my $first_exon = $ptrans   ->start_exon();
	  my $filler;

	  # last exon
	  if( $last_exon->end_phase != 0 ) {
	      $filler = substr($last_exon->seq->seq, $last_exon->seq->length - $last_exon->end_phase);
	  } 
	  $filler .= 'N' x (3 - $last_exon->end_phase);

	  # do first exon now.

	  if( $first_exon->phase != 0 ) {
	      $filler .= 'N' x $first_exon->phase;
	  } else {
	      $filler .= 'NNN';
	  }
	  $filler .= substr($first_exon->seq->seq, 0, (3 - $first_exon->phase) % 3);

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
  
  $seqstr =~ s/\*$//;

  my $trans_seq = Bio::Seq->new( -seq => $seqstr , '-id' => $self->translation->id() ) ;


  return $trans_seq;
}

=head2 seq

Returns a Bio::Seq object which consists of just
the sequence of the exons concatenated together,
without messing about with padding with N's from
Exon phases like B<dna_seq> does.

=cut

sub seq {
    my( $self ) = @_;
    
    my $transcript_seq_string = '';
    foreach my $ex ($self->each_Exon) {
        $transcript_seq_string .= $ex->seq->seq;
    }
    
    my $seq = Bio::Seq->new(
        -DISPLAY_ID => $self->id,
        -MOLTYPE    => 'dna',
        -SEQ        => $transcript_seq_string,
        );
    return $seq;
}

#PL: prolly should read 'cDNA' for 'mRNA' below:

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
  
  # JGRG - $strand unused, so commented out
  #my $strand = $self->{_trans_exon_array}[0]->strand;

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

    $prev = $exon;  #heikki: this line was missing!  
    $mrna  .= $tmp; 
  }

  my $seq = new Bio::PrimarySeq(-SEQ => $mrna,-ID => $self->id);

  return $seq;
}

=head2 contig_dna

  Title   : contig_dna
  Usage   : $tran->contig_dna($dna);
   Function: Sets the dna sequence of the contig
  Returns : Bio::Seq
  Args    : Bio::Seq

=cut

# PL: why is this here? Is this OK when using VirtualContigs? 
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

 Usage : <internal function> translates a coherent transcript.  Uncoherent
           transcripts need to be broken up with split to partial first.
           (coherent means: having the same intron phase, i.e., all
           contiguous exon pairs that don not have phase differences).

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
#       print STDERR "Exon phase " . $exon->id ." " . $exon->phase . "\t" . $exon->start . "\t" . $exon->end . " " .$exon->strand. " ".$exon->entire_seq->id ."\n";
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
       my $exseq = new Bio::PrimarySeq(-SEQ => $tstr , '-id' => 'dummy' , -moltype => 'dna');
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


   my $temp_seq = Bio::Seq->new( -SEQ => $tstr , '-id' => 'temp', -moltype => 'dna' );
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
  return new Bio::Seq(-SEQ => $tmpseq);
}

  
sub pep_coords {
    my $self = shift;

    # for mapping the peptide coords back onto the dna sequence
    # it would be handy to have a list of the peptide start end coords
    # for each exon
  
    my @starts;
    my @ends;
  
    my $fullpep = $self->translate()->seq;

    
    foreach my $ex ($self->translateable_exons) {

	my $tex=$ex->translate;
	

	my $pep=$tex->seq;
	$pep =~ s/X$//g;
	
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
 Usage   : $start_exon = $transcript->start_exon;
 Returns : The first exon in the transcript.
 Args    : NONE


=cut

sub start_exon{
   my ($self,@args) = @_;

   return ${$self->{'_trans_exon_array'}}[0];

}

=head2 end_exon

 Title   : end_exon
 Usage   : $end_exon = $transcript->end_exon;
 Returns : The last exon in the transcript.
 Args    : NONE


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

# sneaky web only function...
sub gene_is_known {
    my ($self,$value) = @_;
    
    if( defined $value ) {
         $self->{'_web_hack_gene_is_known'} = $value;
    }
    
    return $self->{'_web_hack_gene_is_known'};
}

=head2 rna_pos

  Title   : rna_pos
  Usage   : $loc = $feat->dna_seq(23456)
  Function: Translates genomic coordinates into mRNA coordinates
  Returns : integer
  Args    : integer, genomic location

=cut

sub rna_pos {
    my ($self, $loc) = @_;

    my $start = $self->start_exon->start;
    #test that loc is within  mRNA
    return undef if $loc < $start;
    return undef if $loc >= $self->end_exon->end;

    my $mrna = 1;

    my $prev = undef;
    foreach my $exon ($self->each_Exon) {
	
	my $tmp = length $exon->seq->seq();
	#$tmp -= $exon->phase if not $prev;

	# we now have to figure out if the phase is compatible. If it
	# is not, we need to add some stuff in...

	if( $prev ) {
	    if( $prev->end_phase != $exon->phase ) {
		if( $prev->end_phase == 0 ) {
		    if( $exon->phase == 1 ) {
			$mrna += 2;
		    }

		    if( $exon->phase == 2 ) {
			$mrna += 1;
		    }
		} elsif ( $prev->end_phase == 1 ) {
		    if( $exon->phase == 0 ) {
			$mrna += 2;
		    }
		    
		    if( $exon->phase == 2 ) {
			$mrna += 1;
		    }
		} elsif ( $prev->end_phase == 2 ) {
		    if( $exon->phase == 0 ) {
			$mrna += 1;
		    }
		    
		    if( $exon->phase == 1 ) {
			$mrna += 2;
		    }
		} else {
		    $self->warn("Impossible phases in calculating fixing stuff");
		}
	    }
	} # end of if previous is there

	if ($loc < $exon->end) {
	    return $loc - $exon->start + $mrna ;
	}
	$mrna  += $tmp;
	$prev = $exon;
    }
    #return $mrna;
}

=head2 dn_length

  Title   : dna_length
  Usage   : $loc = $feat->dna_length;
  Function: return the length of the transcript''s DNA
  Returns : integer
  Args    : nn

=cut

sub dna_length {
     my ($self) = @_;

     # # not setting:
     # if( defined $value ) {
     #      $self->{'_dna_length'} = $value;
     # }

     if (! defined $self->{'_dna_length'}) { 
         # get from dna_seq;
         $self->{'_dna_length'} = $self->dna_seq->length;         
     }
     return $self->{'_dna_length'};
}

sub gene_name {
     my ($self,$value) = @_;
     
     if( defined $value ) {
          $self->{'_web_hack_gene_name'} = $value;
     }
     
     return $self->{'_web_hack_gene_name'};
}

=head2 description

 Title   : description
 Usage   : $obj->description($newval)
 Function: 
 Returns : value of description
 Args    : newvalue (optional)


=cut

sub description{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'description'} = $value;
    }
    return $obj->{'description'};

}



=head2 each_Exon_in_context

 Title   : each_Exon_in_context
 Usage   : @exons = $t->each_Exon_in_context($vc->id)
 Function: returns exons with this particular context (aka seqname)
 Example :
 Returns : 
 Args    :


=cut

sub each_Exon_in_context{
   my ($self,$context) = @_;

   my @exons;

   foreach my $exon ( $self->each_Exon ) {
       if( $exon->seqname eq $context ) {
	   push(@exons,$exon);
       }
   }

   return @exons;
}

=head2 is_start_exon_in_context

 Title   : is_start_exon_in_context
 Usage   : if( $t->is_start_exon_in_context($vc->id) ==0 ) {
               # transcript runs off this VC
 Function: returns 1 or 0 depending whether the start exon is 
           in this context or not
 Example :
 Returns : 
 Args    :


=cut

sub is_start_exon_in_context{
   my ($self,$context) = @_;

   if( $self->start_exon->seqname eq $context ) {
       return 1;
   } else {
       return 0;
   }

}


=head2 is_end_exon_in_context

 Title   : is_end_exon_in_context
 Usage   : if( $t->is_end_exon_in_context($vc->id) ==0 ) {
               # transcript runs off this VC
 Function: returns 1 or 0 depending whether the end exon is 
           in this context or not
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub is_end_exon_in_context{
   my ($self,$context) = @_;

   if( $self->end_exon->seqname eq $context ) {
       return 1;
   } else {
       return 0;
   }

}


=head2 strand_in_context

 Title   : strand_in_context
 Usage   : $strand = $t->strand_in_context($vc->id)
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub strand_in_context{
   my ($self,$context) = @_;

   my @exons = $self->each_Exon_in_context($context);

   if( scalar(@exons) == 0 ) {
       $self->warn("There are no exons in this context. Bad to call this - returning strand 0 from strand_in_context on Transcript");
       return 0;
   }

   return $exons[0]->strand;
}




1;
