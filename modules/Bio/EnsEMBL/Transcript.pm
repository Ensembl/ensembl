#
# BioPerl module for Transcript
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
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

     # Returns an array of Exon objects
     my @exons = @{$tran->get_all_Exons}     
     # Returns the peptide translation of the exons as a Bio::Seq
     my $pep   = $tran->translate()       
     # Sorts exons into order (forward for + strand, reverse for - strand)
     $tran->sort()                        

=head1 CONTACT

Email questions to the ensembl developer mailing list <ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Transcript;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Intron;
use Bio::EnsEMBL::Translation;
use Bio::DBLinkContainerI;
use Bio::EnsEMBL::TranscriptI;
use Bio::EnsEMBL::Mapper;

@ISA = qw(Bio::EnsEMBL::Root Bio::DBLinkContainerI Bio::EnsEMBL::TranscriptI);
# new() is inherited from Bio::Root::Object

sub new {
  my($class,@args) = @_;

  if( ref $class ) { 
      $class = ref $class;
  }

  my $self = {};
  bless $self,$class;

  $self->{'_trans_exon_array'} = [];

  # set stuff in self from @args
  foreach my $a (@args) {
    $self->add_Exon($a);
  }


  return $self; # success - we hope!
}



=head2 get_all_DBLinks

  Arg [1]    : 
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub get_all_DBLinks {
  my $self = shift;

  if( !defined $self->{'_db_link'} ) {
    $self->{'_db_link'} = [];
    if( defined $self->adaptor ) {
      $self->adaptor->db->get_DBEntryAdaptor->fetch_all_by_Transcript($self);
    }
  } 
  
  return $self->{'_db_link'};
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

   unless(defined $value && ref $value && 
	  $value->isa('Bio::Annotation::DBLink') ) {
     $self->throw("This [$value] is not a DBLink");
   }

   if( !defined $self->{'_db_link'} ) {
     $self->{'_db_link'} = [];
   }

   push(@{$self->{'_db_link'}},$value);
}



sub dbID {
   my $self = shift;
   
   if( @_ ) {
      my $value = shift;
      $self->{'dbID'} = $value;
    }
    return $self->{'dbID'};

}

=head2 external_db

 Title   : external_db
 Usage   : $ext_db = $obj->external_db();
 Function: external_name if available
 Returns : the external db link for this transcript
 Args    : new external db (optional)

=cut

sub external_db {
  my ($self, $arg ) = @_;

  if( defined $arg ) {
    $self->{'_external_db'} = $arg;
  }

  return $self->{'_external_db'};
}


=head2 external_name

 Title   : external_name
 Usage   : $ext_name = $obj->external_name();
 Function: external_name if available
 Example : 
 Returns : the external name of this transcript
 Args    : new external name (optional)

=cut

sub external_name {
  my ($self, $arg ) = @_;

  if( defined $arg ) {
    $self->{'_external_name'} = $arg;
  }

  return $self->{'_external_name'};
}


sub is_known {
  my $self = shift;
  if( defined $self->external_name() && $self->external_name() ne '' ) {
    return 1;
  } else {
    return 0;
  }
}


sub type {
  my ($self, $type) = @_;

  if(defined $type) {
    $self->{'_type'} = $type;
  }

  return $self->{'_type'};
}

sub adaptor {
   my $self = shift;
   
   if( @_ ) {
      my $value = shift;
      $self->{'adaptor'} = $value;
    }
    return $self->{'adaptor'};

}

=head2 _translation_id

 Title   : _translation_id
 Usage   : $obj->_translation_id($newval)
 Function: 
 Returns : translation objects dbID
 Args    : newvalue (optional)


=cut

sub _translation_id {
   my $self = shift;
   
   if( @_ ) {
      my $value = shift;
      $self->{'_translation_id'} = $value;
    }
    return $self->{'_translation_id'};

}


=head2 translation

 Title   : translation
 Usage   : $obj->translation($newval)
 Function: 
 Returns : value of translation
 Args    : newvalue (optional)


=cut

sub translation {
  my $self = shift;
  if( @_ ) {
    my $value = shift;
    if( ! ref $value || !$value->isa('Bio::EnsEMBL::Translation') ) {
      $self->throw("This [$value] is not a translation");
    }
    $self->{'translation'} = $value;
  } else {
    if( ! defined $self->{'translation'} &&
	defined $self->_translation_id() ) {
      $self->{'translation'} = 
	$self->adaptor->db->get_TranslationAdaptor->fetch_by_dbID( 
					    $self->_translation_id(), $self );
    }
  }
  return $self->{'translation'};
}

=head2 start

 Description: it returns the start coordinate of the lef-most exon, i.e.
              the 5prime exon in the forward strand and the 3prime exon in the reverse strand

=cut


sub start {
  my $self = shift;
  my $arg = shift;
  
  my $strand;
  my $start;
  if( defined $arg ) {
    $self->{'_start'} = $arg;
  } elsif(!  defined $self->{'_start'} ) {

    $strand = $self->start_Exon->strand();
    if( $strand == 1 ) {
      $start = $self->start_Exon->start();
    } else {
      $start = $self->end_Exon->start();
    }
    $self->{'_start'} = $start;
  }
  
  return $self->{'_start'};
}


sub end {
  my $self = shift;
  my $arg = shift;

  my $strand;
  my $end
;
  if( defined $arg ) {
    $self->{'_end'} = $arg;
  } elsif( ! defined $self->{'_end'} ) {
    $strand = $self->start_Exon->strand();
    if( $strand == 1 ) {
      $end = $self->end_Exon->end();
    } else {
      $end = $self->start_Exon->end();
    }
    $self->{'_end'} = $end;
  }
  
  return $self->{'_end'};
}


=head2 spliced_seq

  Args       : none
  Example    : none
  Description: retrieves all Exon sequences and concats them together. No phase padding magic is 
               done, even if phases dont align.
  Returntype : txt
  Exceptions : none
  Caller     : general

=cut

sub spliced_seq {
  my ( $self ) = @_;
  
  my $seq_string = "";
  for my $ex ( @{$self->get_all_Exons()} ) {
    $seq_string .= $ex->seq()->seq();
  }

  return $seq_string;
}


=head2 translateable_seq

  Args       : none
  Example    : none
  Description: returns a string with the translateable part of the
               Sequence. It magically pads the exon sequences with
               N if the phases of the Exons dont align
  Returntype : txt
  Exceptions : none
  Caller     : general

=cut

sub translateable_seq {
  my ( $self ) = @_;

  my $mrna = "";
  my $lastphase = 0;
  my $first = 1;

  foreach my $exon (@{$self->get_all_translateable_Exons()}) {

    my $phase = 0;
    if (defined($exon->phase)) {
      $phase = $exon->phase;
    }
    
    # startpadding is needed if MONKEY_EXONS are on
    if( $first && (! defined $ENV{'MONKEY_EXONS'}) ) {
      $mrna .= 'N' x $phase;
      $first = 0;
    }

    if( $phase != $lastphase && ( defined $ENV{'MONKEY_EXONS'})) {
      # endpadding for the last exon
      if( $lastphase == 1 ) {
	$mrna .= 'NN';
      } elsif( $lastphase == 2 ) {
	$mrna .= 'N';
      }
      #startpadding for this exon
      $mrna .= 'N' x $phase;
    }
    $mrna .= $exon->seq->seq();
    $lastphase = $exon->end_phase();
  }
  return $mrna;
}



sub coding_start {
  my $self = shift;
  my $arg = shift;

  my $strand;
  my $start;

  
  if( defined $arg ) {
    $self->{'coding_start'} = $arg;
  } elsif( ! defined $self->{'coding_start'} && 
	   defined $self->translation() ) {
    $strand = $self->translation()->start_Exon->strand();
    if( $strand == 1 ) {
      $start = $self->translation()->start_Exon->start();
      $start += ( $self->translation()->start() - 1 );
    } else {
      $start = $self->translation()->end_Exon->end();
      $start -= ( $self->translation()->end() - 1 );
    }
    $self->{'coding_start'} = $start;
  }

  return $self->{'coding_start'};
}


sub coding_end {
  my $self = shift;
  my $arg = shift;

  my $strand;
  my $end;
  
  if( defined $arg ) {
    $self->{'coding_end'} = $arg;
  } elsif( ! defined $self->{'coding_end'} && 
	   defined $self->translation() ) {
    $strand = $self->translation()->start_Exon->strand();
    if( $strand == 1 ) {
      $end = $self->translation()->end_Exon->start();
      $end += ( $self->translation()->end() - 1 );
    } else {
      $end = $self->translation()->start_Exon->end();
      $end -= ( $self->translation()->start() - 1 );
    }
    $self->{'coding_end'} = $end;
  }
  
  return $self->{'coding_end'};
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
   unless(defined $exon && ref $exon && $exon->isa("Bio::EnsEMBL::Exon") ) {
     $self->throw("[$exon] is not a Bio::EnsEMBL::Exon!");
   }

   #invalidate the start, end and strand - they may need to be recalculated
   $self->{'_start'} = undef;
   $self->{'_end'} = undef;
   $self->{'_strand'} = undef;

   push(@{$self->{'_trans_exon_array'}},$exon);
}



=head2 get_all_Exons

  Arg [1]    : none
  Example    : my @exons = @{$transcript->get_all_Exons()};
  Description: Returns an listref of the exons in this transcipr in order.
               i.e. the first exon in the listref is the 5prime most exon in 
               the transcript.
  Returntype : a list reference to Bio::EnsEMBL::Exon objects
  Exceptions : none
  Caller     : general

=cut

sub get_all_Exons {
   my ($self) = @_;

   return $self->{'_trans_exon_array'};
}



=head2 length


    my $t_length = $transcript->length

Returns the sum of the length of all the exons in
the transcript.

=cut

sub length {
    my( $self ) = @_;
    
    my $length = 0;
    foreach my $ex (@{$self->get_all_Exons}) {
        $length += $ex->length;
    }
    return $length;
}



=head2 get_all_Introns

  Args       : none
  Example    : @introns = @{$transcript->get_all_Introns()};
  Description: Returns an listref of Bio::EnsEMBL::Intron objects.  The result 
               is not cached in any way, so calling each_Intron multiple times
               will create new Intron objects (although they will, of course, 
               have the same properties).
  Returntype : list reference to Bio::EnsEMBL::Intron objects
  Exceptions : none
  Caller     : general

=cut

sub get_all_Introns {
    my( $self ) = @_;
    
    my @exons = @{$self->get_all_Exons};
    my $last = @exons - 1;
    my( @int );
    for (my $i = 0; $i < $last; $i++) {
        my $intron = Bio::EnsEMBL::Intron->new;
        $intron->upstream_Exon  ($exons[$i]    );
        $intron->downstream_Exon($exons[$i + 1]);
        push(@int, $intron);
    }
    return \@int;
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

   $self->{'_start'} = undef;
   $self->{'_end'} = undef;
   $self->{'_strand'} = undef;

   $self->{'_trans_exon_array'} = [];
}



=head2 five_prime_utr and three_prime_utr

    my $five_prime  = $transcrpt->five_prime_utr
        or warn "No five prime UTR";
    my $three_prime = $transcrpt->three_prime_utr
        or warn "No three prime UTR";

These methods return a B<Bio::Seq> object
containing the sequence of the five prime or
three prime UTR, or undef if there isn't a UTR.

Both method throw an exception if there isn't a
translation attached to the transcript object.

=cut

sub five_prime_utr {
    my( $self ) = @_;
    
    my $translation = $self->translation
        or $self->throw("No translation attached to transcript object");
    my $start_exon_id   = $translation->start_Exon->stable_id;
    my $t_start         = $translation->start;
    
    my $seq_string = '';
    foreach my $ex (@{$self->get_all_Exons}) {
        if (((defined $ex->stable_id)&&($ex->stable_id eq $start_exon_id))  
	                        # The criteria for the world with stable_ids.
	    ||($ex == $translation->start_Exon)) {
	                        # The criteria for the genebuild where stable_ids are not always about.
            my $start   = $ex->start;
            my $end     = $ex->end;
            my $strand  = $ex->strand;
            
            if ($strand == 1) {
                $end   = $start + $t_start - 2;
            } else {
                $start = $end   - $t_start + 2;
            }
            
            if ($start <= $end) {
                my $utr_exon = Bio::EnsEMBL::Exon->new;
                $utr_exon->start($start);
                $utr_exon->end($end);
                $utr_exon->strand($strand);
                $utr_exon->attach_seq($ex->entire_seq);
                $seq_string .= $utr_exon->seq->seq;
            }
            last;   # At end of UTR
        } else {
            $seq_string .= $ex->seq->seq;
        }
    }
    
    if ($seq_string) {
        my $seq = Bio::Seq->new;
        $seq->id($self->stable_id . '-five_prime_UTR');
        $seq->seq($seq_string);
        return $seq;
    } else {
        return;
    }
}

sub three_prime_utr {
    my( $self ) = @_;
    
    my $translation = $self->translation
        or $self->throw("No translation attached to transcript object");
    my $end_exon_id   = $translation->end_Exon->stable_id;
    my $t_end         = $translation->end;
    
    my $seq_string = '';
    my $in_utr = 0;
    foreach my $ex (@{$self->get_all_Exons}) {
        if ($in_utr) {
            $seq_string .= $ex->seq->seq;
        }
        elsif ((defined $ex->stable_id)&&($ex->stable_id eq $end_exon_id)
	                        # The criteria for the world with stable_ids.
	       ||($ex == $translation->end_Exon)) {
	                        # The criteria for the genebuild where stable_ids are not always about.
            $in_utr = 1;
            my $start   = $ex->start;
            my $end     = $ex->end;
            my $strand  = $ex->strand;
            
            if ($strand == 1) {
                $start = $start + $t_end;
            } else {
                $end   = $end   - $t_end;
            }
            
            if ($start <= $end) {
                my $utr_exon = Bio::EnsEMBL::Exon->new;
                $utr_exon->start($start);
                $utr_exon->end($end);
                $utr_exon->strand($strand);
                $utr_exon->attach_seq($ex->entire_seq);
                $seq_string .= $utr_exon->seq->seq;
            }
        }
    }
    
    if ($seq_string) {
        my $seq = Bio::Seq->new;
        $seq->id($self->stable_id . '-three_prime_UTR');
        $seq->seq($seq_string);
        return $seq;
    } else {
        return;
    }
}


=head2 get_all_translateable_Exons

  Args       : none
  Example    : none
  Description: Returns a list of exons that translate with the
               start and end exons truncated to the CDS regions.
               Will not work correctly if Exons are Sticky. 
  Returntype : listref Bio::EnsEMBL::Exon
  Exceptions : If there is no Translation object
  Caller     : Genebuild, $self->translate()

=cut


sub get_all_translateable_Exons {
  my ( $self ) = @_;

  my $translation = $self->translation
    or $self->throw("No translation attached to transcript object");
  my $start_exon      = $translation->start_Exon;
  my $end_exon        = $translation->end_Exon;
  my $t_start         = $translation->start;
  my $t_end           = $translation->end;

  my( @translateable );

  foreach my $ex (@{$self->get_all_Exons}) {

    if ($ex ne $start_exon and ! @translateable) {
      next;   # Not yet in translated region
    }

    my $length  = $ex->length;
        
    my $adjust_start = 0;
    my $adjust_end = 0;
    # Adjust to translation start if this is the start exon
    if ($ex == $start_exon ) {
      if ($t_start < 1 or $t_start > $length) {
	$self->throw("Translation start '$t_start' is outside exon $ex length=$length");
      }
      $adjust_start = $t_start - 1;
    }
        
    # Adjust to translation end if this is the end exon
    if ($ex == $end_exon) {
      if ($t_end < 1 or $t_end > $length) {
	$self->throw("Translation end '$t_end' is outside exon $ex length=$length");
      }
      $adjust_end = $t_end - $length;
    }

    # Make a truncated exon if the translation start or
    # end causes the coordinates to be altered.
    if ($adjust_end || $adjust_start) {
      my $newex = $ex->adjust_start_end( $adjust_start, $adjust_end );

      push( @translateable, $newex );
    } else {
      push(@translateable, $ex);
    }
        
    # Exit the loop when we've found the last exon
    last if $ex eq $end_exon;
  }
  return \@translateable;
}


sub translateable_exons {
    my( $self ) = @_;
  
    $self->warn( "Please use get_all_translateable_Exons(). Careful as it returns listref." );
    
    return @{$self->get_all_translateable_Exons()};
}




=head2 translate

  Args       : none
  Example    : none
  Description: return the peptide (plus eventuel stop codon) for this transcript.
               Does N padding of non phase matching exons. It uses translateable_seq
               internally. 
  Returntype : Bio::Seq
  Exceptions : If no Translation is set in this Transcript
  Caller     : general

=cut

sub translate {
  my ($self) = @_;

  my $mrna = $self->translateable_seq();
  my $display_id;

  if( defined $self->translation->stable_id ) {
    $display_id = $self->translation->stable_id;
  } elsif ( defined $self->temporary_id ) {
    $display_id = $self->temporary_id;
  } else {
    $display_id = $self->translation->dbID;
  }
	
  $mrna =~ s/TAG$|TGA$|TAA$//i;
  # the above line will remove the final stop codon from the mrna
  # sequence produced if it is present, this is so any peptide produced
  # won't have a terminal stop codon
  # if you want to have a terminal stop codon either comment this line out
  # or call translatable seq directly and produce a translation from it
  
  my $peptide = Bio::Seq->new( -seq => $mrna,
			       -moltype => "dna",
			       -id => $display_id );
  
 
  
  
  return $peptide->translate;
}

=head2 seq

Returns a Bio::Seq object which consists of just
the sequence of the exons concatenated together,
without messing about with padding with N\'s from
Exon phases like B<dna_seq> does.

=cut

sub seq {
    my( $self ) = @_;
    
    my $transcript_seq_string = '';
    foreach my $ex (@{$self->get_all_Exons}) {
#        $transcript_seq_string .= $ex->seq;
        $transcript_seq_string .= $ex->seq->seq;
    }
    
    my $seq = Bio::Seq->new(
        -DISPLAY_ID => $self->stable_id,
        -MOLTYPE    => 'dna',
        -SEQ        => $transcript_seq_string,
        );

    return $seq;
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
  my @exons = @{$self->get_all_Exons()};

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




  
sub pep_coords {
    my $self = shift;

    # for mapping the peptide coords back onto the dna sequence
    # it would be handy to have a list of the peptide start end coords
    # for each exon
  
    my ($p,$f,$l) = caller;
    $self->warn("$f:$l  Calls to pep_coords should no longer be necessary. Please use pep2genomic");
    my @starts;
    my @ends;
  
    my $fullpep = $self->translate()->seq;

    
    foreach my $ex ($self->translateable_exons) {

	my $tex=$ex->translate;
	

	my $pep=$tex->seq;
	$pep =~ s/X$//g;
	
	my $start = index($fullpep,$pep) + 1;
	
	my $end = $start + CORE::length($pep) - 1;
    
	push(@starts,$start);
	push(@ends,$end);
	
    }

    return \@starts,\@ends;
}



=head1 pep2genomic

  Arg  1   : integer start - relative to peptide
  Arg  2   : integer end   - relative to peptide

  Function : Provides a list of Bio::EnsEMBL::SeqFeatures which
             is the genomic coordinates of this start/end on the peptide

  Returns  : list of Bio::EnsEMBL::SeqFeature

=cut

sub convert_peptide_coordinate_to_contig {
  my ($self,$start,$end) = @_;
  my ($p,$f,$l) = caller;
  $self->warn("$f:$l convert_peptide_coo... is deprecated. Please use pep2genomic");

  $self->pep2genomic( $start, $end );
}

sub pep2genomic {
  my ($self,$start,$end) = @_;

  if( !defined $end ) {
    $self->throw("Must call with start/end");
  }

  # move start end into translate cDNA coordinates now.
  # much easier!
  $start = 3* $start-2;
  $end   = 3* $end;

  return $self->cdna2genomic( $start, $end );
}


sub cdna2genomic {
  my ($self,$start,$end) = @_;

  my $mapper;
  my @out;

  if( !defined $end ) {
    $self->throw("Must call with start/end");
  }

  $mapper = $self->_get_cdna_coord_mapper();
  
  my @mapped_coords = $mapper->map_coordinates( $self, $start, $end, 1, "cdna" );

  for my $coord ( @mapped_coords ) {
    if( $coord->isa( "Bio::EnsEMBL::Mapper::Coordinate" ) ) {
      my $sf = Bio::EnsEMBL::SeqFeature->new
	( -start => $coord->start(),
	  -end => $coord->end(),
	  -strand => $coord->strand()
	);
      $sf->contig( $coord->id() );
      push( @out, $sf );
    }
  } 
	  
  return @out;
}
  
sub _get_cdna_coord_mapper {
  my ( $self ) = @_;

  if( defined $self->{'_exon_coord_mapper'} ) {
    return $self->{'_exon_coord_mapper'};
  } 
  
  #
  # the mapper is loaded with OBJECTS in place of the IDs !!!!
  #  the objects are the contigs in the exons
  #

  my $mapper;
  $mapper = Bio::EnsEMBL::Mapper->new( "cdna", "genomic" );
  my @exons = @{$self->get_all_translateable_Exons() };
  my $start = 1;
  for my $exon ( @exons ) {
    $exon->load_genomic_mapper( $mapper, $self, $start );
    $start += $exon->length;
  }
  $self->{'_exon_coord_mapper'} = $mapper;
  return $mapper;
}

  

sub find_coord {
  my ($self,$coord,$type) = @_;
 
  my ($p,$f,$l) = caller;
  $self->warn("$f:$l find_coord is deprecated. Use pep2genomic");

  my $count = 0;
  my @exons = @{$self->get_all_Exons};
  my $end   = $#exons;
  my $dna;

  my ($starts,$ends) = $self->pep_coords;
  my $strand = $exons[0]->strand;

  # $starts and $ends are array refs containing the _peptide_ coordinates
  # of each exon. We may have 1 missing residue that spans an intron.
  # We ignore these.

  if ($strand == 1) {
    foreach my $ex (@{$self->get_all_Exons}) {
      
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

    foreach my $ex (@{$self->get_all_Exons}) {
      
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



=head2 start_Exon

 Title   : start_Exon
 Usage   : $start_exon = $transcript->start_Exon;
 Returns : The first exon in the transcript.
 Args    : NONE

=cut

sub start_Exon{
   my ($self,@args) = @_;

   return ${$self->{'_trans_exon_array'}}[0];
}

=head2 end_Exon

 Title   : end_exon
 Usage   : $end_exon = $transcript->end_Exon;
 Returns : The last exon in the transcript.
 Args    : NONE

=cut

sub end_Exon{
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

=head2 rna_pos

  Title   : rna_pos
  Usage   : $loc = $feat->dna_seq(23456)
  Function: Translates genomic coordinates into mRNA coordinates
            ARNE: padding probably not correct
  Returns : integer
  Args    : integer, genomic location

=cut

sub rna_pos {
    my ($self, $loc) = @_;

    my $start = $self->start_exon->start;
    #test that loc is within  mRNA
    return undef if $loc < $start;
    return undef if $loc >= $self->end_Exon->end;

    my $mrna = 1;

    my $prev = undef;
    foreach my $exon (@{$self->get_all_Exons}) {
	
	my $tmp = CORE::length( $exon->seq->seq());
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

=head2 dna_length

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


=head2 version

 Title   : version
 Usage   : $obj->version()
 Function: 
 Returns : value of version
 Args    : 

=cut

sub version{

    my ($self,$value) = @_;
    

    if( defined $value ) {
      $self->{'_version'} = $value;
    }

    if( exists $self->{'_version'} ) {
      return $self->{'_version'};
    }

    $self->_get_stable_entry_info();

    return $self->{'_version'};

}


=head2 stable_id

 Title   : stable_id
 Usage   : $obj->stable_id
 Function: 
 Returns : value of stable_id
 Args    : 


=cut

sub stable_id{

    my ($self,$value) = @_;
    

    if( defined $value ) {
      $self->{'_stable_id'} = $value;
      return;
    }

    if( exists $self->{'_stable_id'} ) {
      return $self->{'_stable_id'};
    }

    $self->_get_stable_entry_info();

    return $self->{'_stable_id'};

}

sub _get_stable_entry_info {
   my $self = shift;

   if( !defined $self->adaptor ) {
     return undef;
   }

   $self->adaptor->get_stable_entry_info($self);

}

=head2 temporary_id

 Title   : temporary_id
 Usage   : $obj->temporary_id($newval)
 Function: Temporary ids are used for Genscan predictions - which should probably
           be moved over to being stored inside the gene tables anyway. Bio::EnsEMBL::TranscriptFactory use this
 Example : 
 Returns : value of temporary_id
 Args    : newvalue (optional)


=cut

sub temporary_id{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'temporary_id'} = $value;
    }
    return $obj->{'temporary_id'};

}


sub finex_string {
   my ($self) = @_;

   my $finex;

   if ($self->stable_id ne "") {
     $finex = $self->stable_id;
   } else {
     $finex = $self->dbID;
   }

   $finex .= " ";

   my @exons = @{$self->get_all_Exons};

   $finex .= scalar(@exons) . " ";

   if ($exons[0]->strand == 1) {
      @exons = sort {$a->start <=> $b->start} @exons;
   } else {
      @exons = sort {$b->start <=> $a->start} @exons;
   }


   my $found_start = 0;
   my $found_end   = 0;

   foreach my $exon (@exons) {
     my $length = $exon->length;

     
     if ($exon == $self->translation->start_Exon &&
	 $exon == $self->translation->end_Exon) {
       $length = $self->translation->end - $self->translation->start + 1;

       $found_start = 1;
       $found_end   = 1;

       $finex .= $exon->phase . ":" . $exon->end_phase . ":" . $length . " ";

     } elsif ($exon == $self->translation->start_Exon) {
       $length = $exon->length - $self->translation->start + 1;
       $found_start = 1;

       $finex .= $exon->phase . ":" . $exon->end_phase . ":" . $length . " ";

     } elsif ($exon == $self->translation->end_Exon) {
       $length = $self->translation->end;
       $found_end = 1;

       $finex .= $exon->phase . ":" . $exon->end_phase . ":" . $length . " ";
       
     } elsif ($found_start == 1 && $found_end == 0) {
       $length = $exon->length;

       $finex .= $exon->phase . ":" . $exon->end_phase . ":" . $length . " ";
     }
   }

   $finex =~ s/\ $//;

   return $finex;
}


=head2 transform

  Arg  1    : hashref $old_new_exon_map
              a hash that maps old to new exons for a whole gene
  Function  : maps transcript in place to different coordinate system,
              It does so by replacing its old exons with new ones.
  Returntype: none
  Exceptions: none
  Caller    : Gene->transform()

=cut


sub transform {
  my $self = shift;
  my $href_exons = shift;
  my @mapped_list_of_exons;

  foreach my $exon (@{$self->get_all_Exons()}) {
    # the old exon was successfully remapped then store the new exon
    if ( exists $$href_exons{$exon} ) {
      push @mapped_list_of_exons, $$href_exons{$exon};
    }
    # but for the case where the exon was unable to be mapped, as it
    # was outside the bounds of the slice, include the original exon.
    else {
      push @mapped_list_of_exons, $exon;
    }
  }

  # flush the old list of exons
  $self->{'_trans_exon_array'} = [];

  # attach the new list of exons to the transcript
  push @{$self->{'_trans_exon_array'}},@mapped_list_of_exons;

  if( defined $self->{'translation'} ) {
    $self->translation->transform( $href_exons );
  }

  #invalidate the current start, end, strand - they need to be recalculated
  $self->{'_start'} = undef;
  $self->{'_end'} = undef;
  $self->{'_strand'} = undef;
  $self->{'_exon_coord_mapper'} = undef;
}



##########################################################
#
# sub DEPRECATED METHODS FOLLOW
#
##########################################################

sub translateable_dna {
   my ($self) = @_;

   my ($p,$f,$l) = caller;
   $self->warn("$f:$l translateable_dna isdeprecated. Please use translateable_seq (returns string)");

   my $tstr = $self->translateable_seq();
   my $temp_seq = Bio::Seq->new( -SEQ => $tstr , '-id' => 'temp', -moltype => 'dna' );

   return $temp_seq;
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

  $self->warn( "dna_seq is deprecated. Use spliced_seq to get concat off exon sequences as string." );

  my $mrna = $self->spliced_seq();
  my $seq = new Bio::PrimarySeq(-SEQ => $mrna,-ID => $self->stable_id);

  return $seq;
}

sub _translate_coherent{
   my ($self,$debug) = @_;
   
  $self->throw( "_translate_coherent is gone" );
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
  my ( $self ) = @_;
  my ($p,$f,$l) = caller;

  $self->warn( "$f:$l split_transcript_to_Partial deprecated, no replacement" );
  
}                                       # split_Transcript_to_Partial


=head2 is_partial

 Title   : is_partial
 Usage   : $obj->is_partial($newval)
 Function: 
 Example : 
 Returns : value of is_partial
 Args    : newvalue (optional)


=cut

sub is_partial{
   my ($self,$value) = @_;
   my ($p,$f,$l) = caller;

   $self->warn( "$f:$l is_partial definatly completely deprecated, just dont use it" );

   if( defined $value) {
      $self->{'is_partial'} = $value;
    }
    return $self->{'is_partial'};

}

=head2 id

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use stable_id or dbID instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub id {
   my $self = shift;
   my $value = shift;

   my ($p,$f,$l) = caller;
   $self->warn("$f:$l id deprecated. Please choose from stable_id or dbID, or, for predictions, temporary id");

  # catch adaptorless genes
  if( defined $value ) {
    $self->warn("$f:$l stable ids are loaded separately and dbIDs are generated on writing. Ignoring set value $value");
    return;
  }

   if( defined $self->stable_id ) {
     return $self->stable_id();
   } elsif ( defined $self->temporary_id ) {
     return $self->temporary_id; 
   } else {
     return $self->dbID;
   }

}


=head2 each_DBLink

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use get_all_DBLinks instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub each_DBLink {
   my ($self,@args) = @_;

   $self->warn("each_DBLink has been renamed get_all_DBLinks");
   return $self->get_all_DBLinks(@args);
}


=head2 contig_dna

  Title   : DEPRECATED contig_dna
  Usage   : DEPRECATED $tran->contig_dna($dna);
  Function: DEPRECATED Sets the dna sequence of the contig
  Returns : DEPRECATED Bio::Seq
  Args    : DEPRECATED Bio::Seq

=cut

# PL: why is this here? Is this OK when using VirtualContigs? 
sub contig_dna {
  my ($self,$dna) = @_;

  $self->warn("Transcript->contig_dna() is deprecated. \n");

  if (defined($dna)) {
    $self->{_contig_dna} = $dna;
  }

  return $self->{_contig_dna};
}


=head2 gene_is_known

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut
# sneaky web only function...
sub gene_is_known {
    my ($self,$value) = @_;
  
    $self->warn("Transcript->gene_is_known is deprecated. \n");
  
    if( defined $value ) {
         $self->{'_web_hack_gene_is_known'} = $value;
    }
    
    return $self->{'_web_hack_gene_is_known'};
}



=head2 each_Exon

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use get_all_Exons instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub each_Exon{
   my ($self) = @_;

   my ($p,$f,$l) = caller;
   $self->warn("$f:$l each_Exon is deprecated. Please use get_all_Exons");

   return @{$self->{'_trans_exon_array'}};
}



=head2 each_Exon_in_context

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use get_all_Exons and search through for seqname 
               yourself
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub each_Exon_in_context{
   my ($self,$context) = @_;

  $self->warn("Transcript->each_Exon_in_context is deprecated, use " .
              "get_all_Exons and search through them for matching seqnames " .
              "yourself");

   my @exons;

   foreach my $exon ( @{$self->get_all_Exons} ) {
       if( $exon->seqname eq $context ) {
	   push(@exons,$exon);
       }
   }

   return @exons;
}



=head2 is_start_exon_in_context

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub is_start_exon_in_context{
   my ($self,$context) = @_;

  $self->warn("Transcript->is_start_exon_in_context is deprecated");

   if( $self->start_exon->seqname eq $context ) {
       return 1;
   } else {
       return 0;
   }

}



=head2 is_end_exon_in_context

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub is_end_exon_in_context{
   my ($self,$context) = @_;

  $self->warn("Transcript->is_end_Exon_in_context is deprecated");

   if( $self->end_exon->seqname eq $context ) {
       return 1;
   } else {
       return 0;
   }

}


=head2 strand_in_context

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub strand_in_context{
   my ($self,$context) = @_;

  $self->warn("Transcript->strand_in_context is deprecated");

   my @exons = $self->each_Exon_in_context($context);

   if( scalar(@exons) == 0 ) {
	print STDERR "TRANSCRIPT: ", $self->stable_id,"\n";

       $self->warn("There are no exons in this context. Bad to call this - returning strand 0 from strand_in_context on Transcript");
       return 0;
   }

   return $exons[0]->strand;
}



=head2 exon_dna

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use $exon->seq->seq instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub exon_dna {
  my ($self,$exon) = @_;

  $self->warn("You are calling exon_dna - you really don't need this method.".
	      "\$exon->seq->seq should work (if not, attach the contig to the".
	      "exon objects correctly");

  my $tmpseq = $self->contig_dna->str($exon->start,$exon->end);
  
  if ($exon->strand == -1) {
    $tmpseq =~ tr/ATGCatgc/TACGtacg/;
    $tmpseq = reverse($tmpseq);
  }
  return new Bio::Seq(-SEQ => $tmpseq);
}


=head2 first_exon

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use start_exon instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

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


=head2 last_exon

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use end_exon instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

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

=head2 get_Exon_by_dbID

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use get_all_Exons instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_Exon_by_dbID {
   my ($self,$exonid) = @_;

   $self->warn("get_Exon_by_dbID is deprecated. Search through the exons " .
	       "returned from get_all_Exons yourself\n");

   if( !defined $exonid ) {
      $self->throw("Must call with exonid");
   }

   # not nice - linear search
   foreach my $exon ( @{$self->get_all_Exons()} ) {

      if( $exon->dbID eq $exonid ) {
          return $exon;
      }
   }

   return undef;
}


=head2 each_Intron

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use get_all_Introns instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub each_Intron {
  my $self = shift;
  $self->warn('each_Intron has been renamed get_all_Introns');
  return $self->get_all_Introns();
}



=head2 number

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use scalar @{transcript->get_all_Exons}; instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub number{
   my ($self,@args) = @_;

   $self->warn("number is deprecated use: " .
	       "scalar @{transcript->get_all_Exons}; \n");

   return scalar(@{$self->{'_trans_exon_array'}});   
}



=head2 gene_name

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub gene_name {
     my ($self,$value) = @_;
     
     $self->warn("Bio::EnsEMBL::Transcript::gene_name is deprecated do not use"
		 . "\n" . caller);

     if( defined $value ) {
          $self->{'_web_hack_gene_name'} = $value;
     }
     
     return $self->{'_web_hack_gene_name'};
}


=head2 start_exon

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use start_Exon instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub start_exon {
  my ($self, @args) = @_;

  $self->warn("start_exon has been renamed start_Exon\n");
  return $self->start_Exon(@args);
}


=head2 end_exon

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use end_Exon instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub end_exon {
  my ($self, @args) = @_;

  $self->warn("end_exon has been renamed end_Exon\n");
  return $self->end_Exon(@args);
}

1;

__END__;
