
#
# BioPerl module for DB::VirtualPrimarySeq
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DB::VirtualPrimarySeq - Object that "pretends" to be a PrimarySeq without actual phyiscally having the DNA.

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...
package Bio::EnsEMBL::DB::VirtualPrimarySeq;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::PrimarySeqI;

@ISA = qw(Bio::Root::Object Bio::PrimarySeqI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my ($vmap,$un,$clone,$length) = $self->_rearrange([qw(VMAP UN CLONE LENGTH)],@args);

  if (!$vmap->isa('Bio::EnsEMBL::DB::VirtualMap')) {
      $self->throw("$vmap is not a Bio::EnsEMBL::DB::VirtualMap!");
  }
  if (! defined $un) {
      $self->throw("Need to provide a unique number for the VirtualPrimarySeq id");
  }
  if (! defined $length) {
      $self->throw("Need to pass on the length of the vc");
  }
  
  $self->_vmap($vmap);

  my $id="virtual_contig_$un";
  $self->id($id);
  $self->primary_id($id);
  $self->display_id($id);
  $self->accession_number($id);
  $self->length($length);

  if (defined $clone) {
      $self->_clone_map($clone);
  }

  my $make = $self->SUPER::_initialize;

# set stuff in self from @args
 return $make; # success - we hope!
}

=head2 display_id

 Title   : display_id
 Usage   : $obj->display_id($newval)
 Function: get/set method for the dna id
 Returns : value of display_id
 Args    : newvalue (optional) 

=cut

sub display_id {
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'_display_id'} = $value;
    }
    return $self->{'_display_id'};
} 

=head2 primary_id

 Title   : primary_id
 Usage   : $obj->primary_id($newval);
 Function: get/set method for the primary id
 Returns : value of primary id
 Args    : newvalue (optional)


=cut

sub primary_id {
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'_primary_id'} = $value;
    }
    return $self->{'_primary_id'};
} 

=head2 accession_number

 Title   : accession_number
 Usage   : $obj->accession_number($newval)
 Function: get/set method for the dna id
 Returns : value of accession_number
 Args    : newvalue (optional) 

=cut

sub accession_number {
    my ($self) = @_;
    return $self->display_id;
}

=head2 seq

 Title   : seq
 Usage   : $string    = $obj->seq()
 Function: Returns the sequence as a string of letters.
 Returns : A scalar
 Args    : none

=cut

sub seq {
   my ($self) = @_;
   
   my $seq = $self->_seq_cache();
   
   if( defined $seq ) {
       return $seq;
   }
   
   # we have to move across the map, picking up the sequences,
   # truncating them and then adding them into the final product.
   my @map_contigs=$self->_vmap->get_all_MapContigs;
   
   my $seq_string;
   my $last_point = 1;
   
   # if there is a left overhang, add it 
   
   if( $self->_vmap->left_overhang() > 0 ) {
       $seq_string = 'N' x $self->_vmap->left_overhang();
   }
   
   #Goes through each MapContig
   foreach my $mc ( @map_contigs ) {
       my $tseq = $mc->contig->primary_seq();
       
       if( $mc->start != ($last_point+1) ) {
	   
           # Tony: added a throw here - if we get negative numbers of inserted N's
	   
	   my $no = $mc->start - $last_point;
	   
           if ($no < 0){
	       $self->throw("Error. Trying to insert negative number ($no) of N\'s into contig sequence");
           }
	   
	   $seq_string .= 'N' x $no;
	   $last_point += $no;
       } 
       
       my $trunc;
       my $end;
       
       if( $self->_clone_map == 1 ) {
	   $end = $mc->contig->length;
       } else {
	   if($mc->rightmost_end) {
	       $end = $mc->rightmost_end;
	   } else {
	       if( $mc->orientation == 1 ) {
		   $end = $mc->contig->golden_end;
	       } else {
		   $end = $mc->contig->golden_start;
	       }
	   }
       }
       if( $mc->orientation == 1 ) {
	   $trunc = $tseq->subseq($mc->start_in,$end);
       } else {
	   my $subseq = $tseq->subseq($end,$mc->start_in);
	   $trunc = $self->revcom($subseq);
       }
       $seq_string .= $trunc;
       $last_point += length($trunc);
   }
   
   # if there is a right overhang, add it 
   if( $self->_vmap->right_overhang() > 0 ) {
       $seq_string .= 'N' x $self->_vmap->right_overhang();
   }
   
   $self->_seq_cache($seq);
   return $seq_string;
}

=head2 subseq

 Title   : subseq
 Usage   : $substring = $obj->subseq(10,40);
 Function: returns the subseq from start to end, where the first base
           is 1 and the number is inclusive, ie 1-2 are the first two
           bases of the sequence
           Start cannot be larger than end but can be equal
 Returns : a string
 Args    : start and end scalars

=cut

sub subseq{
   my ($self,$start,$end) = @_;
   
   if( $start > $end ){
       $self->throw("in subseq, start [$start] has to be greater than end [$end]");
   }
   if( $start <= 0 || $end > $self->length ) {
       $self->throw("You have to have start positive and length less than the total length of sequence");
   }
   
   my ($start_rc,$start_rc_pos)=$self->_vmap->vcpos_to_rcpos($start);
   my ($end_rc,$end_rc_pos)=$self->_vmap->vcpos_to_rcpos($end);
   my $mc;
   my $start_gap=0;
   my $end_gap=0;

   if ($start_rc ne 'N') {
       print STDERR "START: contig ".$start_rc->id." and pos. $start_rc_pos\n";
       $mc=$self->_vmap->get_MapContig($start_rc->id);  
   }
   else {
       print STDERR "Start vc position in gap...\n";
       $start_gap=1;
   }

   if ($end_rc ne 'N') {
       print STDERR "END: contig ".$end_rc->id." and pos. $end_rc_pos\n";
   }
   else {
       print STDERR "End vc position in gap...\n";
       $end_gap=1;
   }

   #Let's deal straight away with the simplest case, i.e. start and end in 
   #the same RawContig. This is a very powerful way to reduce the time it 
   #takes to do small subseqs (for example for exons!)
   
   #If start and end RawContig identical, just do a subseq and complement it 
   #if the orientation of the RawContig in the VirtualContig is -1
   my $subseq;
   if ((!$start_gap) && (!$end_gap)) {
       if ($start_rc->id eq $end_rc->id) {
	   print STDERR "Using the new fast VirtualPrimarySeq method to retrieve sequence!\n";
	   
	   if ($mc->orientation == 1) {
	       my $seq=$start_rc->primary_seq->subseq($start_rc_pos,$end_rc_pos);
	       return $seq;
	   }
	   else {
	       my $seq=$start_rc->primary_seq->subseq($end_rc_pos,$start_rc_pos);
	       #return $seq;
	       return $self->revcom($seq);
	   }
       }
   }
       
   #If the start and end RawContig are different, then it gets a bit more complicated...
   #else {

       #First of all we get the seq from the start in the start RawContig to its golden_end
       #Again, if the orientation is negative we do it the other way around and revcom it!
       
       #if ($mc->orientation == 1) {
	   #$subseq=$start_rc->primary_seq->subseq($start_rc_pos,$start_rc->golden_end);
       #}
       #else {
	   #my $seq=$start_rc->primary_seq->subseq($mc->contig->golden_end,$start_rc_pos);
	   #$subseq=$self->revcom($seq);
       #}

       #Then loop through each MapContig (note: they are given back sorted by start in vc)
       #my $before_start=1;
       #foreach my $mc ($self->get_all_MapContigs) {
	   
	   #If this MapContig is before our start contig, skip it
	   #if (($mc->contig->id ne $start_rc->id)&& ($before_start=1)){
	       #next;
	   #}

	   #Then find the start contig...
	   #elsif ($mc->contig->id eq $start_rc->id) {
	       #$before_start=0;
	       #next;
	   #}
	   
	   #...and go through the rest of the contigs
	   
           #If we find the end contig, we add the last piece of 
           #subseq, and get out of the loop
	   #elsif ($mc->contig->id eq $end_rc->id) {
	       #if ($mc->orientation == 1) {
		   #$subseq.=$end_rc->primary_seq->subseq($end_rc_pos,$end_rc->golden_end);
	       #}
	       #else {
		   #my $seq=$end_rc->primary_seq->subseq($end_rc->golden_end,$end_rc_pos);
		   #$subseq.=$self->revcom($seq);
	       #} 
	       #last;
	   #}
	   
	   #If it is one of the intermediate contigs, add its whole seq to the subseq
	   #else {
	       #if ($mc->orientation == 1) {
		   #$subseq.=$mc->contig->primary_seq->subseq($mc->contig->golden_start,$mc->contig->golden_end);
	       #}
	       #else {
		   #my $seq=$mc->contig->primary_seq->subseq($mc->contig->golden_end,$mc->contig->golden_start);
		   #$subseq.=$self->revcom($seq);
	      #} 
	   #}
       #}
   #} 
   
   #Need to do this properly...
   $start--;
   return substr $self->seq, $start, ($end-$start);
   #return $subseq;
}

=head2 moltype

 Title   : moltype
 Usage   : if( $obj->moltype eq 'dna' ) { /Do Something/ }
 Function: Returns the type of sequence 
 Returns : dna
 Args    : none


=cut

sub moltype{
   return "dna";
}

=head2 id

 Title   : id
 Usage   : $id = $seq->id()
 Function: maps to display id
 Returns : display id
 Args    : none


=cut

sub id {
   my ($self)= @_;

   return $self->display_id();
}

=head2 length

 Title   : length
 Usage   : $len = $seq->length()
 Function: Returns the length of the sequence
 Returns : scalar
 Args    : none


=cut

sub length {
    my ($self,$value)= @_;

    if( defined $value) {
	$self->{'_length'} = $value;
    }
    return $self->{'_length'};    
}

=head2 can_call_new

 Title   : can_call_new
 Usage   : if( $obj->can_call_new ) {
             $newobj = $obj->new( %param );
	 }
 Function: indicates that this object can call the ->new method
 Example :
 Returns : 1 or 0
 Args    :


=cut

sub can_call_new{

   return 0;
}

=head2 _seq_cache

 Title   : _seq_cache
 Usage   : $obj->_seq_cache
 Function: get/set method for the seq cache
 Returns : value of _seq_cache
 Args    : new value (optional)

=cut

sub _seq_cache{
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'_seq_cache'} = $value;
    }
    return $self->{'_seq_cache'};
} 

=head2 _vmap

 Title   : _vmap
 Usage   : $obj->_vmap($newval);
 Function: get/set method for the virtual map
 Returns : Bio::EnsEMBL::DB::VirtualMap object
 Args    : newvalue (optional)


=cut

sub _vmap{
    my ($obj,$value) = @_;
    
    if( defined $value) {
	$obj->{'_vmap'} = $value;
    }
    return $obj->{'_vmap'};    
}

=head2 _clone_map

 Title   : _clone_map
 Usage   : $obj->_clone_map($newval);
 Function: get/set method for _clone_map
 Returns : value of _clone_map
 Args    : newvalue (optional)


=cut

sub _clone_map{
    my ($obj,$value) = @_;
    
    if( defined $value) {
	$obj->{'_clone_map'} = $value;
    }
    return $obj->{'_clone_map'};    
}

=head2 revcom

 Title   : revcom
 Usage   : $obj->revcom($newval);
 Function: get/set method for revcom
 Returns : value of revcom
 Args    : newvalue (optional)


=cut

sub revcom{
    my ($self,$str)=@_;
    $str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
    my $revcom = CORE::reverse $str;
    return $revcom;
}
1;
