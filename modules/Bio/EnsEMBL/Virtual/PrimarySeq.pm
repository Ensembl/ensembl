
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

Bio::EnsEMBL::Virtual::PrimarySeq - Object that "pretends" to be a PrimarySeq without actual phyiscally having the DNA.

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
package Bio::EnsEMBL::Virtual::PrimarySeq;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::PrimarySeqI;

@ISA = qw(Bio::Root::Object Bio::PrimarySeqI);

sub new {
    my ($class,@args) = @_;
    
    my $self = {};
    bless $self,$class;
    
    my ($vmap,$id,$length) = 
	$self->_rearrange([qw( 
			       VMAP
			       ID
			       )],@args);
    
    if (!$vmap->isa('Bio::EnsEMBL::Virtual::Map')) {
	$self->throw("$vmap is not a Bio::EnsEMBL::DB::VirtualMap!");
    }
    if (! defined $id) {
	$self->throw("Need to provide a unique number for the VirtualPrimarySeq id");
    }
    $self->_vmap($vmap);
    $self->id($id);
    $self->primary_id($id);
    $self->display_id($id);
    $self->accession_number($id);
    return $self;
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
   return $self->subseq(1,$self->length);
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
   
   #print STDERR "Looking at $start..$end\n";

   if( $start > $end ){
       $self->throw("in subseq, start [$start] has to be greater than end [$end]");
   }
   if( $start <= 0) {
       $self->throw("In subseq start must be positive");
   }
   if( $end > $self->length ) {
       $self->throw("Cannot ask for more than length of sequence $end vs".$self->length);
   }


   # I have no doubt that there is an easier way of writing this.
   # believe me, I *thought* this was the easier way of writing this.
   # remember that everything can be reversed as well. ;)

   # cases to worry about are

   # starting conditons:
   # start-end in one contig (optimise this for fast retrieval), return
   # start-end in gap, return 
   # all-gap virtual contigs, return 
   # start in gap, end in first contig, return
   # start in contig, end in first gap, return
   # start in gap, end further way - going into main loop

   # main loop adds gap and contig pieces

   # end conditions:
   # end in gap before last contig
   # end in contig
   # end in gap after last contig

       
   my @mapcontigs=$self->_vmap->each_MapContig();

   my $start_contig=shift(@mapcontigs);

   if( !defined $start_contig ) {
       # all gap contig!
       return 'N' x ($end - $start +1);
   }

   while ($start_contig->end < $start) {
       $start_contig = shift(@mapcontigs);
   }
   
   # could be in the middle of a gap
   if( $start_contig->start > $end ) {
       #print STDERR "start gap\n";
       return 'N' x ($end - $start +1);
   }

   # check to see if this is simple, start-end in one contig
   if( $start >= $start_contig->start && $end <= $start_contig->end ) {
       # map straight away.
       #print STDERR "simple map\n";
       if( $start_contig->orientation == 1 ) {
	   return $start_contig->contig->primary_seq->subseq($start_contig->rawcontig_start + ($start - $start_contig->start),$start_contig->rawcontig_start + ($end - $start_contig->start));
       } else {
	   my $temp = $start_contig->contig->primary_seq->subseq($start_contig->rawcontig_end - ($end - $start_contig->start),$start_contig->rawcontig_end - ($start - $start_contig->start));
	   $temp =~ tr/ATGCNatgcn/TACGNtacgn/;
	   $temp = reverse $temp;
	   return $temp;
       }
   }

   my $seqstr = "";
       
   # ok end is > than start. See if start is actually in contig
   if( $start < $start_contig->start ) {
       #print STDERR "start in gap before contig ... honest\n";

       #nope. Got some N's to put in
       $seqstr .= 'N' x ($start_contig->start - $start);
       # now put in the rest of this contig


       #print STDERR "Looking at $end vs ",$start_contig->end,"\n";

       # of course, end could be in this contig. Bugger.
       if( $end <= $start_contig->end ) {
	   #print STDERR "start in gap before: end in contig\n";

	   if( $start_contig->orientation == 1 ) {
	       $seqstr .= $start_contig->contig->primary_seq->subseq($start_contig->rawcontig_start,$start_contig->rawcontig_start + ($end - $start_contig->start));
	       return $seqstr;
	   } else {
	       my $temp = $start_contig->contig->primary_seq->subseq($start_contig->rawcontig_end-($end - $start_contig->start),$start_contig->rawcontig_end);
	       $temp =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
	       $temp = reverse $temp;
	       $seqstr .= $temp;
	       return $seqstr;
	   }
       } else {
	   $seqstr .= $start_contig->seq;
       }
   } else {
       # start is in the middle of a contig
       # print STDERR "start in middle of a contig\n";

       if( $start_contig->orientation == 1 ) {
	   $seqstr .= $start_contig->contig->primary_seq->subseq($start_contig->rawcontig_start + ($start - $start_contig->start),$start_contig->rawcontig_end);
       } else {
	   my $temp = $start_contig->contig->primary_seq->subseq($start_contig->rawcontig_start,$start_contig->rawcontig_end - ($start - $start_contig->start));
	   $temp =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
	   $temp = reverse $temp;
	   $seqstr .= $temp;
       }
   }
	   

   
   my $previous = $start_contig;
   my $current;
   #print STDERR "About to enter loop...\n";
   while( ($current = shift @mapcontigs) ) {
       #print STDERR "Looking at ",$current->end," vs ",$end,"\n";
       if( $end <= $current->end ) {
	   last;
       }
       
       # check to add Ns
       if( $previous->end+1 != $current->start ) {
	   $seqstr .= 'N' x ($current->start - $previous->end -1);
       }
       
       # add sequence
       $seqstr .= $current->seq;
       $previous = $current;
   }
   # end of sequence
   if( !defined $current ) {
	   # last contig was beside a gap. Sneaky
	   $seqstr .= 'N' x ($end - $previous->end);
	   return $seqstr;
	   
   }

   # last contig
   
   # end could be not inside, in which case only add N's
   if( $end < $current->start ) {
       $seqstr .= 'N' x ($end - $previous->end);
   } else {
       # ok. Add the remainder

       # could be that there are some N's to add first
       if( $previous->end+1 != $current->start ) {
	   $seqstr .= 'N' x ($current->start - $previous->end -1);
       }

       # add in remainder of this contig
       if( $current->orientation == 1 ) {
	   $seqstr .= $current->contig->primary_seq->subseq($current->rawcontig_start,$current->rawcontig_start + ($end - $current->start));
       } else {
	   #print STDERR "$end vs",$current->start,"\n";
	   
	   my $temp = $current->contig->primary_seq->subseq($current->rawcontig_end - ($end - $current->start),$current->rawcontig_end);
	   $temp =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
	   $temp = reverse $temp;
	   $seqstr .= $temp;
       }
   }

   return $seqstr;
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
    my ($self)= @_;
    return $self->_vmap->length();

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


=head2 _revcom

 Title   : _revcom
 Usage   : $obj->_revcom();



=cut
    
sub _revcom{
    my ($self,$str)=@_;
    $str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
    my $revcom = CORE::reverse $str;
    return $revcom;
}
1;
