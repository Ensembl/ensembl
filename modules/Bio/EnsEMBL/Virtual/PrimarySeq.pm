
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
package Bio::EnsEMBL::DB::VirtualPrimarySeq;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::PrimarySeqI;

@ISA = qw(Bio::Root::Object Bio::PrimarySeqI);

sub new {
    my ($class,$vmap,$length,$id) = @_;
    
    my $self = {};
    bless $self,$class;

    if (!$vmap->isa('Bio::EnsEMBL::Virtual::Map')) {
	$self->throw("$vmap is not a Bio::EnsEMBL::DB::VirtualMap!");
    }
    if (! defined $id) {
	$self->throw("Need to provide a unique number for the VirtualPrimarySeq id");
    }
    if (! defined $length) {
	$self->throw("Need to pass on the length of the vc");
    }
  
    $self->_vmap($vmap);

    $self->id($id);
    $self->primary_id($id);
    $self->display_id($id);
    $self->accession_number($id);
    $self->length($length);

    if (defined $clone) {
	$self->_clone_map($clone);
    }
    
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

   $self->throw("Needs to be reimplemented in subseq");

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

   $self->throw("Not implemented yet!");

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
