
#
# BioPerl module for DBPrimarySeq
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

DBPrimarySeq - A lightweight DB connected PrimarySeq Object

=head1 SYNOPSIS

First get a DBPrimarySeq object for a RawContig:
my $seq=$contig->DB_Primary_Seq;

You can then apply all the usual PrimarySeq methods to this object
my $seq->seq;

=head1 DESCRIPTION

This object is intended as a lightweight alternative to the standard 
PrimarySeq object, it is not meant as a replacement for the usual PrimarySeq.
It is primarily used when the sequence to be dealt with is too large to be 
held in memory efficiently. All the methods of PrimarySeqI are implemented, 
but rather than holding data in mmemory, they chat to the mysql database to 
get the data.

=head1 CONTACT

e-mail: elia@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::DBPrimarySeq;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;


@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my($dna_id,$dbh) =
      $self->_rearrange([qw(DNA 
                            DB_HANDLE
			    )],
			@args);

  my $make = $self->SUPER::_initialize;

  if( !defined $dna_id) {
      $self->throw("You must provide a dna id to create a DBPrimarySeq object!");
  }
  if( !defined $dbh) {
      $self->throw("You must provide a database handle to create a DBPrimarySeq object!");
  }
  
  $self->dna_id($dna_id);
  $self->db_handle($dbh);

# set stuff in self from @args
 return $make; # success - we hope!
}

=head2 dna_id

 Title   : dna_id
 Usage   : $obj->dna_id($newval)
 Function: get/set method for the dna id
 Example : 
 Returns : value of dna_id
 Args    : newvalue (optional) 

=cut

sub dna_id{
    my ($obj,$value) = @_;
    
    if( defined $value) {
	$obj->{'dna_id'} = $value;
    }
    return $obj->{'dna_id'};
}

=head2 contig_id

 Title   : contig_id
 Usage   : $obj->contig_id($newval)
 Function: get/set method for the contig internal id
 Example : 
 Returns : value of contig internal id
 Args    : newvalue (optional) 
=cut

sub contig_id{
    my ($self,$value) = @_;
    
    my $id=$self->dna_id;
    
    my $sth=$self->db_handle->prepare("SELECT internal_id FROM contig WHERE dna = $id");
    $sth->execute(); 
    
    my($contig_id) = $sth->fetchrow
	or $self->throw("No contig id for dna " .$id);
    
    return $contig_id;
}

=head2 display_id

 Title   : display_id
 Usage   : $id_string = $obj->display_id();
 Function: returns the EnsEMBL contig id
           note that most other id methods are mapped to this one
 Returns : A string
 Args    : None


=cut

sub display_id {
   my ($self) = @_;

   my $id=$self->contig_id;
   
   my $sth=$self->db_handle->prepare("SELECT id FROM contig WHERE internal_id = $id");
   $sth->execute(); 
   
   my($display_id) = $sth->fetchrow
       or $self->throw("No EnsEMBL id for contig internal id " .$id);

   return $display_id;
}

=head2 db_handle

 Title   : db_handle
 Usage   : $obj->_db_handle($newval)
 Function: get/set methof for the db_handle
 Returns : value of _db_handle
 Args    : newvalue (optional) 

=cut

sub db_handle{
   my ($self,$value) = @_;
   
   if( defined $value) {
       $self->{'_db_handle'} = $value;
   }
   return $self->{'_db_handle'};
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
   
   my $id=$self->dna_id;
   
   my $sth=$self->db_handle->prepare("SELECT sequence FROM dna WHERE id = $id");
   $sth->execute();
   my($str) = $sth->fetchrow
       or $self->throw("No DNA sequence for dna id " . $id);
   
   return $str;
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

   # remove one from start, and then length is end-start
   
   $start--;

   my $id=$self->dna_id;
   my $length= $end-$start;
   
   my $sth=$self->db_handle->prepare("SELECT SUBSTRING(sequence,$start,$length) FROM dna WHERE id = $id");
   $sth->execute(); 
   
   my($subseq) = $sth->fetchrow
       or $self->throw("Could not fetch substr of dna " .$id);
   
   return $subseq;
}

=head2 moltype

 Title   : moltype
 Usage   : if( $obj->moltype eq 'dna' ) { /Do Something/ }
 Function: Returns the type of sequence 
 Returns : dna
 Args    : none


=cut

sub moltype{
   my ($self,@args) = @_;

   return "dna";
}

=head2 id

 Title   : id
 Usage   : $id = $seq->id()
 Function: maps to display id
 Returns : display id
 Args    : none


=cut

sub  id {
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

sub  length {
   my ($self)= @_;

   my $id=$self->dna_id;
   
   my $sth=$self->db_handle->prepare("SELECT length(sequence) FROM dna WHERE id = $id");
   $sth->execute(); 
   
   my($length) = $sth->fetchrow
       or $self->throw("Could not determine length of dna " .$id);
   
   return $length;
}
