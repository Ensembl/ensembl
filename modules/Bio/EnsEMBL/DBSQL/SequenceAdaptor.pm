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

SequenceAdaptor - produce sequence strings from locations

=head1 SYNOPSIS



=head1 DESCRIPTION


=head1 CONTACT

e-mail: elia@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


# Let the code begin...
package Bio::EnsEMBL::DBSQL::SequenceAdaptor;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object
use Bio::EnsEMBL::DBSQL::Adaptor;


@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);



=head2 fetch_by_contig_id_start_end

  Arg  1    : int rawContigdbID
  Arg  2    : int startBasePair
  Arg  3    : int endBasePair
    a -1 means until the end
  Function  : retrieves the dna string from the database from the 
              given RawContig internal id.
  Returntype: txt
  Exceptions: endBasePair should be less or equal than length of contig
  Caller    : Bio::EnsEMBL::RawContig::seq(), RawContig::subseq()

=cut

sub fetch_by_contig_id_start_end {
  my ( $self, $contig_id, $start, $end ) = @_;
  my $sth;

  if( $end == -1 ) { 
    $sth = $self->prepare( "SELECT SUBSTRING( sequence, $start )
                            FROM dna d, contig c 
                            WHERE d.dna_id = c.dna_id 
                            AND c.contig_id = $contig_id" );

   $sth->execute(); 
   
   my($subseq) = $sth->fetchrow
     or $self->throw("Could not fetch substr of contig " .$contig_id );
   

}


=head2 fetch_by_Slice_start_end

  Arg  1    : Bio::EnsEMBL::Slice slice
              The slice from which you want the sequence
  Arg  2    : int startBasePair 
              count from 1
  Arg  3    : int endBasePair 
              count from 1, -1 is last one
  Function  : retrieves from db the sequence for this slice
              uses AssemblyMapper to find the assembly
  Returntype: txt
  Exceptions: endBasePair should be less or equal to length of slice
  Caller    : Bio::EnsEMBL::Slice::seq(), Slice::subseq()

=cut


sub fetch_by_Slice_start_end {
   my ( $self, $slice, $start, $end ) = @_;


}


=head2 fetch_by_assembly_location

  Arg [1]   : none, txt, int, Bio::EnsEMBL::Example formal_parameter_name
    Additional description lines
    list, listref, hashref
  Function  : testable description
  Returntype: none, txt, int, float, Bio::EnsEMBL::Example
  Exceptions: none
  Caller    : object::methodname or just methodname
  Example   :  ( optional )

=cut

sub fetch_by_assembly_location {
   my $self = shift;
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
       $self->throw("in subseq, start [$start] cannot be greater than end [$end]");
   }

   my $add = 0;
   if( $end > $self->length ) {
       print STDERR ("TROUBLE - $end greater than length ".$self->length);
       $add = $end-$self->length;
       $end = $self->length;
   }

   #if( $start <= 0 || $end > $self->length ) {
   #    $self->throw("You have to have start positive and length less than the total length of sequence - calling $start:$end vs".$self->length);
   #}
   
   my $id=$self->dna_id;
   my $length= $end-$start+1;
   
   my $sth=$self->db_handle->prepare("SELECT SUBSTRING(sequence,$start,$length) FROM dna WHERE dna_id = $id");
   $sth->execute(); 
   
   my($subseq) = $sth->fetchrow
       or $self->throw("Could not fetch substr of dna " .$id);
   
   $subseq .= 'N' x $add;
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

sub  length {
   my ($self)= @_;

   if( defined $self->_length() ) {
       return $self->_length();
   }

   my $id=$self->dna_id;
   
   my $sth=$self->db_handle->prepare("SELECT length(sequence) FROM dna WHERE dna_id = $id");
   $sth->execute(); 
   
   my($length) = $sth->fetchrow
       or $self->throw("Could not determine length of dna " .$id);

   $self->_length($length);

   return $length;
}

=head2 _length

 Title   : _length
 Usage   : $obj->_length($newval)
 Function: 
 Returns : value of _length
 Args    : newvalue (optional)


=cut

sub _length{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_length'} = $value;
    }
    return $obj->{'_length'};

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

=head2 desc

 Title   : desc
 Usage   : $obj->desc($newval)
 Function: 
 Example : 
 Returns : value of desc
 Args    : newvalue (optional)


=cut

sub desc {
   my ($self,$value) = @_;
   if( defined $value && $value ne '' ) {
       $self->{'desc'} = $value;
   } 
   return $self->{'desc'} || '';
}

1;
