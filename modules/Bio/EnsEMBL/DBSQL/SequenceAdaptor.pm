#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::SequenceAdaptor
#
# Cared for by Arne Stabenau <stabenau@ebi.ac.uk>
#
# Copyright EMBL/EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::SequenceAdaptor - produce sequence strings from locations

=head1 SYNOPSIS

$seq_adptr = $database_adaptor->get_SequenceAdaptor();
$dna = $seq_adptr->fetch_by_RawContig_start_end_strand($contig, 1, 1000, -1);

=head1 DESCRIPTION

An adaptor for the retrieval of sequences of DNA from the database

=head1 CONTACT

Arne Stabenau - stabenau@ebi.ac.uk
Elia Stupka - elia@fugu-sg.org

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


package Bio::EnsEMBL::DBSQL::SequenceAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate);


@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 fetch_by_RawContig_start_end_strand

  Description: DEPRECATED use fetch_by_Slice_start_end_strand instead

=cut

sub fetch_by_RawContig_start_end_strand {
  deprecate('Use fetch_by_Slice_start_end_strand instead.');
  fetch_by_Slice_start_end_strand(@_);
}




=head2 fetch_by_Slice_start_end_strand

  Arg  [1]   : Bio::EnsEMBL::Slice slice
               The slice from which you want the sequence
  Arg  [2]   : int startBasePair 
               count from 1
  Arg  [3]   : int endBasePair 
               count from 1, -1 is last one
  Arg  [4]   : int strand 
               1, -1
  Example    : $dna = $seq_adptr->fetch_by_Slice_start_end_strand($slice, 1, 
                                                                  1000, -1);
  Description: retrieves from db the sequence for this slice
               uses AssemblyMapper to find the assembly
  Returntype : string 
  Exceptions : endBasePair should be less or equal to length of slice 
  Caller     : Bio::EnsEMBL::Slice::seq(), Slice::subseq() 

=cut

sub fetch_by_Slice_start_end_strand {
   my ( $self, $slice, $start, $end, $strand ) = @_;

   my $seq= "";

   if(!$slice || !ref($slice) || !$slice->isa("Bio::EnsEMBL::Slice")) {
     throw("Slice argument is required.");
   }

   if( $end == -1 ) {
     $end = $slice->end() - $slice->start() + 1;
   }

   if($start > $end) {
     throw("Start must be less than or equal to end.");
   }

   $strand ||= 1;

   #get a new slice that spans the exact region to retrieve dna from
   my $right_expand  = $end - $slice->length(); #negative is fine
   my $left_expand   = 1 - $start; #negative is fine

   $slice = $slice->expand($left_expand, $right_expand);

   # we need to project this slice onto the sequence coordinate system
   # if it is not already in it

   my $csa = $self->db->get_CoordSystemAdaptor();
   my $seqlevel = $csa->fetch_sequence_level();

   my @projection;
   if($slice->coord_system->equals($seqlevel)) {
     #create a fake projection to the entire length of the same slice
     @projection = ([1, $slice->length, $slice]);
   } else {
     @projection = @{$slice->project($seqlevel->name(), $seqlevel->version())};
   }

   my $slice_adaptor = $self->db()->get_SliceAdaptor();
   my $sth = $self->prepare(
               "SELECT SUBSTRING( d.sequence, ?, ?)
                FROM dna d
                WHERE d.seq_region_id = ?");
   my $seq;
   my $total = 0;
   my $tmp_seq;

   #fetch sequence from each of the sequence regions projected onto
   foreach my $segment (@projection) {
     my ($start, $end, $seq_slice) = @$segment;

     #check for gaps between segments and pad them with Ns
     my $gap = $start - $total - 1;
     if($gap) {
       $seq .= 'N' x $gap;
     }

     my $seq_region_id = $slice_adaptor->get_seq_region_id($seq_slice);

     $sth->execute($seq_slice->start, $seq_slice->length(), $seq_region_id);
     $sth->bind_columns(\$tmp_seq);
     $sth->fetch();

     #reverse compliment on negatively oriented slices
     if($seq_slice->strand == -1) {
       $self->_reverse_comp(\$tmp_seq);
     }

     $seq .= $tmp_seq;

     $total = $end;
   }

   $sth->finish();

   #check for any remaining gaps at the end
   my $gap = $slice->length - $total;
   if($gap) {
     $seq .= 'N' x $gap;
   }

   #if they asked for the negative slice strand revcomp the whole thing
   $self->_reverse_compliment(\$seq) if($strand == -1);

   return \$seq;
}



=head2 fetch_by_assembly_location

  Description: DEPRECATED use fetch_by_Slice_start_end_strand() instead.

=cut

sub fetch_by_assembly_location {
   my ( $self, $chrStart, $chrEnd, 
        $strand, $chrName, $assemblyType ) = @_;

   deprecate('Use fetch_by_Slice_start_end_strand() instead');

   my $slice_adaptor = $self->db->get_SliceAdaptor();
   my $slice = $slice_adaptor->fetch_by_region('toplevel', $chrName,
                                               $chrStart, $chrEnd,
                                               $strand);

   return $self->fetch_by_Slice_start_end_strand($slice,1, $slice->length,1);
}



=head2 store

  Arg [1]    : string $sequence the dna sequence to be stored in the database
  Arg [2]    : string $date create date to be associated with the dna sequence
               to be stored.
  Example    : $dbID = $seq_adaptor->store('ACTGGGTACCAAACAAACACAACA', $date);
  Description: stores a dna sequence in the databases dna table and returns the
               database identifier for the new record.
  Returntype : int
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::RawContigAdaptor::store

=cut

sub store {
  my ($self, $seq_region_id, $sequence, $date) = @_;

  if(!$seq_region_id) {
    throw('seq_region_id is required');
  }

  $sequence = uc($sequence);

  my $statement = $self->prepare(
        "INSERT INTO dna(seq_region_id, sequence,created) " .
        "VALUES(?, ?, FROM_UNIXTIME(?))");

  my $rv = $statement->execute($seq_region_id, $sequence, $date);
  $self->throw("Failed to insert dna $sequence") unless $rv;

  my $id = $statement->{'mysql_insertid'};

  $statement->finish();

  return $id;
}




=head2 _reverse_comp

  Arg [1]    : reference to a string $seqref
  Example    : $self->_reverse_comp(\$seqref);
  Description: Does an in place reverse compliment of a passed in string
               reference.  The string is passed by reference
               rather than by value for memory efficiency.
  Returntype : none
  Exceptions : none
  Caller     : internal

=cut

sub _reverse_comp {
  my $self = shift;
  my $seqref = shift;

  $$seqref = reverse( $$seqref );
  $$seqref =~
    tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;

  return undef;
}


1;
