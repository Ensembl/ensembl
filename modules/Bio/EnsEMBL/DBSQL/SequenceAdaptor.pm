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
use Bio::EnsEMBL::DBSQL::BaseAdaptor;


@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);



=head2 fetch_by_contig_id_start_end_strand

  Arg  1    : int rawContigdbID
  Arg  2    : int startBasePair
  Arg  3    : int endBasePair
    a -1 means until the end
  Arg  4    : int strand
              -1, 1 are possible values
  Function  : retrieves the dna string from the database from the 
              given RawContig internal id.
  Returntype: txt
  Exceptions: endBasePair should be less or equal than length of contig
              returns undef if query fails
  Caller    : Bio::EnsEMBL::RawContig::seq(), RawContig::subseq()

=cut


  
sub fetch_by_contig_id_start_end_strand {
  my ( $self, $contig_id, $start, $end, $strand ) = @_;
  my $sth;
  
  if( $start < 1 ) {
    $self->throw( "Wrong parameters" );
  }

  my $query;

  if( $end == -1 ) { 

    $query = "SELECT c.length, SUBSTRING( d.sequence, $start )
                FROM dna d, contig c 
               WHERE d.dna_id = c.dna_id 
                 AND c.contig_id = $contig_id";
      
  } else {
    my $length = $end - $start + 1;
    if( $length < 1 ) {
      $self->throw( "Wrong parameters" );
    }
    
    $query = "SELECT c.length, SUBSTRING( d.sequence, $start, $length )
                FROM dna d, contig c 
               WHERE d.dna_id = c.dna_id 
                 AND c.contig_id = $contig_id";    
  }

  if( defined $self->db()->dnadb() ) {
    $sth = $self->db()->dnadb()->prepare( $query );
  } else {
    $sth = $self->prepare( $query );
  }


  $sth->execute();

  if( my $aref = $sth->fetchrow_arrayref() ) {
    my ( $length, $seq ) = @$aref;
    $seq =~ s/\s//g;
    if( $strand == -1 ) {
      return $self->_reverse_comp( $seq );
    } else {
      return $seq;
    }
  } else {
    return undef;
  }
}

  


=head2 fetch_by_Slice_start_end_strand

  Arg  1    : Bio::EnsEMBL::Slice slice
              The slice from which you want the sequence
  Arg  2    : int startBasePair 
              count from 1
  Arg  3    : int endBasePair 
              count from 1, -1 is last one
  Arg  4    : int strand 
              1, -1
  Function  : retrieves from db the sequence for this slice
              uses AssemblyMapper to find the assembly
  Returntype: txt
  Exceptions: endBasePair should be less or equal to length of slice
  Caller    : Bio::EnsEMBL::Slice::seq(), Slice::subseq()

=cut


sub fetch_by_Slice_start_end_strand {
   my ( $self, $slice, $start, $end, $strand ) = @_;

   my $seq= "";

   if( !$slice ){
     $self->throw("need a slice to work\n");
   }

   unless 
     ($slice->isa("Bio::EnsEMBL::Slice")) {
       $self->throw("$slice isn't a slice");
     }
   
   if( $end == -1 ) {
     $end = $slice->chr_end() - $slice->chr_start() + 1;
   }

   # need to check the strand'edness of the slice as this
   # affects the direction in which the dna seq is grabbed
   if ( $slice->strand == 1 ) {
     $seq = $self->fetch_by_assembly_location
       (
	$slice->chr_start()+$start-1,
	$slice->chr_start()+$end-1,
	$strand,
	$slice->chr_name(),
	$slice->assembly_type() 
       );
   }
   elsif ( $slice->strand == -1 ) {
     $seq = $self->fetch_by_assembly_location
       (
	$slice->chr_end()-$end+1,
	$slice->chr_end()-$start+1,
	$strand,
	$slice->chr_name(),
	$slice->assembly_type() 
       );
   }
   else {
     $self->throw("Incorrect strand set on slice $slice");
   }
   return $seq;
}



=head2 fetch_by_assembly_location

  Arg   1   : int $chrStart
  Arg   2   : int $chrEnd
  Arg   3   : int $strand
  Arg   4   : txt $chrName
  Arg   5   : txt $assemblyType
  Function  : retrieve specified sequence from db. Using AssemblyMapper. Gaps are 
              filled with N
  Returntype: txt
  Exceptions: Wrong parameters give undef as result
  Caller    : general, fetch_by_Slice_start_end_strand

=cut


sub fetch_by_assembly_location {
   my ( $self, $chrStart, $chrEnd, 
        $strand, $chrName, $assemblyType ) = @_;

   my $mapper = $self->db->get_AssemblyMapperAdaptor->fetch_by_type($assemblyType);
   # $mapper->register_region($chrName,$chrStart,$chrEnd);
   
   my @coord_list = $mapper->map_coordinates_to_rawcontig
     ( $chrName, $chrStart, $chrEnd, $strand );
   
   # for each of the pieces get sequence
   my $seq = "";
   for my $segment ( @coord_list ) {
     if( $segment->isa( "Bio::EnsEMBL::Mapper::Coordinate" )) {

       my $contig_seq = $self->fetch_by_contig_id_start_end_strand
	 ( $segment->id(),
	   $segment->start(),
	   $segment->end(),
	   $segment->strand() );

       $seq .= $contig_seq;
	   
     } else {
       # its a gap
       my $length = $segment->end() - $segment->start() + 1;
       $seq .= "N" x $length;
     }
   }
   
   return $seq;
}



=head2 _reverse_comp

  Arg  1    : txt $dna_sequence
  Function  : build reverse complement string
  Returntype: txt
  Exceptions: none
  Caller    : private to this module

=cut

sub _reverse_comp {
  my $self = shift;
  my $seq = shift;
  
  $_ = reverse( $seq );
  tr/CGTAcgta/GCATgcat/;
  return $_;
}




1;
