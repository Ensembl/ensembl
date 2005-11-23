#
# Ensembl module for Bio::EnsEMBL::DBSQL::SequenceAdaptor
#
# Cared for by Arne Stabenau <stabenau@ebi.ac.uk>
#
# Copyright Ensembl
#
# You may distribute this module under the same terms as perl itself

=head1 NAME

Bio::EnsEMBL::DBSQL::SequenceAdaptor - produce sequence strings from locations

=head1 SYNOPSIS

$seq_adptr = $database_adaptor->get_SequenceAdaptor();
$dna = ${$seq_adptr->fetch_by_Slice_start_end_strand($slice, 1, 1000, -1);}

=head1 DESCRIPTION

An adaptor for the retrieval of DNA sequence from the EnsEMBL database

=head1 CONTACT

Post questions/comments to the EnsEMBL development list:
ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::SequenceAdaptor;

use vars qw(@ISA @EXPORT);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate);
use Bio::EnsEMBL::Utils::Sequence  qw(reverse_comp);
use Bio::EnsEMBL::Utils::Cache;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

our $SEQ_CHUNK_PWR   = 18; # 2^18 = approx. 250KB
our $SEQ_CACHE_SZ    = 5;
our $SEQ_CACHE_MAX   = (2 ** $SEQ_CHUNK_PWR) * $SEQ_CACHE_SZ;

@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});

=head2 new

  Arg [1]    : none
  Example    : my $sa = $db_adaptor->get_SequenceAdaptor();
  Description: Constructor.  Calls superclass constructor and initialises
               internal cache structure.
  Returntype : Bio::EnsEMBL:DBSQL::SequenceAdaptor
  Exceptions : none
  Caller     : DBAdaptor::get_SequenceAdaptor
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  # use an LRU cache to limit the size
  my %seq_cache;
  tie(%seq_cache, 'Bio::EnsEMBL::Utils::Cache', $SEQ_CACHE_SZ);

  $self->{'seq_cache'} = \%seq_cache;

  return $self;
}


=head2 fetch_by_Slice_start_end_strand

  Arg  [1]   : Bio::EnsEMBL::Slice slice
               The slice from which you want the sequence
  Arg  [2]   : (optional) int startBasePair 
               The start base pair relative to the start of the slice. Negative
               values or values greater than the length of the slice are fine.
               default = 1
  Arg  [3]   : (optional) int endBasePair
               The end base pair relative to the start of the slice. Negative
               values or values greater than the length of the slice are fine,
               but the end must be greater than or equal to the start
               count from 1
               default = the length of the slice
  Arg  [4]   : (optional) int strand 
               1, -1
               default = 1
  Example    : $dna = $seq_adptr->fetch_by_Slice_start_end_strand($slice, 1, 
                                                                  1000, -1);
  Description: retrieves from db the sequence for this slice
               uses AssemblyMapper to find the assembly
  Returntype : string 
  Exceptions : endBasePair should be less or equal to length of slice 
  Caller     : Bio::EnsEMBL::Slice::seq(), Slice::subseq() 
  Status     : Stable

=cut

sub fetch_by_Slice_start_end_strand {
   my ( $self, $slice, $start, $end, $strand ) = @_;

   if(!ref($slice) || !$slice->isa("Bio::EnsEMBL::Slice")) {
     throw("Slice argument is required.");
   }

   $start = 1 if(!defined($start));

   $end = $slice->end() - $slice->start() + 1 if(!defined($end));

   throw("Start must be less than or equal to end.") if($start > $end);

   $strand ||= 1;

   #get a new slice that spans the exact region to retrieve dna from
   my $right_expand  = $end - $slice->length(); #negative is fine
   my $left_expand   = 1 - $start; #negative is fine

   if($right_expand || $left_expand) {
     $slice = $slice->expand($left_expand, $right_expand);
   }

   #retrieve normalized 'non-symlinked' slices
   #this allows us to support haplotypes and PARs
   my $slice_adaptor = $slice->adaptor();
   my @symproj=@{$slice_adaptor->fetch_normalized_slice_projection($slice)};

   if(@symproj == 0) {
     throw('Could not retrieve normalized Slices. Database contains ' .
           'incorrect assembly_exception information.');
   }

   #call this method again with any slices that were 'symlinked' to by this
   #slice
   if(@symproj != 1 || $symproj[0]->[2] != $slice) {
     my $seq;
     foreach my $segment (@symproj) {
       my $symlink_slice = $segment->[2];
       #get sequence from each symlinked area
       $seq .= ${$self->fetch_by_Slice_start_end_strand($symlink_slice,
                                                        1,undef,1)};
     }
     if($strand == -1) {
       reverse_comp(\$seq);
     }
     return \$seq;
   }

   # we need to project this slice onto the sequence coordinate system
   # even if the slice is in the same coord system, we want to trim out
   # flanking gaps (if the slice is past the edges of the seqregion)
   my $csa = $self->db->get_CoordSystemAdaptor();
   my $seqlevel = $csa->fetch_sequence_level();

   my @projection=@{$slice->project($seqlevel->name(), $seqlevel->version())};

   my $seq = '';
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

     $tmp_seq = ${$self->_fetch_seq($seq_region_id,
                                    $seq_slice->start, $seq_slice->length())};

     #reverse compliment on negatively oriented slices
     if($seq_slice->strand == -1) {
       reverse_comp(\$tmp_seq);
     }

     $seq .= $tmp_seq;

     $total = $end;
   }

   #check for any remaining gaps at the end
   my $gap = $slice->length - $total;
   if($gap) {
     $seq .= 'N' x $gap;
   }

   #if the sequence is too short it is because we came in with a seqlevel
   #slice that was partially off of the seq_region.  Pad the end with Ns
   #to make long enough
   if(length($seq) != $slice->length()) {
     $seq .= 'N' x ($slice->length() - length($seq));
   }

   #if they asked for the negative slice strand revcomp the whole thing
   reverse_comp(\$seq) if($strand == -1);

   return \$seq;
}



sub _fetch_seq {
  my $self          = shift;
  my $seq_region_id = shift;
  my $start         = shift;
  my $length           = shift;

  my $cache = $self->{'seq_cache'};

  if($length < $SEQ_CACHE_MAX) {
    my $chunk_min = ($start-1) >> $SEQ_CHUNK_PWR;
    my $chunk_max = ($start + $length - 1) >> $SEQ_CHUNK_PWR;

    # piece together sequence from cached component parts

    my $entire_seq = undef;
    for(my $i = $chunk_min; $i <= $chunk_max; $i++) {
      if($cache->{"$seq_region_id:$i"}) {
        $entire_seq .= $cache->{"$seq_region_id:$i"};
      } else {
        # retrieve uncached portions of the sequence

        my $sth = $self->prepare
          ("SELECT SUBSTRING( d.sequence, ?, ?)
            FROM dna d
            WHERE d.seq_region_id = ?");

        my $tmp_seq;

        my $min = ($i << $SEQ_CHUNK_PWR) + 1;

	$sth->bind_param(1,$min,SQL_INTEGER);
	$sth->bind_param(2,1 << $SEQ_CHUNK_PWR,SQL_INTEGER);
	$sth->bind_param(3,$seq_region_id,SQL_INTEGER);

        $sth->execute();
        $sth->bind_columns(\$tmp_seq);
        $sth->fetch();
        $sth->finish();

        # always give back uppercased sequence so it can be properly softmasked
        $entire_seq .= uc($tmp_seq);
        $cache->{"$seq_region_id:$i"} = $tmp_seq;
      }
    }

    # return only the requested portion of the entire sequence
    my $min =  ($chunk_min    << $SEQ_CHUNK_PWR) + 1;
    my $max = ($chunk_max+1) << $SEQ_CHUNK_PWR;
    my $seq = substr($entire_seq, $start-$min, $length);

    return \$seq;
  } else {

    # do not do any caching for requests of very large sequences
    my $sth = $self->prepare
      ("SELECT SUBSTRING( d.sequence, ?, ?)
        FROM dna d
        WHERE d.seq_region_id = ?");

    my $tmp_seq;

    $sth->bind_param(1,$start,SQL_INTEGER);
    $sth->bind_param(2,$length,SQL_INTEGER);
    $sth->bind_param(3,$seq_region_id,SQL_INTEGER);

    $sth->execute();
    $sth->bind_columns(\$tmp_seq);
    $sth->fetch();
    $sth->finish();

    # always give back uppercased sequence so it can be properly softmasked
    $tmp_seq = uc($tmp_seq);

    return \$tmp_seq;
  }
}


=head2 store

  Arg [1]    : int $seq_region_id the id of the sequence region this dna
               will be associated with.
  Arg [2]    : string $sequence the dna sequence to be stored 
               in the database.  Note that the sequence passed in will be
               converted to uppercase.
  Example    : $seq_adaptor->store(11, 'ACTGGGTACCAAACAAACACAACA');
  Description: stores a dna sequence in the databases dna table and returns the
               database identifier for the new record.
  Returntype : none
  Exceptions : throw if the database insert fails
  Caller     : sequence loading scripts
  Status     : Stable

=cut

sub store {
  my ($self, $seq_region_id, $sequence) = @_;

  if(!$seq_region_id) {
    throw('seq_region_id is required');
  }

  $sequence = uc($sequence);

  my $statement = 
    $self->prepare("INSERT INTO dna(seq_region_id, sequence) VALUES(?,?)");

  $statement->bind_param(1,$seq_region_id,SQL_INTEGER);
  $statement->bind_param(2,$sequence,SQL_LONGVARCHAR);
  $statement->execute();

  $statement->finish();

  return;
}




=head2 fetch_by_assembly_location

  Description: DEPRECATED use fetch_by_Slice_start_end_strand() instead.

=cut

sub fetch_by_assembly_location {
   my ( $self, $chrStart, $chrEnd, 
        $strand, $chrName, $assemblyType ) = @_;

   deprecate('Use fetch_by_Slice_start_end_strand() instead');

   my $csa = $self->db->get_CoordSystem();
   my $top_cs = @{$csa->fetch_all};

   my $slice_adaptor = $self->db->get_SliceAdaptor();
   my $slice = $slice_adaptor->fetch_by_region($top_cs->name(), $chrName,
                                               $chrStart, $chrEnd,
                                               $strand, $top_cs->version);

   return $self->fetch_by_Slice_start_end_strand($slice,1, $slice->length,1);
}


=head2 fetch_by_RawContig_start_end_strand

  Description: DEPRECATED use fetch_by_Slice_start_end_strand instead

=cut

sub fetch_by_RawContig_start_end_strand {
  deprecate('Use fetch_by_Slice_start_end_strand instead.');
  fetch_by_Slice_start_end_strand(@_);
}




1;
