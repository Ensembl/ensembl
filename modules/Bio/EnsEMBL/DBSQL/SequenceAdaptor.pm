=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::SequenceAdaptor - produce sequence strings from locations

=head1 SYNOPSIS

  my $sa = $registry->get_adaptor( 'Human', 'Core', 'Sequence' );

  my $dna =
    ${ $sa->fetch_by_Slice_start_end_strand( $slice, 1, 1000, -1 ) };

=head1 DESCRIPTION

An adaptor for the retrieval of DNA sequence from the EnsEMBL database

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::SequenceAdaptor;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw deprecate);
use Bio::EnsEMBL::Utils::Sequence  qw(reverse_comp);
use Bio::EnsEMBL::Utils::Scalar qw( assert_ref );
use DBI qw/:sql_types/;
use base qw(Bio::EnsEMBL::DBSQL::BaseAdaptor Bio::EnsEMBL::DBSQL::BaseSequenceAdaptor);
our @EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});

=head2 new

  Arg [1]    : none
  Example    : my $sa = $db_adaptor->get_SequenceAdaptor();
  Description: Constructor. Calls superclass constructor and initialises
               internal cache structure.
  Returntype : Bio::EnsEMBL::DBSQL::SequenceAdaptor
  Exceptions : none
  Caller     : DBAdaptor::get_SequenceAdaptor
  Status     : Stable

=cut

sub new {
  my ($caller, $db, $chunk_power, $cache_size) = @_;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new($db);
  $self->_init_seq_instance($chunk_power, $cache_size);
  $self->_populate_seq_region_edits();
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
  Description: Retrieves from db the sequence for this slice
               uses AssemblyMapper to find the assembly
  Returntype : string 
  Exceptions : endBasePair should be less or equal to length of slice 
  Caller     : Bio::EnsEMBL::Slice::seq(), Slice::subseq() 
  Status     : Stable

=cut

sub fetch_by_Slice_start_end_strand {
   my ( $self, $slice, $start, $end, $strand ) = @_;

   if(!ref($slice) || !($slice->isa("Bio::EnsEMBL::Slice") or $slice->isa('Bio::EnsEMBL::LRGSlice')) ) {
     throw("Slice argument is required.");
   }

   $start = 1 if(!defined($start));

 
   if ( ( !defined($end) || $start > $end || $start < 0 || $end < 0 || $slice->start> $slice->end ) && $slice->is_circular ) {
         
       if ( !defined($end) || ($start > $end ) ) {
     return $self->_fetch_by_Slice_start_end_strand_circular( $slice, $start, $end, $strand );
       }

       if ( defined($end) && ($end < 0) ) {
     $end += $slice->seq_region_length;
       }
       
       if ($start < 0) {
           $start += $slice->seq_region_length;
       }

       if($slice->start> $slice->end) {
           return $self->_fetch_by_Slice_start_end_strand_circular( $slice, $slice->start, $slice->end, $strand );
       }
  }
        
  if ( ( !defined($end) ) && (not $slice->is_circular) ) {
           $end = $slice->end() - $slice->start() + 1;
  }

  if ( $start > $end ) {
      throw("Start must be less than or equal to end.");
  }

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

     if(!defined $tmp_seq) {
         throw('No sequence found for seq_region '.$seq_region_id.':'.$seq_slice->start);
     }

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

   if(defined($self->{_rna_edits_cache}) and defined($self->{_rna_edits_cache}->{$slice->get_seq_region_id})){
     $self->_rna_edit($slice,\$seq);
   }

   #if they asked for the negative slice strand revcomp the whole thing
   reverse_comp(\$seq) if($strand == -1);

   return \$seq;
}

=head2 can_access_Slice

  Description : Returns 1 since we can access any Slice's data

=cut

sub can_access_Slice {
  return 1;
}

sub _fetch_by_Slice_start_end_strand_circular {
  my ( $self, $slice, $start, $end, $strand ) = @_;

  assert_ref( $slice, 'Bio::EnsEMBL::Slice' );
  
  $strand ||= 1;
  if ( !defined($start) ) {
    $start ||= 1;
  }

  if ( !defined($end) ) {
      $end = $slice->end() - $slice->start() + 1;
  }

  if ( $start > $end && $slice->is_circular() ) {
    my ($seq, $seq1, $seq2);

    my $midpoint = $slice->seq_region_length - $slice->start + 1;
    $seq1 = ${ $self->_fetch_by_Slice_start_end_strand_circular( $slice, 1,  $midpoint, 1 )};
    $seq2 = ${ $self->_fetch_by_Slice_start_end_strand_circular( $slice, $midpoint + 1, $slice->length(), 1 )};

    $seq = $slice->strand > 0 ? "$seq1$seq2" : "$seq2$seq1";

    reverse_comp( \$seq ) if ( $strand == -1 );

    return \$seq;
  }

  # Get a new slice that spans the exact region to retrieve dna from
  my $right_expand = $end - $slice->length();    #negative is fine
  my $left_expand  = 1 - $start;                 #negative is fine

  if ( $right_expand || $left_expand ) {
    $slice =
        $slice->strand > 0
      ? $slice->expand( $left_expand,  $right_expand )
      : $slice->expand( $right_expand, $left_expand );
  }

  # Retrieve normalized 'non-symlinked' slices.  This allows us to
  # support haplotypes and PARs.
  my $slice_adaptor = $slice->adaptor();
  my @symproj =
    @{ $slice_adaptor->fetch_normalized_slice_projection($slice) };

  if ( @symproj == 0 ) {
    throw(   'Could not retrieve normalized Slices. Database contains '
           . 'incorrect assembly_exception information.' );
  }

  # Call this method again with any slices that were 'symlinked' to by
  # this slice.
  if ( @symproj != 1 || $symproj[0]->[2] != $slice ) {
    my $seq;
    foreach my $segment (@symproj) {
      my $symlink_slice = $segment->[2];

      # Get sequence from each symlinked area.
      $seq .= ${
        $self->fetch_by_Slice_start_end_strand( $symlink_slice, 1,
                                                undef, 1 ) };
    }
    if ( $strand == -1 ) {
      reverse_comp( \$seq );
    }

    return \$seq;
  }

  # We need to project this slice onto the sequence coordinate system
  # even if the slice is in the same coord system, we want to trim out
  # flanking gaps (if the slice is past the edges of the seqregion).
  my $csa      = $self->db->get_CoordSystemAdaptor();
  my $seqlevel = $csa->fetch_sequence_level();

  my @projection =
    @{ $slice->project( $seqlevel->name(), $seqlevel->version() ) };

  my $seq   = '';
  my $total = 0;
  my $tmp_seq;

  # Fetch sequence from each of the sequence regions projected onto.
  foreach my $segment (@projection) {
    my ( $start, $end, $seq_slice ) = @{$segment};

    # Check for gaps between segments and pad them with Ns
    my $gap = $start - $total - 1;
    if ($gap) {
      $seq .= 'N' x $gap;
    }

    my $seq_region_id = $slice_adaptor->get_seq_region_id($seq_slice);

    $tmp_seq = ${
      $self->_fetch_seq( $seq_region_id, $seq_slice->start(),
                         $seq_slice->length() ) };

    # Reverse compliment on negatively oriented slices.
    if ( $seq_slice->strand == -1 ) {
      reverse_comp( \$tmp_seq );
    }

    $seq .= $tmp_seq;

    $total = $end;
  }

  # Check for any remaining gaps at the end.
  my $gap = $slice->length() - $total;

  if ($gap) {
    $seq .= 'N' x $gap;
  }

  # If the sequence is too short it is because we came in with a
  # seqlevel slice that was partially off of the seq_region.  Pad the
  # end with Ns to make long enough
  if ( length($seq) != $slice->length() ) {
    $seq .= 'N' x ( $slice->length() - length($seq) );
  }

  if ( defined( $self->{_rna_edits_cache} )
       && defined(
            $self->{_rna_edits_cache}->{ $slice->get_seq_region_id } ) )
  {
    $self->_rna_edit( $slice, \$seq );
  }

  return \$seq;
} ## end sub _fetch_by_Slice_start_end_strand_circular



=head2 _rna_edit

  Description : Performs within sequence region editting when 
                the underlying sequence is incorrect. Used by LRGs.

=cut

sub _rna_edit {
  my $self  = shift;
  my $slice = shift;
  my $seq   = shift; #reference to string

  my $s_start = $slice->start;   # substr start at 0 , but seq starts at 1 (so no -1 here)
  my $s_end = $s_start+length($$seq) - 1; # But we do need -1 here to keep the coords closed

  foreach my $edit (@{$self->{_rna_edits_cache}->{$slice->get_seq_region_id}}){
    # seq_region_attrib for an exit looks like
    # <start> <end> <edit>
    # 2568 2569 GT
    my ($start, $end, $txt) = split (/\s+/, $edit);

    # check that RNA edit is not outside the requested region : happens quite often with LRG regions
    next if ($end < $s_start);
    next if ($s_end < $start);
    # Length of the edit ($txt)
    my $edit_length = length($txt);

    # If the edit isn't fully encompassed by the slice, we need to extract the
    # edit's sub-sequence that we're patching on to the sequence
    if($start < $s_start || $end > $s_end) {
	my $edit_offset;
	# Find the offset for the start of the edit, eg.
	#      seq slice:   TC
	#      LRG edit :  CG
	# Offset should be 1 to extract starting at the G
	# Formula: offset = max(s_start - start, 0)
	$edit_offset = ($s_start - $start) < 0 ? 0 : $s_start - $start;

	# Find the length of the overlapping piece of the edit.
	# This needs to take in to account the offset
	#      seq slice:   TC
	#      LRG edit :    GA
	# Where we're want the length to be 1 from an existing offset of 1 to
	# extract just G
	# Forumla: length = length(edit) - offset - max(end - s_end, 0)
	$edit_length = length($txt) - $edit_offset - ($end - $s_end < 0 ? 0 : $end - $s_end);

	# Extract the overlapping edit
	$txt = substr($txt, $edit_offset, $edit_length);
    }

    # Apply the patch, we don't want negative offsets as that's totally wrong.
    # Do a max(start - s_start, 0) to prevent this
    substr($$seq,($start-$s_start < 0 ? 0 : $start-$s_start),$edit_length, $txt);
  }
  return;
}

=head2 _fetch_raw_seq

  Description : Communicates with the database to fetch back sequence

=cut

sub _fetch_raw_seq {
  my ($self, $id, $start, $length) = @_;
  my $sql = <<'SQL';
SELECT UPPER(SUBSTR(d.sequence, ?, ?))
FROM dna d
WHERE d.seq_region_id =?
SQL
  my $seq = $self->dbc()->sql_helper()->execute_single_result(
    -SQL => $sql, 
    -PARAMS => [[$start, SQL_INTEGER], [$length, SQL_INTEGER], [$id, SQL_INTEGER]],
    -NO_ERROR => 1
  );
  return \$seq;
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

=head2 remove   
   
  Arg [1]    : int $seq_region_id the id of the sequence region this dna   
               is associated with.   
  Example    : $seq_adaptor->remove(11);   
  Description: removes a dna sequence for a given seq_region_id    
  Returntype : none    
  Exceptions : throw if the database delete fails    
  Caller     : Internal    
  Status     : Stable    
   
=cut   
   
sub remove {   
  my ($self, $seq_region_id) = @_;   
   
  if(!$seq_region_id) {    
    throw('seq_region_id is required');    
  }    
   
  my $statement =    
    $self->prepare("DELETE FROM dna WHERE seq_region_id = ?");   
   
  $statement->bind_param(1,$seq_region_id,SQL_INTEGER);    
  $statement->execute();   
   
  $statement->finish();    
   
  return;    
}    

=head2 _populate_seq_region_edits

  Description:  Query the database for any _rna_edit attributes attached to a seq region

=cut

sub _populate_seq_region_edits {
  my ($self) = @_;
  my $sql;
  my @params = ('_rna_edit');
  if($self->db()->is_multispecies()) {
    $sql = <<'SQL';
select sra.seq_region_id, sra.value 
from seq_region_attrib sra 
join attrib_type using (attrib_type_id) 
join seq_region s using (seq_region_id)
join coord_system cs using (coord_system_id)
where code =?
and species_id =?
SQL
    push(@params, $self->db()->species_id());
  }
  else {
    $sql = <<'SQL';
select sra.seq_region_id, sra.value 
from seq_region_attrib sra join attrib_type using (attrib_type_id) 
where code = ?
SQL
  }

  my $mapper = sub {
    my ($row, $array) = @_;
    my ($seq_region_id, $value) = @{$row};
    my ($start, $end, $substring) = split (/\s+/, $value);
    my $edit_length = ($end - $start) + 1;
    my $substring_length = length($substring);
    if($edit_length != $substring_length) {
      throw "seq_region_id $seq_region_id has an attrib of type '_rna_edit' (value '$value'). Edit length ${edit_length} is not the same as the replacement's length ${substring_length}. Please fix. We only support substitutions via this mechanism";
    }
    if(defined $array) {
      push(@{$array}, $value);
      return;
    }
    return [$value];
  };

  my $edits = $self->dbc()->sql_helper->execute_into_hash(-SQL => $sql, -PARAMS => ['_rna_edit'], -CALLBACK => $mapper);
  $self->{_rna_edits_cache} = $edits if %{$edits};
  return;
}

1;
