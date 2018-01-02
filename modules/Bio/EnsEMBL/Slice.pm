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

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Slice - Arbitary Slice of a genome

=head1 SYNOPSIS

  $sa = $db->get_SliceAdaptor;

  $slice =
    $sa->fetch_by_region( 'chromosome', 'X', 1_000_000, 2_000_000 );

  # get some attributes of the slice
  my $seqname = $slice->seq_region_name();
  my $start   = $slice->start();
  my $end     = $slice->end();

  # get the sequence from the slice
  my $seq = $slice->seq();

  # get some features from the slice
  foreach $gene ( @{ $slice->get_all_Genes } ) {
    # do something with a gene
  }

  foreach my $feature ( @{ $slice->get_all_DnaAlignFeatures } ) {
    # do something with dna-dna alignments
  }

=head1 DESCRIPTION

A Slice object represents a region of a genome.  It can be used to retrieve
sequence or features from an area of interest.

NOTE: The Slice is defined by its Strand, but normal behaviour for get_all_*
methods is to return Features on both Strands. 

=head1 METHODS

=cut

package Bio::EnsEMBL::Slice;
use vars qw(@ISA);
use strict;

use Bio::PrimarySeqI;


use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning stack_trace_dump);
use Bio::EnsEMBL::RepeatMaskedSlice;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::ProjectionSegment;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Iterator;
use Bio::EnsEMBL::DBSQL::MergedAdaptor;

use Bio::EnsEMBL::StrainSlice;
#use Bio::EnsEMBL::IndividualSlice;
#use Bio::EnsEMBL::IndividualSliceFactory;
use Bio::EnsEMBL::Mapper::RangeRegistry;
use Bio::EnsEMBL::SeqRegionSynonym;
use Scalar::Util qw(weaken isweak);

# use Data::Dumper;

my $registry = "Bio::EnsEMBL::Registry";

@ISA = qw(Bio::PrimarySeqI);


=head2 new

  Arg [...]  : List of named arguments
               Bio::EnsEMBL::CoordSystem COORD_SYSTEM
               string SEQ_REGION_NAME,
               int    START,
               int    END,
               int    SEQ_REGION_LENGTH, (optional)
               string SEQ (optional)
               int    STRAND, (optional, defaults to 1)
               Bio::EnsEMBL::DBSQL::SliceAdaptor ADAPTOR (optional)
  Example    : $slice = Bio::EnsEMBL::Slice->new(-coord_system => $cs,
                                                 -start => 1,
                                                 -end => 10000,
                                                 -strand => 1,
                                                 -seq_region_name => 'X',
                                                 -seq_region_length => 12e6,
                                                 -adaptor => $slice_adaptor);
  Description: Creates a new slice object.  A slice represents a region
               of sequence in a particular coordinate system.  Slices can be
               used to retrieve sequence and features from an area of
               interest in a genome.

               Coordinates start at 1 and are inclusive.  Negative
               coordinates or coordinates exceeding the length of the
               seq_region are permitted.  Start must be less than or equal.
               to end regardless of the strand.

               Slice objects are immutable. Once instantiated their attributes
               (with the exception of the adaptor) may not be altered.  To
               change the attributes a new slice must be created.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : throws if start, end, coordsystem or seq_region_name not specified or not of the correct type
  Caller     : general, Bio::EnsEMBL::SliceAdaptor
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  #new can be called as a class or object method
  my $class = ref($caller) || $caller;

  my ($seq, $coord_system, $seq_region_name, $seq_region_length,
      $start, $end, $strand, $adaptor, $empty) =
        rearrange([qw(SEQ COORD_SYSTEM SEQ_REGION_NAME SEQ_REGION_LENGTH
                      START END STRAND ADAPTOR EMPTY)], @_);

  #empty is only for backwards compatibility
  if ($empty) {
    deprecate(   "Creation of empty slices is no longer needed"
               . "and is deprecated" );
     my $self = bless( { 'empty' => 1 }, $class );
    $self->adaptor($adaptor);
    return $self;
  }

  if ( !defined($seq_region_name) ) {
    throw('SEQ_REGION_NAME argument is required');
  }
  if ( !defined($start) ) { throw('START argument is required') }
  if ( !defined($end) )   { throw('END argument is required') }

  ## if ( $start > $end + 1 ) {
  ##   throw('start must be less than or equal to end+1');
  ## }

  if ( !defined($seq_region_length) ) { $seq_region_length = $end }

  if ( $seq_region_length <= 0 ) {
    throw('SEQ_REGION_LENGTH must be > 0');
  }

  if ( defined($seq) && CORE::length($seq) != ( $end - $start + 1 ) ) {
    throw('SEQ must be the same length as the defined LENGTH not '
        . CORE::length($seq)
        . ' compared to '
        . ( $end - $start + 1 ) );
  }

  if(defined($coord_system)) {
   if(!ref($coord_system) || !$coord_system->isa('Bio::EnsEMBL::CoordSystem')){
     throw('COORD_SYSTEM argument must be a Bio::EnsEMBL::CoordSystem');
   }
   if($coord_system->is_top_level()) {
     throw('Cannot create slice on toplevel CoordSystem.');
   }
  } else {
   warning("Slice without coordinate system");
   #warn(stack_trace_dump());
  }

  $strand ||= 1;

  if($strand != 1 && $strand != -1) {
    throw('STRAND argument must be -1 or 1');
  }

  if(defined($adaptor)) {
    if(!ref($adaptor) || !$adaptor->isa('Bio::EnsEMBL::DBSQL::SliceAdaptor')) {
      throw('ADAPTOR argument must be a Bio::EnsEMBL::DBSQL::SliceAdaptor');
    }
  }

  my $self = bless {'coord_system'      => $coord_system,
                'seq'               => $seq,
                'seq_region_name'   => $seq_region_name,
                'seq_region_length' => $seq_region_length,
                'start'             => int($start),
                'end'               => int($end),
                'strand'            => $strand}, $class;

  $self->adaptor($adaptor);

  return $self;

}

=head2 new_fast

  Arg [1]    : hashref to be blessed
  Description: Construct a new Bio::EnsEMBL::Slice using the hashref.
  Exceptions : none
  Returntype : Bio::EnsEMBL::Slice
  Caller     : general
  Status     : Stable

=cut


sub new_fast {
  my $class = shift;
  my $hashref = shift;
  my $self = bless $hashref, $class;
  weaken($self->{adaptor})  if ( ! isweak($self->{adaptor}) );
  return $self;
}

=head2 adaptor

  Arg [1]    : (optional) Bio::EnsEMBL::DBSQL::SliceAdaptor $adaptor
  Example    : $adaptor = $slice->adaptor();
  Description: Getter/Setter for the slice object adaptor used
               by this slice for database interaction.
  Returntype : Bio::EnsEMBL::DBSQL::SliceAdaptor
  Exceptions : thorws if argument passed is not a SliceAdaptor
  Caller     : general
  Status     : Stable

=cut

sub adaptor{
   my $self = shift;

   if(@_) {
     my $ad = shift;
     if(defined($ad)) {
       if(!ref($ad) || !$ad->isa('Bio::EnsEMBL::DBSQL::SliceAdaptor')) {
         throw('Argument must be a Bio::EnsEMBL::DBSQL::SliceAdaptor');
       }
     }
     weaken($self->{'adaptor'} = $ad);
   }

   return $self->{'adaptor'};
}



=head2 seq_region_name

  Arg [1]    : none
  Example    : $seq_region = $slice->seq_region_name();
  Description: Returns the name of the seq_region that this slice is on. For
               example if this slice is in chromosomal coordinates the
               seq_region_name might be 'X' or '10'.

               This function was formerly named chr_name, but since slices can
               now be on coordinate systems other than chromosomal it has been
               changed.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub seq_region_name {
  my $self = shift;
  return $self->{'seq_region_name'};
}

=head2 seq_region_start

  Example    : $seq_region_start = $slice->seq_region_start();
  Description: Returns the start position of this slice relative to the
               start of the sequence region that it was created on.
               Since slices are always in genomic coordinates this is
               an alias to start()
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub seq_region_start {
    my $self = shift;
    return $self->start();
}

=head2 seq_region_end

  Example    : $seq_region_end = $slice->seq_region_end();
  Description: Returns the end position of this slice relative to the
               start of the sequence region that it was created on.
               Since slices are always in genomic coordinates this is
               an alias to end()
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub seq_region_end {
    my $self = shift;
    return $self->end();
}

=head2 seq_region_length

  Arg [1]    : none
  Example    : $seq_region_length = $slice->seq_region_length();
  Description: Returns the length of the entire seq_region that this slice is
               on. For example if this slice is on a chromosome this will be
               the length (in basepairs) of the entire chromosome.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub seq_region_length {
  my $self = shift;
  return $self->{'seq_region_length'};
}


=head2 coord_system

  Arg [1]    : none
  Example    : print $slice->coord_system->name();
  Description: Returns the coordinate system that this slice is on.
  Returntype : Bio::EnsEMBL::CoordSystem
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub coord_system {
  my $self = shift;
  return $self->{'coord_system'};
}

=head2 source

  Arg [1]    : (optional) String $value
  Example    : print $slice->source();
  Description: Returns the source this slice is coming from
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub source {
  my $self = shift;
  $self->{'source'} = shift if (@_);
  return $self->{'source'};
}

=head2 coord_system_name

  Arg [1]    : none
  Example    : print $slice->coord_system_name()
  Description: Convenience method.  Gets the name of the coord_system which
               this slice is on.
               Returns undef if this Slice does not have an attached
               CoordSystem.
  Returntype: string or undef
  Exceptions: none
  Caller    : general
  Status     : Stable

=cut

sub coord_system_name {
  my $self = shift;
  my $csystem = $self->{'coord_system'};
  return ($csystem) ? $csystem->name() : undef;
}


=head2 centrepoint

  Arg [1]    : none
  Example    : $cp = $slice->centrepoint();
  Description: Returns the mid position of this slice relative to the
               start of the sequence region that it was created on.
               Coordinates are inclusive and start at 1.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub centrepoint {
  my $self = shift;
  return ($self->{'start'}+$self->{'end'})/2;
}

=head2 start

  Arg [1]    : none
  Example    : $start = $slice->start();
  Description: Returns the start position of this slice relative to the
               start of the sequence region that it was created on.
               Coordinates are inclusive and start at 1.  Negative coordinates
               or coordinates exceeding the length of the sequence region are
               permitted.  Start is always less than or equal to end
               regardless of the orientation of the slice.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub start {
  my $self = shift;
  return $self->{'start'};
}



=head2 end

  Arg [1]    : none
  Example    : $end = $slice->end();
  Description: Returns the end position of this slice relative to the
               start of the sequence region that it was created on.
               Coordinates are inclusive and start at 1.  Negative coordinates
               or coordinates exceeding the length of the sequence region are
               permitted.  End is always greater than or equal to start
               regardless of the orientation of the slice.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub end {
  my $self = shift;
  return $self->{'end'};
}



=head2 strand

  Arg [1]    : none
  Example    : $strand = $slice->strand();
  Description: Returns the orientation of this slice on the seq_region it has
               been created on
  Returntype : int (either 1 or -1)
  Exceptions : none
  Caller     : general, invert
  Status     : Stable

=cut

sub strand{
  my $self = shift;
  return $self->{'strand'};
}





=head2 name

  Arg [1]    : none
  Example    : my $results = $cache{$slice->name()};
  Description: Returns the name of this slice. The name is formatted as a colon
               delimited string with the following attributes:
               coord_system:version:seq_region_name:start:end:strand

               Slices with the same name are equivalent and thus the name can
               act as a hash key.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub name {
  my $self = shift;

  my $cs = $self->{'coord_system'};

  return join(':',
              ($cs) ? $cs->name()    : '',
              ($cs) ? $cs->version() : '',
              $self->{'seq_region_name'},
              $self->{'start'},
              $self->{'end'},
              $self->{'strand'});
}



=head2 length

  Arg [1]    : none
  Example    : $length = $slice->length();
  Description: Returns the length of this slice in basepairs
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub length {
  my ($self) = @_;

  my $length = $self->{'end'} - $self->{'start'} + 1;
 
  if ( $self->{'start'} > $self->{'end'} && $self->is_circular() ) {    
     $length += $self->{'seq_region_length'};
  }

  return $length;
}

=head2 is_reference
  Arg        : none
  Example    : my $reference = $slice->is_reference()
  Description: Returns 1 if slice is a reference  slice else 0
  Returntype : int
  Caller     : general
  Status     : At Risk

=cut

sub is_reference {
  my ($self) = @_;

  if ( !defined( $self->{'is_reference'} ) ) {
    $self->{'is_reference'} =
      $self->adaptor()->is_reference( $self->get_seq_region_id() );
  }

  return $self->{'is_reference'};
}

=head2 is_toplevel
  Arg        : none
  Example    : my $top = $slice->is_toplevel()
  Description: Returns 1 if slice is a toplevel slice else 0
  Returntype : int
  Caller     : general
  Status     : At Risk

=cut

sub is_toplevel {
  my ($self) = @_;

  if ( !defined( $self->{'toplevel'} ) ) {
    $self->{'toplevel'} =
      $self->adaptor()->is_toplevel( $self->get_seq_region_id() );
  }

  return $self->{'toplevel'};
}

=head2 has_karyotype
  Arg        : none
  Example    : my $top = $slice->has_karyotype()
  Description: Returns 1 if slice is part of the karyotype else 0
  Returntype : int
  Caller     : general
  Status     : At Risk

=cut

sub has_karyotype {
  my ($self) = @_;

  if ( !defined( $self->{'karyotype'} ) ) {
    $self->{'karyotype'} =
      $self->adaptor()->has_karyotype( $self->get_seq_region_id() );
  }

  return $self->{'karyotype'};
}

=head2 karyotype_rank
  Arg        : none
  Example    : my $rank = $slice->karyotype_rank()
  Description: Returns the numeric ranking in the karyotype. Otherwise 0 is returned
  Returntype : int
  Caller     : general
  Status     : At Risk

=cut

sub karyotype_rank {
  my ($self) = @_;
  if(! defined( $self->{karyotype_rank})) {
    my $rank = $self->adaptor()->get_karyotype_rank($self->get_seq_region_id());
    $self->{karyotype_rank} = $rank if $rank;
  }
  return $self->{karyotype_rank} || 0;
}

=head2 is_circular
  Arg        : none
  Example    : my $circ = $slice->is_circular()
  Description: Returns 1 if slice is a circular slice else 0
  Returntype : int
  Caller     : general
  Status     : Stable

=cut

sub is_circular {
  my ($self) = @_;
  my $adaptor = $self->adaptor();
  return 0 if ! defined $adaptor;
  if (! exists $self->{'circular'}) {
    my $id = $adaptor->get_seq_region_id($self);
    $self->{circular} = $adaptor->is_circular($id);
  }
  return $self->{circular};
}

=head2 invert

  Arg [1]    : none
  Example    : $inverted_slice = $slice->invert;
  Description: Creates a copy of this slice on the opposite strand and
               returns it.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub invert {
  my $self = shift;

  # make a shallow copy of the slice via a hash copy and flip the strand
  my %s = %$self;
  $s{'strand'} = $self->{'strand'} * -1;

  # reverse compliment any attached sequence
  reverse_comp(\$s{'seq'}) if($s{'seq'});

  # bless and return the copy
  return  bless \%s, ref $self;
}



=head2 seq

  Arg [1]    : none
  Example    : print "SEQUENCE = ", $slice->seq();
  Description: Returns the sequence of the region represented by this
               slice formatted as a string.
               Only available for the default coord_system
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub seq {
  my $self = shift;

  # special case for in-between (insert) coordinates
  return '' if($self->start() == $self->end() + 1);

  return $self->{'seq'} if($self->{'seq'});

  if($self->adaptor()) {
    my $seqAdaptor = $self->adaptor()->db()->get_SequenceAdaptor();
    return ${$seqAdaptor->fetch_by_Slice_start_end_strand($self,1,undef,1)};
  }

  # no attached sequence, and no db, so just return Ns
  return 'N' x $self->length();
}



=head2 subseq

  Arg  [1]   : int $startBasePair
               relative to start of slice, which is 1.
  Arg  [2]   : int $endBasePair
               relative to start of slice.
  Arg  [3]   : (optional) int $strand
               The strand of the slice to obtain sequence from. Default
               value is 1.
  Description: returns string of dna sequence
  Returntype : txt
  Exceptions : end should be at least as big as start
               strand must be set
  Caller     : general
  Status     : Stable

=cut

sub subseq {
  my ( $self, $start, $end, $strand ) = @_;

  if ( $end+1 < $start ) {
    throw("End coord + 1 is less than start coord");
  }

  # handle 'between' case for insertions
  return '' if( $start == $end + 1);

  $strand = 1 unless(defined $strand);

  if ( $strand != -1 && $strand != 1 ) {
    throw("Invalid strand [$strand] in call to Slice::subseq.");
  }
  my $subseq;
  if($self->adaptor){
    my $seqAdaptor = $self->adaptor->db->get_SequenceAdaptor();
    $subseq = ${$seqAdaptor->fetch_by_Slice_start_end_strand
      ( $self, $start,
        $end, $strand )};
  } else {
    ## check for gap at the beginning and pad it with Ns
    if ($start < 1) {
      $subseq = "N" x (1 - $start);
      $start = 1;
    }
    $subseq .= substr ($self->seq(), $start-1, $end - $start + 1);
    ## check for gap at the end and pad it with Ns
    if ($end > $self->length()) {
      $subseq .= "N" x ($end - $self->length());
    }
    reverse_comp(\$subseq) if($strand == -1);
  }
  return $subseq;
}

=head2 sub_Slice_Iterator

  Arg[1]      : int The chunk size to request
  Example     : my $i = $slice->sub_Slice_Iterator(60000); 
                while($i->has_next()) { warn $i->next()->name(); }
  Description : Returns an iterator which batches subslices of this Slice 
                in the requested chunk size
  Returntype  : Bio::EnsEMBL::Utils::Iterator next() will return the next
                 chunk of Slice
  Exceptions  : None

=cut

sub sub_Slice_Iterator {
  my ($self, $chunk_size) = @_;
  throw "Need a chunk size to divide the slice by" if ! $chunk_size;
  my $here = 1;
  my $end = $self->length();
  my $iterator_sub = sub {
    while($here <= $end) {
      my $there = $here + $chunk_size - 1;
      $there = $end if($there > $end); 
      my $slice = $self->sub_Slice($here, $there);
      $here = $there + 1;
      return $slice;
    }
    return;
  };
  return Bio::EnsEMBL::Utils::Iterator->new($iterator_sub);
}

=head2 assembly_exception_type

  Example     : $self->assembly_exception_type(); 
  Description : Returns the type of slice this is. If it is reference then you
                will get 'REF' back. Otherwise you will get the first
                element from C<get_all_AssemblyExceptionFeatures()>. If no
                assembly exception exists you will get an empty string back.
  Returntype  : String
  Exceptions  : None
  Caller      : Public
  Status      : Beta

=cut

sub assembly_exception_type {
  my ($self) = @_;
  my $type = q{};
  if($self->is_reference()) {
    $type = 'REF';
  }
  else {
    my $assembly_exceptions = $self->get_all_AssemblyExceptionFeatures();
    if(@{$assembly_exceptions}) {
      $type = $assembly_exceptions->[0]->type();
    }
  }
  return $type;
}

=head2 is_chromosome

  Example			: print ($slice->is_chromosome()) ? 'I am a chromosome' : 'Not one'; 
  Description	: A chromosome is a slice with a karyotype rank assigned to it.
  Returntype 	: Boolean indicates if the current object is a chromosome
  Exceptions 	: None

=cut

sub is_chromosome {
  my ($self) = @_;

  return $self->has_karyotype();
}


=head2 get_base_count

  Arg [1]    : none
  Example    : $c_count = $slice->get_base_count->{'c'};
  Description: Retrieves a hashref containing the counts of each bases in the
               sequence spanned by this slice.  The format of the hash is :
               { 'a' => num,
                 'c' => num,
                 't' => num,
                 'g' => num,
                 'n' => num,
                 '%gc' => num }

               All bases which are not in the set [A,a,C,c,T,t,G,g] are
               included in the 'n' count.  The 'n' count could therefore be
               inclusive of ambiguity codes such as 'y'.
               The %gc is the ratio of GC to AT content as in:
               total(GC)/total(ACTG) * 100
               This function is conservative in its memory usage and scales to
               work for entire chromosomes.
  Returntype : hashref
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_base_count {
  my $self = shift;

  my $a = 0;
  my $c = 0;
  my $t = 0;
  my $g = 0;

  my $start = 1;
  my $end;

  my $RANGE = 100_000;
  my $len   = $self->length();

  my $seq;

  while ( $start <= $len ) {
    $end = $start + $RANGE - 1;
    $end = $len if ( $end > $len );

    $seq = $self->subseq( $start, $end );

    $a += $seq =~ tr/Aa//;
    $c += $seq =~ tr/Cc//;
    $t += $seq =~ tr/Tt//;
    $g += $seq =~ tr/Gg//;

    $start = $end + 1;
  }

  my $actg = $a + $c + $t + $g;

  my $gc_content = 0;
  if ( $actg > 0 ) {    # Avoid dividing by 0
    $gc_content = sprintf( "%1.2f", ( ( $g + $c )/$actg )*100 );
  }

  return { 'a'   => $a,
           'c'   => $c,
           't'   => $t,
           'g'   => $g,
           'n'   => $len - $actg,
           '%gc' => $gc_content };
}



=head2 project

  Arg [1]    : string $name
               The name of the coordinate system to project this slice onto
  Arg [2]    : string $version
               The version of the coordinate system (such as 'NCBI34') to
               project this slice onto
  Example    :
    my $clone_projection = $slice->project('clone');

    foreach my $segment (@$clone_projection) {
      my $clone = $segment->to_Slice();
      print $slice->seq_region_name(), ':', $segment->from_start(), '-',
            $segment->from_end(), ' -> ',
            $clone->seq_region_name(), ':', $clone->start(), '-',$clone->end(),
            ':', $clone->strand(), "\n";
    }
  Description: Returns the results of 'projecting' this slice onto another
               coordinate system.  Projecting to a coordinate system that
               the slice is assembled from is analagous to retrieving a tiling
               path.  This method may also be used to 'project up' to a higher
               level coordinate system, however.

               This method returns a listref of triplets [start,end,slice]
               which represents the projection.  The start and end defined the
               region of this slice which is made up of the third value of
               the triplet: a slice in the requested coordinate system.
  Returntype : list reference of Bio::EnsEMBL::ProjectionSegment objects which
               can also be used as [$start,$end,$slice] triplets
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub project {
  my $self = shift;
  my $cs_name = shift;
  my $cs_version = shift;

  throw('Coord_system name argument is required') if(!$cs_name);

  my $slice_adaptor = $self->adaptor();

  if(!$slice_adaptor) {
    warning("Cannot project without attached adaptor.");
    return [];
  }

  if(!$self->coord_system()) {
    warning("Cannot project without attached coord system.");
    return [];
  }


  my $db = $slice_adaptor->db();
  my $csa = $db->get_CoordSystemAdaptor();
  my $cs = $csa->fetch_by_name($cs_name, $cs_version);
  my $slice_cs = $self->coord_system();

  if(!$cs) {
    throw("Cannot project to unknown coordinate system " .
          "[$cs_name $cs_version]");
  }

  # no mapping is needed if the requested coord system is the one we are in
  # but we do need to check if some of the slice is outside of defined regions
  if($slice_cs->equals($cs)) {
    return $self->_constrain_to_region();
  }

  my @projection;
  my $current_start = 1;

  # decompose this slice into its symlinked components.
  # this allows us to handle haplotypes and PARs
  my $normal_slice_proj =
    $slice_adaptor->fetch_normalized_slice_projection($self);
  foreach my $segment (@$normal_slice_proj) {
    my $normal_slice = $segment->[2];

    $slice_cs = $normal_slice->coord_system();

    my $asma = $db->get_AssemblyMapperAdaptor();
    my $asm_mapper = $asma->fetch_by_CoordSystems($slice_cs, $cs);

    # perform the mapping between this slice and the requested system
    my @coords;

    if( defined $asm_mapper ) {
     @coords = $asm_mapper->map($normal_slice->seq_region_name(),
				 $normal_slice->start(),
				 $normal_slice->end(),
				 $normal_slice->strand(),
				 $slice_cs);
    } else {
      $coords[0] = Bio::EnsEMBL::Mapper::Gap->new( $normal_slice->start(),
						   $normal_slice->end());
    }


    # my $last_rank = 0;
    #construct a projection from the mapping results and return it
    foreach my $coord (@coords) {
      my $coord_start  = $coord->start();
      my $coord_end    = $coord->end();
      my $length       = $coord_end - $coord_start + 1;

      if ( $coord_start > $coord_end ) {
        $length =
          $normal_slice->seq_region_length() -
          $coord_start +
          $coord_end + 1;
      }

#      if( $last_rank != $coord->rank){
#	$current_start = 1;
#	print "LAST rank has changed to ".$coord->rank."from $last_rank \n";
#     }
#      $last_rank = $coord->rank;

      #skip gaps
      if($coord->isa('Bio::EnsEMBL::Mapper::Coordinate')) {

        my $coord_cs     = $coord->coord_system();

        # If the normalised projection just ended up mapping to the
        # same coordinate system we were already in then we should just
        # return the original region.  This can happen for example, if we
        # were on a PAR region on Y which refered to X and a projection to
        # 'toplevel' was requested.
        if($coord_cs->equals($slice_cs)) {
          # trim off regions which are not defined
          return $self->_constrain_to_region();
        }
	#create slices for the mapped-to coord system
        my $slice = $slice_adaptor->fetch_by_seq_region_id(
                                                    $coord->id(),
                                                    $coord_start,
                                                    $coord_end,
                                                    $coord->strand());

	my $current_end = $current_start + $length - 1;

	if ($current_end > $slice->seq_region_length() && $slice->is_circular ) {
	    $current_end -= $slice->seq_region_length();
        }

        push @projection, bless([$current_start, $current_end, $slice],
                                "Bio::EnsEMBL::ProjectionSegment");
      }

      $current_start += $length;
    }
  }

  return \@projection;
}


sub _constrain_to_region {
  my $self = shift;

  my $entire_len = $self->seq_region_length();

  #if the slice has negative coordinates or coordinates exceeding the
  #exceeding length of the sequence region we want to shrink the slice to
  #the defined region

  if($self->{'start'} > $entire_len || $self->{'end'} < 1) {
    #none of this slice is in a defined region
    return [];
  }

  my $right_contract = 0;
  my $left_contract  = 0;
  if($self->{'end'} > $entire_len) {
    $right_contract = $entire_len - $self->{'end'};
  }
  if($self->{'start'} < 1) {
    $left_contract = $self->{'start'} - 1;
  }

  my $new_slice;
  if($left_contract || $right_contract) {
      #if slice in negative strand, need to swap contracts
      if ($self->strand == 1) {
	  $new_slice = $self->expand($left_contract, $right_contract);
      }
      elsif ($self->strand == -1) {
	  $new_slice = $self->expand($right_contract, $left_contract);
      }
  } else {
    $new_slice = $self;
  }

  return [bless [1-$left_contract, $self->length()+$right_contract,
                 $new_slice], "Bio::EnsEMBL::ProjectionSegment" ];
}


=head2 expand

  Arg [1]    : (optional) int $five_prime_expand
               The number of basepairs to shift this slices five_prime
               coordinate by.  Positive values make the slice larger,
               negative make the slice smaller.
               coordinate left.
               Default = 0.
  Arg [2]    : (optional) int $three_prime_expand
               The number of basepairs to shift this slices three_prime
               coordinate by. Positive values make the slice larger,
               negative make the slice smaller.
               Default = 0.
  Arg [3]    : (optional) bool $force_expand
               if set to 1, then the slice will be contracted even in the case
               when shifts $five_prime_expand and $three_prime_expand overlap.
               In that case $five_prime_expand and $three_prime_expand will be set
               to a maximum possible number and that will result in the slice
               which would have only 2pbs.
               Default = 0.
  Arg [4]    : (optional) int* $fpref
               The reference to a number of basepairs to shift this slices five_prime
               coordinate by. Normally it would be set to $five_prime_expand.
               But in case when $five_prime_expand shift can not be applied and
               $force_expand is set to 1, then $$fpref will contain the maximum possible
               shift
  Arg [5]    : (optional) int* $tpref
               The reference to a number of basepairs to shift this slices three_prime
               coordinate by. Normally it would be set to $three_prime_expand.
               But in case when $five_prime_expand shift can not be applied and
               $force_expand is set to 1, then $$tpref will contain the maximum possible
               shift
  Example    : my $expanded_slice      = $slice->expand( 1000, 1000);
               my $contracted_slice    = $slice->expand(-1000,-1000);
               my $shifted_right_slice = $slice->expand(-1000, 1000);
               my $shifted_left_slice  = $slice->expand( 1000,-1000);
               my $forced_contracted_slice    = $slice->expand(-1000,-1000, 1, \$five_prime_shift, \$three_prime_shift);

  Description: Returns a slice which is a resized copy of this slice.  The
               start and end are moved outwards from the center of the slice
               if positive values are provided and moved inwards if negative
               values are provided. This slice remains unchanged.  A slice
               may not be contracted below 1bp but may grow to be arbitrarily
               large.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : warning if an attempt is made to contract the slice below 1bp
  Caller     : general
  Status     : Stable

=cut

sub expand {
  my $self              = shift;
  my $five_prime_shift  = shift || 0;
  my $three_prime_shift = shift || 0;
  my $force_expand      = shift || 0;
  my $fpref             = shift;
  my $tpref             = shift;

  if ( $self->{'seq'} ) {
    warning(
       "Cannot expand a slice which has a manually attached sequence ");
    return undef;
  }

  my $sshift = $five_prime_shift;
  my $eshift = $three_prime_shift;

  if ($sshift == 0 && $eshift == 0) { return $self; }

  if ( $self->{'strand'} != 1 ) {
    $eshift = $five_prime_shift;
    $sshift = $three_prime_shift;
  }

  my $new_start = $self->{'start'} - $sshift;
  my $new_end   = $self->{'end'} + $eshift;

  # Wrap around on circular slices
  if (( $new_start <= 0 || $new_start > $self->seq_region_length() || $new_end <= 0 
        || $new_end > $self->seq_region_length() ) && ( $self->is_circular() ) ) {
      
      if ( $new_start <= 0 ) {
        $new_start = $self->seq_region_length() + $new_start;
      }
      if ( $new_start > $self->seq_region_length() ) {
        $new_start -= $self->seq_region_length();
      }
  
      if ( $new_end <= 0 ) {
        $new_end = $self->seq_region_length() + $new_end;
      }
      if ( $new_end > $self->seq_region_length() ) {
        $new_end -= $self->seq_region_length();
      }      
      
  }

  if ( $new_start > $new_end  && (not $self->is_circular() ) ) {
   
      if ($force_expand) {
        # Apply max possible shift, if force_expand is set
        if ( $sshift < 0 ) {
          # if we are contracting the slice from the start - move the
          # start just before the end
          $new_start = $new_end - 1;
          $sshift    = $self->{start} - $new_start;
        }

        if ( $new_start > $new_end ) {
          # if the slice still has a negative length - try to move the
          # end
          if ( $eshift < 0 ) {
            $new_end = $new_start + 1;
            $eshift  = $new_end - $self->{end};
          }
        }
        # return the values by which the primes were actually shifted
        $$tpref = $self->{strand} == 1 ? $eshift : $sshift;
        $$fpref = $self->{strand} == 1 ? $sshift : $eshift;
      }
      if ( $new_start > $new_end ) {
        throw('Slice start cannot be greater than slice end');
      }
  }

  #fastest way to copy a slice is to do a shallow hash copy
  my %new_slice = %$self;
  $new_slice{'start'} = int($new_start);
  $new_slice{'end'}   = int($new_end);

  return bless \%new_slice, ref($self);
} ## end sub expand

=head2 constrain_to_seq_region
  Example    : $new_slice = $slice->expand(1000,10000);
               $new_slice = $new_slice->constrain_to_seq_region();
  Description: Used to prevent overly zealous expand calls going off the end of
               the sequence region. It contracts the start and end where needed
               and produces a slice copy with the tweaked coordinates.
  Returntype : Bio::EnsEMBL::Slice
=cut

sub constrain_to_seq_region {
    my ($self) = @_;
    # circular calculations should already be taken care of
    if ($self->is_circular) {return $self} 
    my $new_start = $self->start;
    my $new_end   = $self->end;
    
    my $seq_region = $self->seq_region_Slice;
    
    if ($new_start < $seq_region->start) {$new_start = $seq_region->start}
    if ($new_end > $seq_region->end) {$new_end = $seq_region->end}
    
    my %new_slice = %$self;
    $new_slice{'start'} = $new_start;
    $new_slice{'end'}   = $new_end;

    return bless \%new_slice, ref($self);
}


=head2 sub_Slice

  Arg   1    : int $start, refers to the start of the subslice relative to the input slice
  Arg   2    : int $end, refers to the end of the subslice relative to the input slice
  Arge [3]   : int $strand
  Example    : none
  Description: Makes another Slice that covers only part of this Slice
               If a Slice is requested which lies outside of the boundaries
               of this function will return undef.  This means that
               behaviour will be consistant whether or not the slice is
               attached to the database (i.e. if there is attached sequence
               to the slice).  Alternatively the expand() method or the
               SliceAdaptor::fetch_by_region method can be used instead.
  Returntype : Bio::EnsEMBL::Slice or undef if arguments are wrong
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub sub_Slice {
  my ( $self, $start, $end, $strand ) = @_;

  my $length = $self->length();

  if( $start < 1 || $start > $length ) {
    # throw( "start argument not valid" );
    return undef;
  }

  if( $end < $start || $end > $length ) {
    # throw( "end argument not valid" )
    return undef;
  }

  my ( $new_start, $new_end, $new_strand, $new_seq );
  if( ! defined $strand ) {
    $strand = 1;
  }

  if( $self->{'strand'} == 1 ) {
    $new_start = $self->{'start'} + $start - 1;
    $new_end = $self->{'start'} + $end - 1;
    $new_strand = $strand;
  } else {
    $new_start = $self->{'end'} - $end + 1;;
    $new_end = $self->{'end'} - $start + 1;
    $new_strand = -$strand;
  }

  if( defined $self->{'seq'} ) {
    $new_seq = $self->subseq( $start, $end, $strand );
  }

  #fastest way to copy a slice is to do a shallow hash copy
  my $new_slice = {%{$self}};
  $new_slice->{'start'} = int($new_start);
  $new_slice->{'end'}   = int($new_end);
  $new_slice->{'strand'} = $new_strand;
  if( $new_seq ) {
    $new_slice->{'seq'} = $new_seq;
  }
  weaken($new_slice->{adaptor});

  return bless $new_slice, ref($self);
}



=head2 seq_region_Slice

  Arg [1]    : none
  Example    : $slice = $slice->seq_region_Slice();
  Description: Returns a slice which spans the whole seq_region which this slice
               is on.  For example if this is a slice which spans a small region
               of chromosome X, this method will return a slice which covers the
               entire chromosome X. The returned slice will always have strand
               of 1 and start of 1.  This method cannot be used if the sequence
               of the slice has been set manually.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : warning if called when sequence of Slice has been set manually.
  Caller     : general
  Status     : Stable

=cut

sub seq_region_Slice {
  my $self = shift;

  if($self->{'seq'}){
    warning("Cannot get a seq_region_Slice of a slice which has manually ".
            "attached sequence ");
    return undef;
  }

  # quick shallow copy
  my $slice;
  %{$slice} = %{$self};
  bless $slice, ref($self);
  weaken($slice->{adaptor});

  $slice->{'start'}  = 1;
  $slice->{'end'}    = $slice->{'seq_region_length'};
  $slice->{'strand'} = 1;

  return $slice;
}


=head2 get_seq_region_id

  Arg [1]    : none
  Example    : my $seq_region_id = $slice->get_seq_region_id();
  Description: Gets the internal identifier of the seq_region that this slice
               is on. Note that this function will not work correctly if this
               slice does not have an attached adaptor. Also note that it may
               be better to go through the SliceAdaptor::get_seq_region_id
               method if you are working with multiple databases since is
               possible to work with slices from databases with different
               internal seq_region identifiers.
  Returntype : int or undef if slices does not have attached adaptor
  Exceptions : warning if slice is not associated with a SliceAdaptor
  Caller     : assembly loading scripts, general
  Status     : Stable

=cut

sub get_seq_region_id {
  my ($self) = @_;
  if($self->adaptor) {
    return $self->adaptor->get_seq_region_id($self);
  } 
  else {
    warning('Cannot retrieve seq_region_id without attached adaptor.');
    return undef;
  }
}

=head2 get_genome_component

  Arg []     : none
  Example    : my $genome_component = $slice->get_genome_component();
  Description: Returns the genome component of the slice
  Returntype : Scalar; the identifier of the genome component of the slice
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_genome_component {
  my $self = shift;
  return $self->adaptor->get_genome_component_for_slice($self);
}

=head2 get_all_Attributes

  Arg [1]    : optional string $attrib_code
               The code of the attribute type to retrieve values for.
  Example    : ($htg_phase) = @{$slice->get_all_Attributes('htg_phase')};
               @slice_attributes    = @{$slice->get_all_Attributes()};
  Description: Gets a list of Attributes of this slice''s seq_region.
               Optionally just get Attrubutes for given code.
  Returntype : listref Bio::EnsEMBL::Attribute
  Exceptions : warning if slice does not have attached adaptor
  Caller     : general
  Status     : Stable

=cut

sub get_all_Attributes {
  my ($self, $attrib_code) = @_;
  if(my $adaptor = $self->_get_CoreAdaptor('Attribute')) {
    my $attribs = $adaptor->fetch_all_by_Slice($self);
    if(defined $attrib_code) {
      my $uc_attrib = uc($attrib_code);
      return [ grep { uc($_->code()) eq $uc_attrib } @{$attribs}];
    }
    return $attribs;
  }
  return [];
}

=head2 get_all_PredictionTranscripts

  Arg [1]    : (optional) string $logic_name
               The name of the analysis used to generate the prediction
               transcripts obtained.
  Arg [2]    : (optional) boolean $load_exons
               If set to true will force loading of all PredictionExons
               immediately rather than loading them on demand later.  This
               is faster if there are a large number of PredictionTranscripts
               and the exons will be used.
  Example    : @transcripts = @{$slice->get_all_PredictionTranscripts};
  Description: Retrieves the list of prediction transcripts which overlap
               this slice with logic_name $logic_name.  If logic_name is
               not defined then all prediction transcripts are retrieved.
  Returntype : listref of Bio::EnsEMBL::PredictionTranscript
  Exceptions : warning if slice does not have attached adaptor
  Caller     : none
  Status     : Stable

=cut

sub get_all_PredictionTranscripts {
   my ($self,$logic_name, $load_exons) = @_;

   if(!$self->adaptor()) {
     warning('Cannot get PredictionTranscripts without attached adaptor');
     return [];
   }
   my $pta = $self->adaptor()->db()->get_PredictionTranscriptAdaptor();
   return $pta->fetch_all_by_Slice($self, $logic_name, $load_exons);
}

=head2 get_all_DnaAlignFeatures

  Arg [1]    : (optional) string $logic_name
               The name of the analysis performed on the dna align features
               to obtain.
  Arg [2]    : (optional) float $score
               The mimimum score of the features to retrieve
  Arg [3]    : (optional) string $dbtype
               The name of an attached database to retrieve the features from
               instead, e.g. 'otherfeatures'.
  Arg [4]    : (optional) float hcoverage
               The minimum hcoverage od the featurs to retrieve
  Example    : @dna_dna_align_feats = @{$slice->get_all_DnaAlignFeatures};
  Description: Retrieves the DnaDnaAlignFeatures which overlap this slice with
               logic name $logic_name and with score above $score.  If
               $logic_name is not defined features of all logic names are
               retrieved.  If $score is not defined features of all scores are
               retrieved. Strand of the Slice is not considered.
  Returntype : listref of Bio::EnsEMBL::DnaDnaAlignFeatures
  Exceptions : warning if slice does not have attached adaptor
  Caller     : general
  Status     : Stable

=cut

sub get_all_DnaAlignFeatures {
  my ($self, $logic_name, $score, $dbtype, $hcoverage) = @_;
  return $self->_get_AlignFeatures($logic_name, $score, $dbtype, $hcoverage, 'DnaAlignFeature');
}

=head2 get_all_ProteinAlignFeatures

  Arg [1]    : (optional) string $logic_name
               The name of the analysis performed on the protein align features
               to obtain.
  Arg [2]    : (optional) float $score
               The mimimum score of the features to retrieve
  Arg [3]    : (optional) string $dbtype
               The name of an attached database to retrieve features from
               instead.
  Arg [4]    : (optional) float hcoverage
               The minimum hcoverage od the featurs to retrieve
  Example    : @dna_pep_align_feats = @{$slice->get_all_ProteinAlignFeatures};
  Description: Retrieves the DnaPepAlignFeatures which overlap this slice with
               logic name $logic_name and with score above $score.  If
               $logic_name is not defined features of all logic names are
               retrieved.  If $score is not defined features of all scores are
               retrieved.
  Returntype : listref of Bio::EnsEMBL::DnaPepAlignFeatures
  Exceptions : warning if slice does not have attached adaptor
  Caller     : general
  Status     : Stable

=cut

sub get_all_ProteinAlignFeatures {
  my ($self, $logic_name, $score, $dbtype, $hcoverage) = @_;
  return $self->_get_AlignFeatures($logic_name, $score, $dbtype, $hcoverage, 'ProteinAlignFeature');
}

=head2 _get_AlignFeatures

  Arg [1]    : (optional) string $logic_name
               The name of the analysis performed on the protein align features
               to obtain.
  Arg [2]    : (optional) float $score
               The mimimum score of the features to retrieve
  Arg [3]    : (optional) string $dbtype
               The name of an attached database to retrieve features from
               instead.
  Arg [4]    : (optional) float hcoverage
               The minimum hcoverage od the featurs to retrieve
  Arg [5]    : string $align_type
               The type of adaptor to retrieve alignments from. Must be an BaseAlignFeature derived
               class
  Description: Generic method which deals with the retrieval of either AlignFeature adaptor and allows
               you to switch the adaptor values are retrieved from.

=cut

sub _get_AlignFeatures {
  my ($self, $logic_name, $score, $dbtype, $hcoverage, $align_type) = @_;
  throw "No align_type given" unless $align_type;
  if(my $adaptor = $self->_get_CoreAdaptor($align_type, $dbtype)) {
    if(defined($score) and defined ($hcoverage)){
      warning "Cannot specify score and hcoverage when querying for $align_type. Using score only";
    }
    if(defined($score)){
      return $adaptor->fetch_all_by_Slice_and_score($self,$score, $logic_name);
    }
    return $adaptor->fetch_all_by_Slice_and_hcoverage($self, $hcoverage, $logic_name);
  }
  return [];
}

=head2 get_all_SimilarityFeatures

  Arg [1]    : (optional) string $logic_name
               the name of the analysis performed on the features to retrieve
  Arg [2]    : (optional) float $score
               the lower bound of the score of the features to be retrieved
  Example    : @feats = @{$slice->get_all_SimilarityFeatures};
  Description: Retrieves all dna_align_features and protein_align_features
               with analysis named $logic_name and with score above $score.
               It is probably faster to use get_all_ProteinAlignFeatures or
               get_all_DnaAlignFeatures if a sepcific feature type is desired.
               If $logic_name is not defined features of all logic names are
               retrieved.  If $score is not defined features of all scores are
               retrieved.
  Returntype : listref of Bio::EnsEMBL::BaseAlignFeatures
  Exceptions : warning if slice does not have attached adaptor
  Caller     : general
  Status     : Stable

=cut

sub get_all_SimilarityFeatures {
  my ($self, $logic_name, $score) = @_;

  my @out = ();

  push @out, @{$self->get_all_ProteinAlignFeatures($logic_name, $score) };
  push @out, @{$self->get_all_DnaAlignFeatures($logic_name, $score) };

  return \@out;
}

=head2 get_all_SimpleFeatures

  Arg [1]    : (optional) string $logic_name
               The name of the analysis performed on the simple features
               to obtain.
  Arg [2]    : (optional) float $score
               The mimimum score of the features to retrieve
  Example    : @simple_feats = @{$slice->get_all_SimpleFeatures};
  Description: Retrieves the SimpleFeatures which overlap this slice with
               logic name $logic_name and with score above $score.  If
               $logic_name is not defined features of all logic names are
               retrieved.  If $score is not defined features of all scores are
               retrieved.
  Returntype : listref of Bio::EnsEMBL::SimpleFeatures
  Exceptions : warning if slice does not have attached adaptor
  Caller     : general
  Status     : Stable

=cut

sub get_all_SimpleFeatures {
  my ($self, $logic_name, $score, $dbtype) = @_;
  if(my $adaptor = $self->_get_CoreAdaptor('SimpleFeature', $dbtype)) {
    return $adaptor->fetch_all_by_Slice_and_score($self, $score, $logic_name);
  }
  return [];
}

=head2 get_all_RepeatFeatures

  Arg [1]    : (optional) string $logic_name
               The name of the analysis performed on the repeat features
               to obtain.
  Arg [2]    : (optional) string/array $repeat_type
               Limits features returned to those of the specified 
               repeat_type. Can specify a single value or an array reference
               to limit by more than one
  Arg [3]    : (optional) string $db
               Key for database e.g. core/vega/cdna/....
  Example    : @repeat_feats = @{$slice->get_all_RepeatFeatures(undef,'Type II Transposons')};
  Description: Retrieves the RepeatFeatures which overlap  with
               logic name $logic_name and with score above $score.  If
               $logic_name is not defined features of all logic names are
               retrieved.
  Returntype : listref of Bio::EnsEMBL::RepeatFeatures
  Exceptions : warning if slice does not have attached adaptor
  Caller     : general
  Status     : Stable

=cut

sub get_all_RepeatFeatures {
  my ($self, $logic_name, $repeat_type, $dbtype) = @_;
  if(my $adaptor = $self->_get_CoreAdaptor('RepeatFeature', $dbtype)) {
    return $adaptor->fetch_all_by_Slice($self, $logic_name, $repeat_type);
  }
  return [];
}

=head2 get_all_LD_values

    Arg [1]     : (optional) Bio::EnsEMBL::Variation::Population $population
    Description : returns all LD values on this slice. This function will only work correctly if the variation
                  database has been attached to the core database. If the argument is passed, will return the LD information
                  in that population
    ReturnType  : Bio::EnsEMBL::Variation::LDFeatureContainer
    Exceptions  : none
    Caller      : contigview, snpview
     Status     : Stable

=cut

sub get_all_LD_values {
  my $self = shift;
  my $population = shift;

  my $ld_adaptor = $self->_get_VariationAdaptor('LDFeatureContainer');
  if($ld_adaptor) {
    return $ld_adaptor->fetch_by_Slice($self,$population);
  }
  return [];
}

=head2 _get_VariationFeatureAdaptor

Shortcut method here because VariationFeature is an often requested
adaptor type.

=cut

sub _get_VariationFeatureAdaptor {
  my ($self, $dbtype) = @_;
  return $self->_get_VariationAdaptor('VariationFeature', $dbtype);
}

=head2 _get_StructuralVariationFeatureAdaptor

Shortcut method here because StructuralVariationFeature is an often requested
adaptor type.

=cut

sub _get_StructuralVariationFeatureAdaptor {
  my ($self, $dbtype) = @_;
  return $self->_get_VariationAdaptor('StructuralVariationFeature', $dbtype);
}

=head2 _get_VariationAdaptor

  Arg  [1]    : String object_type to retrieve an adaptor for
  Arg  [2]    : String dbtype to search for the given adaptor in. Defaults to variation
  Description : Searches for the specified adaptor in the Registry and returns it. Otherwise
                it will return nothing if the adaptor was not found
  ReturnType  : Bio::EnsEMBL::DBSQL::BaseAdaptor derrived instance (specific to variation)
  Exceptions  : none

=cut

sub _get_VariationAdaptor {
  my ($self, $object_type, $dbtype) = @_;
  # very important to do this defaulting since we *have* to assume the variation
  # database is called 'variation'. Using the current group will not work because
  # that will be something like 'core' (most likely), 'ensembl' or 'vega'.
  $dbtype ||= 'variation';

  # Also we do not care about Registry->get_db() for variation DBs
  my $do_not_check_db = 1;
  
  return $self->_get_Adaptor($object_type, $dbtype, $do_not_check_db);
}

=head2 _get_CoreAdaptor

  Arg  [1]    : String object_type to retrieve an adaptor for
  Arg  [2]    : String dbtype to search for the given adaptor in. Defaults to core
  Description : Searches for the specified adaptor in the Registry and returns it. Otherwise
                it will return nothing if the adaptor was not found
  ReturnType  : Bio::EnsEMBL::DBSQL::BaseAdaptor derrived instance (specific to core-like dbs)
  Exceptions  : none

=cut

sub _get_CoreAdaptor {
  my ($self, $object_type, $dbtype) = @_;
  #Simple pass through
  return $self->_get_Adaptor($object_type, $dbtype);
}

=head2 _get_Adaptor

  Arg  [1]    : String object_type to retrieve an adaptor for
  Arg  [2]    : String dbtype to search for the given adaptor in
  Arg  [3]    : Boolean Turn off the checking of Registry->get_db() for your 
                adaptor.
  Description : Searches for the specified adaptor in the Registry and returns it. Otherwise
                it will return nothing if the adaptor was not found. We consult the 
                "special" adaptors held by Bio::EnsEMBL::Registry::get_db() method and then
                fall back to the normal methods of finding an adaptor.

                This method will warn when adaptors are missing but will never through an
                exception. It is up to the calling code to decide how to handle the unavailablity
                of an adaptor.
  ReturnType  : Bio::EnsEMBL::DBSQL::BaseAdaptor derrived instance. Otherwise it returns nothing
  Exceptions  : none

=cut

sub _get_Adaptor {
  my ($self, $object_type, $dbtype, $do_not_check_db) = @_;

  if(!$object_type) {
    warning('Object type is a required parameter');
    return;
  }

  my $adaptor = $self->adaptor();
    
  if(!$adaptor) {
    warning("Cannot get a ${object_type} adaptor without a SliceAdaptor attached to this instance of ".ref($self));
    return;
  }

  my $local_db = $adaptor->db();
  my $species = $local_db->species();

  #First we query for the DBAdaptor using get_db(). This is a deprecated method
  #call but "special" adaptors can be registered via this method. We must
  #consult here 1st to find the possible special adaptor
  if(!$do_not_check_db && $dbtype) {
    my $db = $registry->get_db($local_db, $dbtype);
    if($db) {
      # If we got a return then use this DBAdaptor's species name, group and the given object type.
      # Special adaptors can have different species names
      $adaptor = $registry->get_adaptor($db->species(), $db->group(), $object_type);
    }
    else {
      #Otherwise just use the current species, dbtype and object type
      $adaptor = $registry->get_adaptor($species, $dbtype, $object_type);
    }
  }
  # Otherwise our query group is the one attached to the current adaptor
  else {
    #If not set use the group attached to the local adaptor 
    $dbtype ||= $local_db->group();
    $adaptor = $registry->get_adaptor($species, $dbtype, $object_type);
  }
  return $adaptor if $adaptor;

  warning("No adaptor could be found for the species ${species}, database type ${dbtype} and object type ${object_type}");
  return;
}

=head2 get_all_VariationFeatures

    Args [1]    : (optional) ArrayRef $so_terms
                  SequenceOntology terms to limit the fetch to
    Args [2]    : (optional) boolean $without_children
                  Do not query using the children of the given SO terms 
                  i.e. query using the given terms directly
    Args [3]    : (optional) ArrayRef $included_so 
                  ArrayRef of SequenceOntology which should be queried for
                  without children. This argument allows you to combine SO terms with children
                  from argument 1 with extra non-child SO terms. e.g. you wish to query for
                  all protein_altering_variant (specified in argument 1) variations which 
                  would be defined by child SO terms but also wanted stop_retained_variant linked variations
                  defined by this argument
    Args [4]    : (optional) string $dbtype
                  The dbtype of variation to obtain (i.e. can be different from the "variation" type).
                  This assumes that the extra db has been added to the DBAdaptor under this name (using the
                  DBConnection::add_db_adaptor method).
    Description : Returns all germline variation features on this slice. This function will 
                  only work correctly if the variation database has been attached to the core 
                  database.
                  If $so_terms is specified, only variation features with a consequence type
                  that matches or is an ontological child of any of the supplied terms will
                  be returned
    ReturnType  : listref of Bio::EnsEMBL::Variation::VariationFeature
    Exceptions  : none
    Caller      : contigview, snpview
    Status      : Stable

=cut

sub get_all_VariationFeatures{
  my ($self, $so_terms, $without_children, $included_so, $dbtype) = @_;
  deprecate('get_all_VariationFeatures is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor::fetch_all_by_Slice_SO_terms instead.');
  if (my $vf_adaptor = $self->_get_VariationFeatureAdaptor($dbtype)) {
    return $vf_adaptor->fetch_all_by_Slice_SO_terms($self, $so_terms, $without_children, $included_so);
  }
  return [];
}

=head2 get_all_somatic_VariationFeatures

    Args [1]    : (optional) ArrayRef $so_terms
                  SequenceOntology terms to limit the fetch to
    Args [2]    : (optional) boolean $without_children
                  Do not query using the children of the given SO terms 
                  i.e. query using the given terms directly
    Args [3]    : (optional) ArrayRef $included_so 
                  ArrayRef of SequenceOntology which should be queried for
                  without children. This argument allows you to combine SO terms with children
                  from argument 1 with extra non-child SO terms. e.g. you wish to query for
                  all protein_altering_variant (specified in argument 1) variations which 
                  would be defined by child SO terms but also wanted stop_retained_variant linked variations
                  defined by this argument
    Args [4]    : (optional) string $dbtype
                  The dbtype of variation to obtain (i.e. can be different from the "variation" type).
                  This assumes that the extra db has been added to the DBAdaptor under this name (using the
                  DBConnection::add_db_adaptor method).
    Description : Returns all somatic variation features on this slice. This function will only 
                  work correctly if the variation database has been attached to the core database.
    ReturnType  : listref of Bio::EnsEMBL::Variation::VariationFeature
    Exceptions  : none
    Status      : Stable

=cut

sub get_all_somatic_VariationFeatures {
  my ($self, $so_terms, $without_children, $included_so, $dbtype) = @_;
  deprecate('get_all_somatic_VariationFeatures is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor::fetch_all_somatic_by_Slice_SO_terms instead.');
  if (my $vf_adaptor = $self->_get_VariationFeatureAdaptor($dbtype)) {
    return $vf_adaptor->fetch_all_somatic_by_Slice_SO_terms($self, $so_terms, $without_children, $included_so);
  }
  return [];
}

=head2 get_all_somatic_VariationFeatures_by_source

    Arg [1]     : string $source [optional]
                  The name of the source to query for
    Arg [2]     : string $dbtype [optional]
                  The dbtype of variation to obtain (i.e. can be different from the "variation" type).
                  This assumes that the extra db has been added to the DBAdaptor under this name (using the
                  DBConnection::add_db_adaptor method).
    Description : Returns all somatic variation features, from a defined source name (e.g.'COSMIC'), 
                  on this slice. This function will only work correctly if the variation database
                  has been attached to the core database.
    ReturnType  : listref of Bio::EnsEMBL::Variation::VariationFeature
    Exceptions  : none
    Status      : Stable

=cut

sub get_all_somatic_VariationFeatures_by_source {
  my ($self, $source, $dbtype) = @_;
  deprecate('get_all_somatic_VariationFeatures_by_source is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor::fetch_all_somatic_by_Slice_Source instead.');
  my $constraint = (defined($source)) ? " s.name='$source' " : undef;
  if (my $vf_adaptor = $self->_get_VariationFeatureAdaptor($dbtype)) {
    return $vf_adaptor->fetch_all_somatic_by_Slice_constraint($self, $constraint);
  }
  return [];
}


=head2 get_all_VariationFeatures_with_phenotype

    Arg [1]     : $variation_feature_source [optional]
    Arg [2]     : $phenotype_source [optional]
    Arg [3]     : $phenotype_name [optional]
    Arg [4]     : string $dbtype [optional]
                  The dbtype of variation to obtain (i.e. can be different from the "variation" type).
                  This assumes that the extra db has been added to the DBAdaptor under this name (using the
                  DBConnection::add_db_adaptor method).
    Description : returns all germline variation features on this slice associated with a phenotype.
                  This function will only work correctly if the variation database has been
                  attached to the core database.
                  If $variation_feature_source is set only variations from that source
                  are retrieved.
                  If $phenotype_source is set only variations whose annotations come from
                  $annotation_source will be retrieved.
                  If $phenotype_name is set only variations with that annotation will be retrieved.
                  $phenotype_name can be a phenotype internal dbID.
    ReturnType  : ArrayRef of Bio::EnsEMBL::Variation::VariationFeature
    Exceptions  : none
    Caller      : contigview, snpview
    Status      : Stable

=cut

sub get_all_VariationFeatures_with_phenotype {
  my ($self, $source, $p_source, $phenotype, $dbtype) = @_;
  deprecate('get_all_VariationFeatures_with_phenotype is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor::fetch_all_with_phenotype_by_Slice instead.');
  if (my $vf_adaptor = $self->_get_VariationFeatureAdaptor($dbtype)) {
    return $vf_adaptor->fetch_all_with_phenotype_by_Slice($self, $source, $p_source, $phenotype);
  }
  return [];
}

=head2 get_all_somatic_VariationFeatures_with_phenotype

    Arg [1]     : $variation_feature_source [optional]
    Arg [2]     : $phenotype_source [optional]
    Arg [3]     : $phenotype_name [optional]
    Arg [4]     : string $dbtype [optional]
                  The dbtype of variation to obtain (i.e. can be different from the "variation" type).
                  This assumes that the extra db has been added to the DBAdaptor under this name (using the
                  DBConnection::add_db_adaptor method).
    Description : returns all somatic variation features on this slice associated with a phenotype.
                  (see get_all_VariationFeatures_with_phenotype for further documentation)
    ReturnType  : listref of Bio::EnsEMBL::Variation::VariationFeature
    Exceptions  : none
    Status      : Stable

=cut

sub get_all_somatic_VariationFeatures_with_phenotype {
  my ($self, $source, $p_source, $phenotype, $dbtype) = @_;
  deprecate('get_all_somatic_VariationFeatures_with_phenotype is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor::fetch_all_somatic_with_phenotype_by_Slice instead.');
  if (my $vf_adaptor = $self->_get_VariationFeatureAdaptor($dbtype)) {
    return $vf_adaptor->fetch_all_somatic_with_phenotype_by_Slice($self, $source, $p_source, $phenotype);
  }
  return [];
}

=head2 get_all_VariationFeatures_by_VariationSet

    Arg [1]     : Bio::EnsEMBL:Variation::VariationSet $set
    Arg [2]     : string $dbtype [optional]
                  The dbtype of variation to obtain (i.e. can be different from the "variation" type).
                  This assumes that the extra db has been added to the DBAdaptor under this name (using the
                  DBConnection::add_db_adaptor method).
    Description :returns all variation features on this slice associated with a given set.
                 This function will only work correctly if the variation database has been
                 attached to the core database. 
    ReturnType : listref of Bio::EnsEMBL::Variation::VariationFeature
    Exceptions : none
    Caller     : contigview, snpview
    Status     : Stable

=cut

sub get_all_VariationFeatures_by_VariationSet {
  my ($self, $set, $dbtype) = @_;
  deprecate('get_all_VariationFeatures_by_VariationSet is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor::fetch_all_by_Slice_VariationSet instead.');
  if (my $vf_adaptor = $self->_get_VariationFeatureAdaptor($dbtype)) {
    return $vf_adaptor->fetch_all_by_Slice_VariationSet($self, $set);  
  }
  return [];
}

=head2 get_all_StructuralVariationFeatures

    Arg [1]     : int $include_evidence [optional]
    Arg [2]     : (optional) string $dbtype
                  The dbtype of structural variation to obtain (i.e. can be different from the "variation" type).
                  This assumes that the extra db has been added to the DBAdaptor under this name (using the
                  DBConnection::add_db_adaptor method).
    Description : returns all structural variation features on this slice. This function will only work
                  correctly if the variation database has been attached to the core database.
                  By default, it only returns structural variant features which are not labelled 
                  as "CNV_PROBE".
                  If $include_evidence is set (i.e. $include_evidence=1), structural variation features from 
                  both structural variation (SV) and their supporting structural variations (SSV) will be 
                  returned. By default, it only returns features from structural variations (SV).
    ReturnType  : listref of Bio::EnsEMBL::Variation::StructuralVariationFeature
    Exceptions  : none
    Caller      : contigview, snpview, structural_variation_features
    Status      : Stable

=cut

sub get_all_StructuralVariationFeatures {
  my $self             = shift;
  deprecate('get_all_StructuralVariationFeatures is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor::fetch_all_by_Slice_constraint instead.');
  my $include_evidence = shift;
  my $dbtype           = shift;
  
  if (my $svf_adaptor = $self->_get_StructuralVariationFeatureAdaptor($dbtype)) {
    return $svf_adaptor->fetch_all_by_Slice_constraint($self,undef,$include_evidence);
  }
  return [];
}

=head2 get_all_StructuralVariationFeatures_by_size_range
    Arg [1]     : int $size_min (minimum size of the structural variant)
    Arg [2]     : int $size_max (maximum size of the structural variant) [optional]
    Arg [3]     : int $include_evidence [optional]
    Arg [4]     : (optional) string $dbtype
                  The dbtype of structural variation to obtain (i.e. can be different from the "variation" type).
                  This assumes that the extra db has been added to the DBAdaptor under this name (using the
                  DBConnection::add_db_adaptor method).
    Description : returns all structural variation features overlapping this slice with a size greater than the minimum size 
                  defined in the first argument, and (optional) lesser than the maximun size defined in the second argument. 
                  This function will only work correctly if the variation database has been attached to the core database.
                  By default, it only returns structural variant features which are not labelled 
                  as "CNV_PROBE".
                  If $include_evidence is set (i.e. $include_evidence=1), structural variation features from 
                  both structural variation (SV) and their supporting structural variations (SSV) will be 
                  returned. By default, it only returns features from structural variations (SV).
    ReturnType  : listref of Bio::EnsEMBL::Variation::StructuralVariationFeature
    Exceptions  : none
    Caller      : contigview, snpview, structural_variation_features
    Status      : Stable

=cut

sub get_all_StructuralVariationFeatures_by_size_range {
  my $self             = shift;
  deprecate('get_all_StructuralVariationFeatures_by_size_range is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor::fetch_all_by_Slice_size_range instead.');
  my $size_min         = shift;
  my $size_max         = shift;
  my $include_evidence = shift;
  my $dbtype           = shift;
  
  my $constraint = qq{svf.seq_region_end-svf.seq_region_start>=$size_min};
     $constraint .= qq{ AND svf.seq_region_end-svf.seq_region_start<$size_max } if (defined $size_max);
  
  if (my $svf_adaptor = $self->_get_StructuralVariationFeatureAdaptor($dbtype)) {
    return $svf_adaptor->fetch_all_by_Slice_constraint($self,$constraint,$include_evidence);
  }
  return [];
}


=head2 get_all_somatic_StructuralVariationFeatures_by_size_range
    Arg [1]     : int $size_min (minimum size of the structural variant)
    Arg [2]     : int $size_max (maximum size of the structural variant) [optional]
    Arg [3]     : int $include_evidence [optional]  
    Arg [4]     : (optional) string $dbtype
                  The dbtype of structural variation to obtain (i.e. can be different from the "variation" type).
                  This assumes that the extra db has been added to the DBAdaptor under this name (using the
                  DBConnection::add_db_adaptor method).
    Description : returns all the somatic structural variation features overlapping this slice with a size greater than the minimum 
                  size defined in the first argument, and (optional) lesser than the maximun size defined in the second argument. 
                  This function will only work correctly if the variation database has been attached to the core database.
                  By default, it only returns structural variant features which are not labelled 
                  as "CNV_PROBE".
                  If $include_evidence is set (i.e. $include_evidence=1), structural variation features from 
                  both structural variation (SV) and their supporting structural variations (SSV) will be 
                  returned. By default, it only returns features from structural variations (SV).
    ReturnType  : listref of Bio::EnsEMBL::Variation::StructuralVariationFeature
    Exceptions  : none
    Caller      : contigview, snpview, structural_variation_features
    Status      : Stable

=cut

sub get_all_somatic_StructuralVariationFeatures_by_size_range {
  my $self             = shift;
  deprecate('get_all_somatic_StructuralVariationFeatures_by_size_range is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor::fetch_all_somatic_by_Slice_size_range instead.');
  my $size_min         = shift;
  my $size_max         = shift;
  my $include_evidence = shift;
  my $dbtype           = shift;
    
  my $constraint = qq{svf.seq_region_end-svf.seq_region_start>=$size_min};
     $constraint .= qq{ AND svf.seq_region_end-svf.seq_region_start<$size_max } if (defined $size_max);
  
   if (my $svf_adaptor = $self->_get_StructuralVariationFeatureAdaptor($dbtype)) {
    return $svf_adaptor->fetch_all_somatic_by_Slice_constraint($self,$constraint,$include_evidence);
  }
  return [];
}
=cut


=head2 get_all_StructuralVariationFeatures_by_VariationSet

    Arg [1]     : Bio::EnsEMBL:Variation::VariationSet $set
    Arg [2]     : (optional) string $dbtype
                  The dbtype of structural variation to obtain (i.e. can be different from the "variation" type).
                  This assumes that the extra db has been added to the DBAdaptor under this name (using the
                  DBConnection::add_db_adaptor method).
    Description : returns all structural variation features on this slice associated with a 
                  given set.
                  This function will only work correctly if the variation database has been
                  attached to the core database. 
    ReturnType : listref of Bio::EnsEMBL::Variation::StructuralVariationFeature
    Exceptions : none
    Caller     : contigview, snpview
    Status     : Stable

=cut

sub get_all_StructuralVariationFeatures_by_VariationSet {
  my $self   = shift;
  deprecate('get_all_StructuralVariationFeatures_by_VariationSet is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor::fetch_all_by_Slice_VariationSet instead.');
  my $set    = shift;
  my $dbtype = shift;
  
  if (my $svf_adaptor = $self->_get_StructuralVariationFeatureAdaptor($dbtype)) {
    return $svf_adaptor->fetch_all_by_Slice_VariationSet($self, $set);  
  }
  return [];
}

=head2 get_all_StructuralVariationFeatures_by_Study

    Arg [1]     : Bio::EnsEMBL:Variation::Study $study
    Arg [2]     : (optional) int $include_evidence
    Arg [3]     : (optional) string $dbtype
                  The dbtype of structural variation to obtain (i.e. can be different from the "variation" type).
                  This assumes that the extra db has been added to the DBAdaptor under this name (using the
                  DBConnection::add_db_adaptor method).    
    Description : returns all structural variation features on this slice associated with a 
                  given study.
                  If $include_evidence is set (i.e. $include_evidence=1), structural variation features from 
                  both structural variation (SV) and their supporting structural variations (SSV) will be 
                  returned. By default, it only returns features from structural variations (SV).
                  This function will only work correctly if the variation database has been
                  attached to the core database. 
    ReturnType  : listref of Bio::EnsEMBL::Variation::StructuralVariationFeature
    Exceptions  : none
    Caller      : contigview, snpview
    Status      : Stable

=cut

sub get_all_StructuralVariationFeatures_by_Study {
  my $self  = shift;
  deprecate('get_all_StructuralVariationFeatures_by_Study is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor::fetch_all_by_Slice_Study instead.');
  my $study = shift;
  my $include_evidence = shift;
  my $dbtype = shift;
  
  if (my $svf_adaptor = $self->_get_StructuralVariationFeatureAdaptor($dbtype)) {
    return $svf_adaptor->fetch_all_by_Slice_Study($self, $study, $include_evidence);  
  }
  return [];
}

=head2 get_all_StructuralVariationFeatures_by_source

    Arg [1]     : string $source
    Arg [2]     : int $include_evidence [optional]
    Arg [3]     : (optional) string $dbtype
                  The dbtype of structural variation to obtain (i.e. can be different from the "variation" type).
                  This assumes that the extra db has been added to the DBAdaptor under this name (using the
                  DBConnection::add_db_adaptor method).    
    Description : returns all structural variation features on this slice associated with a 
                  given source name (e.g. DGVa).
                  If $include_evidence is set (i.e. $include_evidence=1), structural variation features from 
                  both structural variation (SV) and their supporting structural variations (SSV) will be 
                  returned. By default, it only returns features from structural variations (SV).
                  This function will only work correctly if the variation database has been
                  attached to the core database. 
    ReturnType  : listref of Bio::EnsEMBL::Variation::StructuralVariationFeature
    Exceptions  : none
    Caller      : contigview, snpview
    Status      : Stable

=cut

sub get_all_StructuralVariationFeatures_by_source {
  my $self             = shift;
  deprecate('get_all_StructuralVariationFeatures_by_source is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor::fetch_all_by_Slice_Source instead.');
  my $source           = shift;
  my $include_evidence = shift;
  my $dbtype           = shift;
    
  my $constraint = (defined($source)) ? " s.name='$source' " : undef;
  
  if (my $svf_adaptor = $self->_get_StructuralVariationFeatureAdaptor($dbtype)) {
    return $svf_adaptor->fetch_all_by_Slice_constraint($self,$constraint,$include_evidence);  
  }
  return [];
}

=head2 get_all_somatic_StructuralVariationFeatures_by_source

    Arg [1]     : string $source
    Arg [2]     : int $include_evidence [optional]
    Arg [3]     : (optional) string $dbtype
                  The dbtype of structural variation to obtain (i.e. can be different from the "variation" type).
                  This assumes that the extra db has been added to the DBAdaptor under this name (using the
                  DBConnection::add_db_adaptor method).      
    Description : returns all somatic structural variation features on this slice associated with a 
                  given source name (e.g. DGVa).
                  If $include_evidence is set (i.e. $include_evidence=1), structural variation features from 
                  both structural variation (SV) and their supporting structural variations (SSV) will be 
                  returned. By default, it only returns features from structural variations (SV).
                  This function will only work correctly if the variation database has been
                  attached to the core database. 
    ReturnType  : listref of Bio::EnsEMBL::Variation::StructuralVariationFeature
    Exceptions  : none
    Caller      : contigview, snpview
    Status      : Stable

=cut

sub get_all_somatic_StructuralVariationFeatures_by_source {
  my $self             = shift;
  deprecate('get_all_somatic_StructuralVariationFeatures_by_source is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor::fetch_all_somatic_by_Slice_Source instead.');
  my $source           = shift;
  my $include_evidence = shift;
  my $dbtype           = shift;
    
  my $constraint = (defined($source)) ? " s.name='$source' " : undef;
  
  if (my $svf_adaptor = $self->_get_StructuralVariationFeatureAdaptor($dbtype)) {
    return $svf_adaptor->fetch_all_somatic_by_Slice_constraint($self,$constraint,$include_evidence);  
  }
  return [];
}



=head2 get_all_somatic_StructuralVariationFeatures

    Arg [1]     : int $include_evidence [optional]
    Arg [2]     : (optional) string $dbtype
                  The dbtype of structural variation to obtain (i.e. can be different from the "variation" type).
                  This assumes that the extra db has been added to the DBAdaptor under this name (using the
                  DBConnection::add_db_adaptor method).  
    Description : returns all somatic structural variation features on this slice. This function will only work
                  correctly if the variation database has been attached to the core database.
                  By default, it only returns somatic structural variant features which are not labelled 
                  as "CNV_PROBE".
                  If $include_evidence is set (i.e. $include_evidence=1), structural variation features from 
                  both structural variation (SV) and their supporting structural variations (SSV) will be 
                  returned. By default, it only returns features from structural variations (SV).
    ReturnType  : listref of Bio::EnsEMBL::Variation::StructuralVariationFeature
    Exceptions  : none
    Caller      : contigview, snpview, structural_variation_features
    Status      : Stable

=cut

sub get_all_somatic_StructuralVariationFeatures {
  my $self             = shift;
  deprecate('get_all_somatic_StructuralVariationFeatures is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor::fetch_all_somatic_by_Slice_constraint instead.');
  my $include_evidence = shift;
  my $dbtype           = shift;
    
  if (my $svf_adaptor = $self->_get_StructuralVariationFeatureAdaptor($dbtype)) {
    return $svf_adaptor->fetch_all_somatic_by_Slice_constraint($self,undef,$include_evidence);
  }
  return [];
}


=head2 get_all_CopyNumberVariantProbeFeatures

    Arg[1]      : string $source [optional]
    Arg [2]     : (optional) string $dbtype
                  The dbtype of structural variation to obtain (i.e. can be different from the "variation" type).
                  This assumes that the extra db has been added to the DBAdaptor under this name (using the
                  DBConnection::add_db_adaptor method).  
    Description : returns all copy number variant probes on this slice. This function will only work
                  correctly if the variation database has been attached to the core database.
                  If $source is set, only CNV probes with that source name will be returned.
                  If $study is set, only CNV probes of that study will be returned.
    ReturnType  : listref of Bio::EnsEMBL::Variation::StructuralVariationFeature
    Exceptions  : none
    Caller      : contigview, snpview, structural_variation_feature
    Status      : At Risk

=cut

sub get_all_CopyNumberVariantProbeFeatures {
  my $self   = shift;
  deprecate('get_all_CopyNumberVariantProbeFeatures is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor::fetch_all_cnv_probe_by_Slice instead.');
  my $dbtype = shift;
  
  if (my $svf_adaptor = $self->_get_StructuralVariationFeatureAdaptor($dbtype)) {
    return $svf_adaptor->fetch_all_cnv_probe_by_Slice($self, @_);
  }
  return [];
}


=head2 get_all_VariationFeatures_by_Population

    Arg [1]     : Bio::EnsEMBL::Variation::Population
    Arg [2]     : $minimum_frequency [optional]
    Arg [3]     : string $dbtype [optional]
                  The dbtype of variation to obtain (i.e. can be different from the "variation" type).
                  This assumes that the extra db has been added to the DBAdaptor under this name (using the
                  DBConnection::add_db_adaptor method).
    Example     : $pop = $pop_adaptor->fetch_by_dbID(659);
                  @vfs = @{$slice->get_all_VariationFeatures_by_Population($pop,$slice)};
    Description : Retrieves all variation features in a slice which are stored for
                  a specified population. If $minimum_frequency is supplied, only
                  variations with a minor allele frequency (MAF) greater than
                  $minimum_frequency will be returned.
    Returntype  : listref of Bio::EnsEMBL::Variation::VariationFeature
    Exceptions  : throw on incorrect argument
    Caller      : general
    Status      : At Risk

=cut

sub get_all_VariationFeatures_by_Population {
  my ($self, $minimum_frequency, $dbtype) = @_;
  deprecate('get_all_VariationFeatures_by_Population is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor::fetch_all_by_Slice_Population instead.');
  if (my $vf_adaptor = $self->_get_VariationFeatureAdaptor($dbtype)) {
    return $vf_adaptor->fetch_all_by_Slice_Population($self, $minimum_frequency);
  }
  return [];
}


=head2 get_all_PhenotypeFeatures
    Arg [1]     : string $type [optional]
    Arg [2]     : string $dbtype [optional]
                  The dbtype of variation to obtain (i.e. can be different from the "variation" type).
                  This assumes that the extra db has been added to the DBAdaptor under this name (using the
                  DBConnection::add_db_adaptor method).
    Description : Returns all phenotype features on this slice. This function will 
                  only work correctly if the variation database has been attached to the core 
                  database.
                  If $type is specified, only phenotypes annotated on the given
                  object type will be returned
    ReturnType  : listref of Bio::EnsEMBL::Variation::PhenotypeFeature
    Exceptions  : none
    Caller      : web
    Status      : Stable

=cut

sub get_all_PhenotypeFeatures {
  my $self   = shift;
  deprecate('get_all_PhenotypeFeatures is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::DBSQL::PhenotypeFeatureAdaptor::fetch_all_by_Slice instead.');
  my $type   = shift;
  my $dbtype = shift;
  
  
  if(my $pf_adaptor = $self->_get_VariationAdaptor('PhenotypeFeature',$dbtype)) {
    if(defined($type)) {
      return $pf_adaptor->fetch_all_by_Slice_type($self, $type);
    }
    else {
      return $pf_adaptor->fetch_all_by_Slice($self);
    }
  }
  return [];
}

=head2 get_all_IndividualSlice

    Args        : none
    Example     : my $individualSlice = $slice->get_by_Population($population);
    Description : Gets the specific Slice for all the individuls in the population
    ReturnType  : listref of Bio::EnsEMB::IndividualSlice
    Exceptions  : none
    Caller      : general

=cut

sub get_all_IndividualSlice{
    my $self = shift;

    my $individualSliceFactory = Bio::EnsEMBL::IndividualSliceFactory->new(
									   -START   => $self->{'start'},
									   -END     => $self->{'end'},
									   -STRAND  => $self->{'strand'},
									   -ADAPTOR => $self->adaptor(),
									   -SEQ_REGION_NAME => $self->{'seq_region_name'},
									   -SEQ_REGION_LENGTH => $self->{'seq_region_length'},
									   -COORD_SYSTEM    => $self->{'coord_system'},
									   );
    return $individualSliceFactory->get_all_IndividualSlice();
}

=head2 get_by_Individual

    Arg[1]      : Bio::EnsEMBL::Variation::Individual $individual
    Example     : my $individualSlice = $slice->get_by_Individual($individual);
    Description : Gets the specific Slice for the individual
    ReturnType  : Bio::EnsEMB::IndividualSlice
    Exceptions  : none
    Caller      : general

=cut

sub get_by_Individual{
    my $self = shift;
    my $individual = shift;

    return Bio::EnsEMBL::IndividualSlice->new(
					  -START   => $self->{'start'},
					  -END     => $self->{'end'},
					  -STRAND  => $self->{'strand'},
					  -ADAPTOR => $self->adaptor(),
#					  -SEQ     => $self->{'seq'},
					  -SEQ_REGION_NAME => $self->{'seq_region_name'},
					  -SEQ_REGION_LENGTH => $self->{'seq_region_length'},
					  -COORD_SYSTEM    => $self->{'coord_system'},
					  -INDIVIDUAL     => $individual);

}

=head2 get_by_strain

    Arg[1]      : string $strain
    Example     : my $strainSlice = $slice->get_by_strain($strain);
    Description : Gets the specific Slice for the strain
    ReturnType  : Bio::EnsEMB::StrainSlice
    Exceptions  : none
    Caller      : general

=cut

sub get_by_strain {
    my $self = shift;
    my $strain_name = shift;
    deprecate('get_by_strain is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::DBSQL::StrainSliceAdaptor::get_by_strain_Slice instead.');
    return Bio::EnsEMBL::StrainSlice->new(
					  -START   => $self->{'start'},
					  -END     => $self->{'end'},
					  -STRAND  => $self->{'strand'},
					  -ADAPTOR => $self->adaptor(),
					  -SEQ     => $self->{'seq'},
					  -SEQ_REGION_NAME => $self->{'seq_region_name'},
					  -SEQ_REGION_LENGTH => $self->{'seq_region_length'},
					  -COORD_SYSTEM    => $self->{'coord_system'},
					  -STRAIN_NAME     => $strain_name);

}

sub calculate_theta {
    my $self = shift;
    my $strains = shift;
    my $feature = shift; #optional parameter. Name of the feature in the Slice you want to calculate

    if(!$self->adaptor()) {
	warning('Cannot get variation features without attached adaptor');
	return 0;
    }
    my $variation_db = $self->adaptor->db->get_db_adaptor('variation');

    unless($variation_db) {
	warning("Variation database must be attached to core database to " .
		"retrieve variation information" );
	return 0;
    }

    #need to get coverage regions for the slice in the different strains
    my $coverage_adaptor = $variation_db->get_ReadCoverageAdaptor;
    my $strain;
    my $differences = [];
    my $slices = [];
    if ($coverage_adaptor){
	my $num_strains = scalar(@{$strains}) +1;
	if (!defined $feature){
	    #we want to calculate for the whole slice
	    push @{$slices}, $self; #add the slice as the slice to calculate the theta value
	}
	else{
	    #we have features, get the slices for the different features
	    my $features = $self->get_all_Exons();
	    map {push @{$slices},$_->feature_Slice} @{$features}; #add the slices of the features
	}
	my $length_regions = 0;
	my $snps = 0;
	my $theta = 0;
	my $last_position = 0;
	#get all the differences in the slice coordinates
	foreach my $strain_name (@{$strains}){
	    my $strain = $self->get_by_strain($strain_name); #get the strainSlice for the strain

	    my $results = $strain->get_all_differences_Slice;
	    push @{$differences}, @{$results} if (defined $results);
	}
	#when we finish, we have, in max_level, the regions covered by all the sample
	#sort the differences by the genomic position
	my @differences_sorted = sort {$a->start <=> $b->start} @{$differences};
	foreach my $slice (@{$slices}){
	    my $regions_covered = $coverage_adaptor->fetch_all_regions_covered($slice,$strains);
	    if (defined $regions_covered){
		foreach my $range (@{$regions_covered}){
		    $length_regions += ($range->[1] - $range->[0]) + 1; #add the length of the genomic region
		    for (my $i = $last_position;$i<@differences_sorted;$i++){
			if ($differences_sorted[$i]->start >= $range->[0] && $differences_sorted[$i]->end <= $range->[1]){
			    $snps++; #count differences in the region
			}
			elsif ($differences_sorted[$i]->end > $range->[1]){
			    $last_position = $i;
			    last;
			}
		    }
		}
		#when all the ranges have been iterated, calculate rho
		#this is an intermediate variable called a in the formula
		#  a = sum i=2..strains 1/i-1
	    }
	}
	my $a = _calculate_a($num_strains);
	$theta = $snps / ($a * $length_regions);
	return $theta;
    }
    else{
	return 0;
    }
}

sub _calculate_a {
    my $max_level = shift;

    my $a = 0;
    for (my $i = 2; $i <= $max_level+1;$i++){
	$a += 1/($i-1);
    }
    return $a;
}

sub calculate_pi {
    my $self = shift;
    my $strains = shift;
    my $feature = shift;

    if(!$self->adaptor()) {
	warning('Cannot get variation features without attached adaptor');
	return 0;
    }
    my $variation_db = $self->adaptor->db->get_db_adaptor('variation');

    unless($variation_db) {
	warning("Variation database must be attached to core database to " .
		"retrieve variation information" );
	return 0;
    }

    #need to get coverage regions for the slice in the different strains
    my $coverage_adaptor = $variation_db->get_ReadCoverageAdaptor;
    my $differences = [];
    my $slices = [];
    if ($coverage_adaptor){
	my $num_strains = scalar(@{$strains}) +1;
	if (!defined $feature){
	    #we want to calculate for the whole slice
	    push @{$slices}, $self; #add the slice as the slice to calculate the theta value
	}
	else{
	    #we have features, get the slices for the different features
	    my $features = $self->get_all_Exons();
	    map {push @{$slices},$_->feature_Slice} @{$features}; #add the slices of the features
	}
	my @range_differences = ();
	my $pi = 0;
	my $regions = 0;
	my $last_position = 0; #last position visited in the sorted list of differences
	my $triallelic = 0;
	my $is_triallelic = 0;
	foreach my $slice (@{$slices}){
	    foreach my $strain_name (@{$strains}){
		my $strain = $slice->get_by_strain($strain_name); #get the strainSlice for the strain
		my $results = $strain->get_all_differences_Slice;
		push @{$differences}, @{$results} if (defined $results);
	    }
	    my @differences_sorted = sort {$a->start <=> $b->start} @{$differences};

	    my $regions_covered = $coverage_adaptor->fetch_all_regions_covered($slice,$strains);
	    #when we finish, we have, in max_level, the regions covered by all the sample
	    #sort the differences
	    if (defined $regions_covered){
		foreach my $range (@{$regions_covered}){
		    for (my $i = $last_position;$i<@differences_sorted;$i++){
			if ($differences_sorted[$i]->start >= $range->[0] && $differences_sorted[$i]->end <= $range->[1]){
			    #check wether it is the same region or different
			    if (!defined $range_differences[0] || ($differences_sorted[$i]->start == $range_differences[0]->start)){
				if (defined $range_differences[0] && ($differences_sorted[$i]->allele_string ne $range_differences[0]->allele_string)){
				    $is_triallelic = 1;
				}
				push @range_differences, $differences_sorted[$i];
			    }
			    else{
				#new site, calc pi for the previous one
				$pi += 2 * (@range_differences/($num_strains)) * ( 1 - (@range_differences/$num_strains));
				if ($is_triallelic) {
				    $triallelic++;
				    $is_triallelic = 0;
				}
				$regions++;
				@range_differences = ();
				#and start a new range
				push @range_differences, $differences_sorted[$i];
			    }
			}
			elsif ($differences_sorted[$i]->end > $range->[1]){
			    $last_position = $i;
			    last;
			}
		    }
		    #calculate pi for last site, if any
		    if (defined $range_differences[0]){
			$pi += 2 * (@range_differences/$num_strains) * ( 1 - (@range_differences/$num_strains));
			$regions++;
		    }
		}
	    }
	    $pi = $pi / $regions; #calculate average pi
	    print "Regions with variations in region $regions and triallelic $triallelic\n\n";
	}
	return $pi;
    }
    else{
	return 0;
    }

}

=head2 get_all_genotyped_VariationFeatures

    Args       : none
    Function   : returns all variation features on this slice that have been genotyped. This function will only work
                correctly if the variation database has been attached to the core database.
    ReturnType : listref of Bio::EnsEMBL::Variation::VariationFeature
    Exceptions : none
    Caller     : contigview, snpview, ldview
    Status     : At Risk
               : Variation database is under development.

=cut

sub get_all_genotyped_VariationFeatures{
  my $self = shift;
  deprecate('get_all_genotyped_VariationFeatures is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor::fetch_all_genotyped_by_Slice instead.');
  if( my $vf_adaptor = $self->_get_VariationFeatureAdaptor) {
    return $vf_adaptor->fetch_all_genotyped_by_Slice($self);
  } 
  return [];
}


=head2 get_all_Genes

  Arg [1]    : (optional) string $logic_name
               The name of the analysis used to generate the genes to retrieve
  Arg [2]    : (optional) string $dbtype
               The dbtype of genes to obtain.  This assumes that the db has
               been added to the DBAdaptor under this name (using the
               DBConnection::add_db_adaptor method).
  Arg [3]    : (optional) boolean $load_transcripts
               If set to true, transcripts will be loaded immediately rather
               than being lazy-loaded on request.  This will result in a
               significant speed up if the Transcripts and Exons are going to
               be used (but a slow down if they are not).
  Arg [4]    : (optional) string $source
               The source of the genes to retrieve.
  Arg [5]    : (optional) string $biotype
               The biotype of the genes to retrieve.
  Example    : @genes = @{$slice->get_all_Genes};
  Description: Retrieves all genes that overlap this slice, including those on
               the reverse strand.
  Returntype : listref of Bio::EnsEMBL::Genes
  Exceptions : none
  Caller     : none
  Status     : Stable

=cut

sub get_all_Genes{
  my ($self, $logic_name, $dbtype, $load_transcripts, $source, $biotype) = @_;
  if(my $adaptor = $self->_get_CoreAdaptor('Gene', $dbtype)) {
    return $adaptor->fetch_all_by_Slice( $self, $logic_name, $load_transcripts, $source, $biotype);
  }
  return [];
}

=head2 get_all_Genes_by_type

  Arg [1]    : string $type
               The biotype of genes wanted.
  Arg [2]    : (optional) string $logic_name
  Arg [3]    : (optional) boolean $load_transcripts
               If set to true, transcripts will be loaded immediately rather
               than being lazy-loaded on request.  This will result in a
               significant speed up if the Transcripts and Exons are going to
               be used (but a slow down if they are not).
  Example    : @genes = @{$slice->get_all_Genes_by_type('protein_coding',
               'ensembl')};
  Description: Retrieves genes that overlap this slice of biotype $type.
               This is primarily used by the genebuilding code when several
               biotypes of genes are used.

               The logic name is the analysis of the genes that are retrieved.
               If not provided all genes will be retrieved instead. Both
               positive and negative strand Genes will be returned.

  Returntype : listref of Bio::EnsEMBL::Genes
  Exceptions : none
  Caller     : genebuilder, general
  Status     : Stable

=cut

sub get_all_Genes_by_type {
  my ($self, $type, $logic_name, $load_transcripts) = @_;
  return $self->get_all_Genes($logic_name, undef, $load_transcripts, undef, $type);
}


=head2 get_all_Genes_by_source

  Arg [1]    : string source
  Arg [2]    : (optional) boolean $load_transcripts
               If set to true, transcripts will be loaded immediately rather
               than being lazy-loaded on request.  This will result in a
               significant speed up if the Transcripts and Exons are going to
               be used (but a slow down if they are not).
  Example    : @genes = @{$slice->get_all_Genes_by_source('ensembl')};
  Description: Retrieves genes that overlap this slice of source $source.
               Strand of the Slice does not affect the result.
  Returntype : listref of Bio::EnsEMBL::Genes
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Genes_by_source {
  my ($self, $source, $load_transcripts) = @_;
  return  $self->get_all_Genes(undef, undef, $load_transcripts, $source);
}

=head2 get_all_Transcripts

  Arg [1]    : (optional) boolean $load_exons
               If set to true exons will not be lazy-loaded but will instead
               be loaded right away.  This is faster if the exons are
               actually going to be used right away.
  Arg [2]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Arg [3]    : (optional) string $db_type
  Example    : @transcripts = @{$slice->get_all_Transcripts)_};
  Description: Gets all transcripts which overlap this slice.  If you want to
               specify a particular analysis or type, then you are better off
               using get_all_Genes or get_all_Genes_by_type and iterating
               through the transcripts of each gene. Strand of the Slice is
               ignored.
  Returntype : reference to a list of Bio::EnsEMBL::Transcripts
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Transcripts {
  my ($self, $load_exons, $logic_name, $dbtype, $source, $biotype) = @_;
  if(my $adaptor = $self->_get_CoreAdaptor('Transcript', $dbtype)) {
    return $adaptor->fetch_all_by_Slice($self, $load_exons, $logic_name, undef, $source, $biotype);
  }
  return [];
}

=head2 get_all_Transcripts_by_type

  Arg [1]    : string $type
               The biotype of transcripts wanted.
  Arg [2]    : (optional) string $logic_name
  Arg [3]    : (optional) boolean $load_exons
               If set to true exons will not be lazy-loaded but will instead
               be loaded right away.  This is faster if the exons are
               actually going to be used right away.

  Example    : @transcripts = @{$slice->get_all_Transcripts_by_type('protein_coding',
               'ensembl')};
  Description: Retrieves transcripts that overlap this slice of biotype $type.
               This is primarily used by the genebuilding code when several
               biotypes of transcripts are used.

               The logic name is the analysis of the transcripts that are retrieved.
               If not provided all transcripts will be retrieved instead. Both
               positive and negative strand transcripts will be returned.

  Returntype : listref of Bio::EnsEMBL::Transcripts
  Exceptions : none
  Caller     : genebuilder, general
  Status     : Stable

=cut

sub get_all_Transcripts_by_type {
  my ($self, $biotype, $logic_name, $load_exons) = @_;
  return $self->get_all_Transcripts($load_exons, $logic_name, undef, undef, $biotype);
}


=head2 get_all_Transcripts_by_source

  Arg [1]    : string source
  Arg [2]    : (optional) boolean $load_exons
               If set to true exons will not be lazy-loaded but will instead
               be loaded right away.  This is faster if the exons are
               actually going to be used right away.
  Example    : @transcripts = @{$slice->get_all_Transcripts_by_source('ensembl')};
  Description: Retrieves transcripts that overlap this slice of source $source.
               Strand of the Slice does not affect the result.
  Returntype : listref of Bio::EnsEMBL::Transcripts
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Transcripts_by_source {
  my ($self, $source, $load_exons) = @_;
  return $self->get_all_Transcripts($load_exons, undef, undef, $source);

}


=head2 get_all_Exons

  Arg [1]    : none
  Example    : @exons = @{$slice->get_all_Exons};
  Description: Gets all exons which overlap this slice.  Note that these exons
               will not be associated with any transcripts, so this may not
               be terribly useful.
  Returntype : reference to a list of Bio::EnsEMBL::Exons
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Exons {
  my $self = shift;
  if(!$self->adaptor()) {
    warning('Cannot get Exons without attached adaptor');
    return [];
  }
  return $self->adaptor->db->get_ExonAdaptor->fetch_all_by_Slice($self);
}

=head2 get_all_KaryotypeBands

  Arg [1]    : none
  Example    : @kary_bands = @{$slice->get_all_KaryotypeBands};
  Description: Retrieves the karyotype bands which this slice overlaps.
  Returntype : listref oif Bio::EnsEMBL::KaryotypeBands
  Exceptions : none
  Caller     : general, contigview
  Status     : Stable

=cut

sub get_all_KaryotypeBands {
  my ($self) = @_;
  if (my $adaptor = $self->_get_CoreAdaptor('KaryotypeBand')) {
    return $adaptor->fetch_all_by_Slice($self);
  }
  return [];
}

=head2 get_repeatmasked_seq

  Arg [1]    : listref of strings $logic_names (optional)
  Arg [2]    : int $soft_masking_enable (optional)
  Arg [3]    : hash reference $not_default_masking_cases (optional, default is {})
               The values are 0 or 1 for hard and soft masking respectively
               The keys of the hash should be of 2 forms
               "repeat_class_" . $repeat_consensus->repeat_class,
                e.g. "repeat_class_SINE/MIR"
               "repeat_name_" . $repeat_consensus->name
                e.g. "repeat_name_MIR"
               depending on which base you want to apply the not default
               masking either the repeat_class or repeat_name. Both can be
               specified in the same hash at the same time, but in that case,
               repeat_name setting has priority over repeat_class. For example,
               you may have hard masking as default, and you may want soft
               masking of all repeat_class SINE/MIR, but repeat_name AluSp
               (which are also from repeat_class SINE/MIR).
               Your hash will be something like {"repeat_class_SINE/MIR" => 1,
                                                 "repeat_name_AluSp" => 0}
  Example    : $rm_slice = $slice->get_repeatmasked_seq();
               $softrm_slice = $slice->get_repeatmasked_seq(['RepeatMask'],1);
  Description: Returns Bio::EnsEMBL::Slice that can be used to create repeat
               masked sequence instead of the regular sequence.
               Sequence returned by this new slice will have repeat regions
               hardmasked by default (sequence replaced by N) or
               or soft-masked when arg[2] = 1 (sequence in lowercase)
               Will only work with database connection to get repeat features.
  Returntype : Bio::EnsEMBL::RepeatMaskedSlice
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_repeatmasked_seq {
    my ($self,$logic_names,$soft_mask,$not_default_masking_cases) = @_;

    return Bio::EnsEMBL::RepeatMaskedSlice->new
      (-START   => $self->{'start'},
       -END     => $self->{'end'},
       -STRAND  => $self->{'strand'},
       -ADAPTOR => $self->adaptor(),
       -SEQ     => $self->{'seq'},
       -SEQ_REGION_NAME => $self->{'seq_region_name'},
       -SEQ_REGION_LENGTH => $self->{'seq_region_length'},
       -COORD_SYSTEM    => $self->{'coord_system'},
       -REPEAT_MASK     => $logic_names,
       -SOFT_MASK       => $soft_mask,
       -NOT_DEFAULT_MASKING_CASES => $not_default_masking_cases);
}



=head2 _mask_features

  Arg [1]    : reference to a string $dnaref
  Arg [2]    : array_ref $repeats
               reference to a list Bio::EnsEMBL::RepeatFeature
               give the list of coordinates to replace with N or with
               lower case
  Arg [3]    : int $soft_masking_enable (optional)
  Arg [4]    : hash reference $not_default_masking_cases (optional, default is {})
               The values are 0 or 1 for hard and soft masking respectively
               The keys of the hash should be of 2 forms
               "repeat_class_" . $repeat_consensus->repeat_class,
                e.g. "repeat_class_SINE/MIR"
               "repeat_name_" . $repeat_consensus->name
                e.g. "repeat_name_MIR"
               depending on which base you want to apply the not default masking either
               the repeat_class or repeat_name. Both can be specified in the same hash
               at the same time, but in that case, repeat_name setting has priority over
               repeat_class. For example, you may have hard masking as default, and
               you may want soft masking of all repeat_class SINE/MIR,
               but repeat_name AluSp (which are also from repeat_class SINE/MIR).
               Your hash will be something like {"repeat_class_SINE/MIR" => 1,
                                                 "repeat_name_AluSp" => 0}
  Example    : none
  Description: replaces string positions described in the RepeatFeatures
               with Ns (default setting), or with the lower case equivalent
               (soft masking).  The reference to a dna string which is passed
               is changed in place.
  Returntype : none
  Exceptions : none
  Caller     : seq
  Status     : Stable

=cut

sub _mask_features {
  my ($self,$dnaref,$repeats,$soft_mask,$not_default_masking_cases) = @_;

  $soft_mask = 0 unless (defined $soft_mask);
  $not_default_masking_cases = {} unless (defined $not_default_masking_cases);

  # explicit CORE::length call, to avoid any confusion with the Slice
  # length method
  my $dnalen = CORE::length($$dnaref);

 REP:foreach my $old_f (@{$repeats}) {
    my $f = $old_f->transfer( $self );
    my $start  = $f->start;
    my $end    = $f->end;
    my $length = ($end - $start) + 1;

    # check if we get repeat completely outside of expected slice range
    if ($end < 1 || $start > $dnalen) {
      # warning("Unexpected: Repeat completely outside slice coordinates.");
      next REP;
    }

    # repeat partly outside slice range, so correct
    # the repeat start and length to the slice size if needed
    if ($start < 1) {
      $start = 1;
      $length = ($end - $start) + 1;
    }

    # repeat partly outside slice range, so correct
    # the repeat end and length to the slice size if needed
    if ($end > $dnalen) {
      $end = $dnalen;
      $length = ($end - $start) + 1;
    }

    $start--;

    my $padstr;
    # if we decide to define masking on the base of the repeat_type, we'll need
    # to add the following, and the other commented line few lines below.
    my $rc_class;
    my $rc_name;

    if ($f->isa('Bio::EnsEMBL::RepeatFeature')) {
      $rc_class = "repeat_class_" . $f->repeat_consensus->repeat_class;
      $rc_name = "repeat_name_" . $f->repeat_consensus->name;
    }

    my $masking_type;
    $masking_type = $not_default_masking_cases->{$rc_class} if (defined $not_default_masking_cases->{$rc_class});
    $masking_type = $not_default_masking_cases->{$rc_name} if (defined $not_default_masking_cases->{$rc_name});

    $masking_type = $soft_mask unless (defined $masking_type);

    if ($masking_type) {
      $padstr = lc substr ($$dnaref,$start,$length);
    } else {
      $padstr = 'N' x $length;
    }
    substr ($$dnaref,$start,$length) = $padstr;
  }
}


=head2 get_all_SearchFeatures

  Arg [1]    : scalar $ticket_ids
  Example    : $slice->get_all_SearchFeatures('BLA_KpUwwWi5gY');
  Description: Retrieves all search features for stored blast
               results for the ticket that overlap this slice
  Returntype : listref of Bio::EnsEMBL::SeqFeatures
  Exceptions : none
  Caller     : general (webby!)
  Status     : Stable

=cut

sub get_all_SearchFeatures {
  my $self = shift;
  my $ticket = shift;
  local $_;
  unless($ticket) {
    throw("ticket argument is required");
  }

  if(!$self->adaptor()) {
    warning("Cannot get SearchFeatures without an attached adaptor");
    return [];
  }

  my $sfa = $self->adaptor()->db()->get_db_adaptor('blast');

  my $offset = $self->start-1;

  my $features = $sfa ? $sfa->get_all_SearchFeatures($ticket, $self->seq_region_name, $self->start, $self->end) : [];

  foreach( @$features ) {
    $_->start( $_->start - $offset );
    $_->end(   $_->end   - $offset );
  };
  return $features;

}

=head2 get_all_AssemblyExceptionFeatures

  Example    : $slice->get_all_AssemblyExceptionFeatures();
  Description: Retrieves all misc features which overlap this slice. If
               a set code is provided only features which are members of
               the requested set are returned.
  Returntype : listref of Bio::EnsEMBL::AssemblyExceptionFeatures
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_AssemblyExceptionFeatures {
  my ($self) = @_;
  if(my $adaptor = $self->_get_CoreAdaptor('AssemblyExceptionFeature')) {
    return $adaptor->fetch_all_by_Slice($self);
  }
  return [];
}



=head2 get_all_MiscFeatures

  Arg [1]    : string $set (optional)
  Arg [2]    : string $database (optional)
  Example    : $slice->get_all_MiscFeatures('cloneset');
  Description: Retrieves all misc features which overlap this slice. If
               a set code is provided only features which are members of
               the requested set are returned.
  Returntype : listref of Bio::EnsEMBL::MiscFeatures
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_MiscFeatures {
  my ($self, $misc_set, $dbtype) = @_;
  if(my $adaptor = $self->_get_CoreAdaptor('MiscFeature', $dbtype)) {
    if($misc_set) {
      return $adaptor->fetch_all_by_Slice_and_set_code($self,$misc_set);
    }
    return $adaptor->fetch_all_by_Slice($self);
  }
  return [];
}

=head2 get_all_MarkerFeatures

  Arg [1]    : (optional) string logic_name
               The logic name of the marker features to retrieve
  Arg [2]    : (optional) int $priority
               Lower (exclusive) priority bound of the markers to retrieve
  Arg [3]    : (optional) int $map_weight
               Upper (exclusive) priority bound of the markers to retrieve
  Example    : my @markers = @{$slice->get_all_MarkerFeatures(undef,50, 2)};
  Description: Retrieves all markers which lie on this slice fulfilling the
               specified map_weight and priority parameters (if supplied).
  Returntype : reference to a list of Bio::EnsEMBL::MarkerFeatures
  Exceptions : none
  Caller     : contigview, general
  Status     : Stable

=cut

sub get_all_MarkerFeatures {
  my ($self, $logic_name, $priority, $map_weight) = @_;
  if(my $adaptor = $self->_get_CoreAdaptor('MarkerFeature')) {
    return $adaptor->fetch_all_by_Slice_and_priority($self, $priority, $map_weight, $logic_name);
  }
  return [];
}


=head2 get_MarkerFeatures_by_Name

  Arg [1]    : string marker Name
               The name (synonym) of the marker feature(s) to retrieve
  Example    : my @markers = @{$slice->get_MarkerFeatures_by_Name('z1705')};
  Description: Retrieves all markers with this ID
  Returntype : reference to a list of Bio::EnsEMBL::MarkerFeatures
  Exceptions : none
  Caller     : contigview, general
  Status     : Stable

=cut

sub get_MarkerFeatures_by_Name {
  my ($self, $name) = @_;
  if(my $adaptor = $self->_get_CoreAdaptor('MarkerFeature')) {
    return $adaptor->fetch_all_by_Slice_and_MarkerName($self, $name);
  }
  return [];
}


=head2 get_all_compara_DnaAlignFeatures

  Arg [1]    : string $qy_species
               The name of the species to retrieve similarity features from
  Arg [2]    : string $qy_assembly
               The name of the assembly to retrieve similarity features from
  Arg [3]    : string $type
               The type of the alignment to retrieve similarity features from
  Arg [4]    : <optional> compara dbadptor to use.
  Example    : $fs = $slc->get_all_compara_DnaAlignFeatures('Mus musculus',
							    'MGSC3',
							    'WGA');
  Description: Retrieves a list of DNA-DNA Alignments to the species specified
               by the $qy_species argument.
               The compara database must be attached to the core database
               for this call to work correctly.  As well the compara database
               must have the core dbadaptors for both this species, and the
               query species added to function correctly.
  Returntype : reference to a list of Bio::EnsEMBL::DnaDnaAlignFeatures
  Exceptions : warning if compara database is not available
  Caller     : contigview
  Status     : Stable

=cut

sub get_all_compara_DnaAlignFeatures {
  my ($self, $qy_species, $qy_assembly, $alignment_type, $compara_db) = @_;

  if(!$self->adaptor()) {
    warning("Cannot retrieve DnaAlignFeatures without attached adaptor");
    return [];
  }

  unless($qy_species && $alignment_type # && $qy_assembly
  ) {
    throw("Query species and assembly and alignmemt type arguments are required");
  }

  if(!defined($compara_db)){
    $compara_db = Bio::EnsEMBL::Registry->get_DBAdaptor("compara", "compara");
  }
  unless($compara_db) {
    warning("Compara database must be attached to core database or passed ".
	    "as an argument to " .
	    "retrieve compara information");
    return [];
  }

  my $dafa = $compara_db->get_DnaAlignFeatureAdaptor;
  return $dafa->fetch_all_by_Slice($self, $qy_species, $qy_assembly, $alignment_type);
}

=head2 get_all_compara_Syntenies

  Arg [1]    : string $query_species e.g. "Mus_musculus" or "Mus musculus"
  Arg [2]    : string $method_link_type, default is "SYNTENY"
  Arg [3]    : <optional> compara dbadaptor to use.
  Description: gets all the compara syntenyies for a specfic species
  Returns    : arrayref of Bio::EnsEMBL::Compara::SyntenyRegion
  Status     : Stable

=cut

sub get_all_compara_Syntenies {
  my ($self, $qy_species, $method_link_type, $compara_db) = @_;

  if(!$self->adaptor()) {
    warning("Cannot retrieve features without attached adaptor");
    return [];
  }

  unless($qy_species) {
    throw("Query species and assembly arguments are required");
  }

  unless (defined $method_link_type) {
    $method_link_type = "SYNTENY";
  }

  if(!defined($compara_db)){
    $compara_db = Bio::EnsEMBL::Registry->get_DBAdaptor("compara", "compara");
  }
  unless($compara_db) {
    warning("Compara database must be attached to core database or passed ".
	    "as an argument to " .
	    "retrieve compara information");
    return [];
  }
  my $gdba = $compara_db->get_GenomeDBAdaptor();
  my $mlssa = $compara_db->get_MethodLinkSpeciesSetAdaptor();
  my $dfa = $compara_db->get_DnaFragAdaptor();
  my $sra = $compara_db->get_SyntenyRegionAdaptor();

  my $this_gdb = $gdba->fetch_by_core_DBAdaptor($self->adaptor()->db());
  my $query_gdb = $gdba->fetch_by_registry_name($qy_species);
  my $mlss;
  if($this_gdb eq $query_gdb) {
    $mlss = $mlssa->fetch_by_method_link_type_GenomeDBs($method_link_type, [$this_gdb]);
  } else {
    $mlss = $mlssa->fetch_by_method_link_type_GenomeDBs($method_link_type, [$this_gdb, $query_gdb]);
  }

  my $cs = $self->coord_system()->name();
  my $sr = $self->seq_region_name();
  my ($dnafrag) = @{$dfa->fetch_all_by_GenomeDB_region($this_gdb, $cs, $sr)};
  return $sra->fetch_all_by_MethodLinkSpeciesSet_DnaFrag($mlss, $dnafrag, $self->start, $self->end);
}

=head2 get_all_Haplotypes

  Arg [1]    : (optional) boolean $lite_flag
               if true lightweight haplotype objects are used
  Example    : @haplotypes = $slice->get_all_Haplotypes;
  Description: Retrieves all of the haplotypes on this slice.  Only works
               if the haplotype adaptor has been attached to the core adaptor
               via $dba->add_db_adaptor('haplotype', $hdba);
  Returntype : listref of Bio::EnsEMBL::External::Haplotype::Haplotypes
  Exceptions : warning is Haplotype database is not available
  Caller     : contigview, general
  Status     : Stable

=cut

sub get_all_Haplotypes {
  my($self, $lite_flag) = @_;

  if(!$self->adaptor()) {
    warning("Cannot retrieve features without attached adaptor");
    return [];
  }

  my $haplo_db = $self->adaptor->db->get_db_adaptor('haplotype');

  unless($haplo_db) {
    warning("Haplotype database must be attached to core database to " .
		"retrieve haplotype information" );
    return [];
  }

  my $haplo_adaptor = $haplo_db->get_HaplotypeAdaptor;

  my $haplotypes = $haplo_adaptor->fetch_all_by_Slice($self, $lite_flag);

  return $haplotypes;
}


sub get_all_DASFactories {
   my $self = shift;
   return [ $self->adaptor()->db()->_each_DASFeatureFactory ];
}

sub get_all_DASFeatures_dsn {
   my ($self, $source_type, $dsn) = @_;

  if(!$self->adaptor()) {
    warning("Cannot retrieve features without attached adaptor");
    return [];
  }
  my @X = grep { $_->adaptor->dsn eq $dsn } $self->adaptor()->db()->_each_DASFeatureFactory;

  return [ $X[0]->fetch_all_Features( $self, $source_type ) ];
}

=head2 get_all_DAS_Features

  Arg [1]    : none
  Example    : $features = $slice->get_all_DASFeatures;
  Description: Retrieves a hash reference to a hash of DAS feature
               sets, keyed by the DNS, NOTE the values of this hash
               are an anonymous array containing:
                (1) a pointer to an array of features;
                (2) a pointer to the DAS stylesheet
  Returntype : hashref of Bio::SeqFeatures
  Exceptions : ?
  Caller     : webcode
  Status     : Stable

=cut
sub get_all_DAS_Features{
  my ($self) = @_;

  $self->{_das_features} ||= {}; # Cache
  $self->{_das_styles} ||= {}; # Cache
  $self->{_das_segments} ||= {}; # Cache
  my %das_features;
  my %das_styles;
  my %das_segments;
  my $slice = $self;

  foreach my $dasfact( @{$self->get_all_DASFactories} ){
    my $dsn = $dasfact->adaptor->dsn;
    my $name = $dasfact->adaptor->name;
#    my $type = $dasfact->adaptor->type;
    my $url = $dasfact->adaptor->url;

 my ($type) = $dasfact->adaptor->mapping;
 if (ref $type eq 'ARRAY') {
   $type = shift @$type;
 }
 $type ||= $dasfact->adaptor->type;
    # Construct a cache key : SOURCE_URL/TYPE
    # Need the type to handle sources that serve multiple types of features

    my $key = join('/', $name, $type);
    if( $self->{_das_features}->{$key} ){ # Use cached
        $das_features{$name} = $self->{_das_features}->{$key};
        $das_styles{$name} = $self->{_das_styles}->{$key};
        $das_segments{$name} = $self->{_das_segments}->{$key};
    } else { # Get fresh data
        my ($featref, $styleref, $segref) = $dasfact->fetch_all_Features( $slice, $type );
        $self->{_das_features}->{$key} = $featref;
        $self->{_das_styles}->{$key} = $styleref;
        $self->{_das_segments}->{$key} = $segref;
        $das_features{$name} = $featref;
        $das_styles{$name} = $styleref;
        $das_segments{$name} = $segref;
    }
  }

  return (\%das_features, \%das_styles, \%das_segments);
}

sub get_all_DASFeatures{
   my ($self, $source_type) = @_;


  if(!$self->adaptor()) {
    warning("Cannot retrieve features without attached adaptor");
    return [];
  }

  my %genomic_features = map { ( $_->adaptor->dsn => [ $_->fetch_all_Features($self, $source_type) ]  ) } $self->adaptor()->db()->_each_DASFeatureFactory;
  return \%genomic_features;

}

sub old_get_all_DASFeatures{
   my ($self,@args) = @_;

  if(!$self->adaptor()) {
    warning("Cannot retrieve features without attached adaptor");
    return [];
  }

  my %genomic_features =
      map { ( $_->adaptor->dsn => [ $_->fetch_all_by_Slice($self) ]  ) }
         $self->adaptor()->db()->_each_DASFeatureFactory;
  return \%genomic_features;

}


=head2 get_all_ExternalFeatures

  Arg [1]    : (optional) string $track_name
               If specified only features from ExternalFeatureAdaptors with
               the track name $track_name are retrieved.
               If not set, all features from every ExternalFeatureAdaptor are
               retrieved.
  Example    : @x_features = @{$slice->get_all_ExternalFeatures}
  Description: Retrieves features on this slice from external feature adaptors
  Returntype : listref of Bio::SeqFeatureI implementing objects in slice
               coordinates
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_ExternalFeatures {
   my ($self, $track_name) = @_;

   if(!$self->adaptor()) {
     warning("Cannot retrieve features without attached adaptor");
     return [];
   }

   my $features = [];

   my $xfa_hash = $self->adaptor->db->get_ExternalFeatureAdaptors;
   my @xf_adaptors = ();

   if($track_name) {
     #use a specific adaptor
     if(exists $xfa_hash->{$track_name}) {
       push @xf_adaptors, $xfa_hash->{$track_name};
     }
   } else {
     #use all of the adaptors
     push @xf_adaptors, values %$xfa_hash;
   }


   foreach my $xfa (@xf_adaptors) {
     push @$features, @{$xfa->fetch_all_by_Slice($self)};
   }

   return $features;
}


=head2 get_all_DitagFeatures

  Arg [1]    : (optional) string ditag type
  Arg [1]    : (optional) string logic_name
  Example    : @dna_dna_align_feats = @{$slice->get_all_DitagFeatures};
  Description: Retrieves the DitagFeatures of a specific type which overlap
               this slice. If type is not defined, all features are retrieved.
               Strandedness of the Slice is ignored.
  Returntype : listref of Bio::EnsEMBL::DitagFeatures
  Exceptions : warning if slice does not have attached adaptor
  Caller     : general
  Status     : Stable

=cut

sub get_all_DitagFeatures {
  my ($self, $type, $logic_name) = @_;
  if(my $adaptor = $self->_get_CoreAdaptor('DitagFeature')) {
    return $adaptor->fetch_all_by_Slice($self, $type, $logic_name);
  }
  return [];
}

# GENERIC FEATURES (See DBAdaptor.pm)

=head2 get_generic_features

  Arg [1]    : (optional) List of names of generic feature types to return.
               If no feature names are given, all generic features are
               returned.
  Example    : my %features = %{$slice->get_generic_features()};
  Description: Gets generic features via the generic feature adaptors that
               have been added via DBAdaptor->add_GenricFeatureAdaptor (if
               any)
  Returntype : Hash of named features.
  Exceptions : none
  Caller     : none
  Status     : Stable

=cut

sub get_generic_features {

  my ($self, @names) = @_;

  if(!$self->adaptor()) {
    warning('Cannot retrieve features without attached adaptor');
    return [];
  }

  my $db = $self->adaptor()->db();

  my %features = ();   # this will hold the results

  # get the adaptors for each feature
  my %adaptors = %{$db->get_GenericFeatureAdaptors(@names)};

  foreach my $adaptor_name (keys(%adaptors)) {

    my $adaptor_obj = $adaptors{$adaptor_name};
    # get the features and add them to the hash
    my $features_ref = $adaptor_obj->fetch_all_by_Slice($self);
    # add each feature to the hash to be returned
    foreach my $feature (@$features_ref) {
      $features{$adaptor_name} = $feature;
    }
  }

  return \%features;

}

=head2 project_to_slice

  Arg [1]    : Slice to project to.
  Example    : my $chr_projection = $clone_slice->project_to_slice($chrom_slice);
                foreach my $segment ( @$chr_projection ){
                  $chr_slice = $segment->to_Slice();
                  print $clone_slice->seq_region_name(). ':'. $segment->from_start(). '-'.
                        $segment->from_end(). ' -> '.$chr_slice->seq_region_name(). ':'. $chr_slice->start().
	                '-'.$chr_slice->end().
                         $chr_slice->strand(). " length: ".($chr_slice->end()-$chr_slice->start()+1). "\n";
                }
  Description: Projection of slice to another specific slice. Needed for where we have multiple mappings
               and we want to state which one to project to.
  Returntype : list reference of Bio::EnsEMBL::ProjectionSegment objects which
               can also be used as [$start,$end,$slice] triplets.
  Exceptions : none
  Caller     : none
  Status     : At Risk

=cut

sub project_to_slice {
  my $self = shift;
  my $to_slice = shift;

  throw('Slice argument is required') if(!$to_slice);

  my $slice_adaptor = $self->adaptor();

  if(!$slice_adaptor) {
    warning("Cannot project without attached adaptor.");
    return [];
  }


  my $mapper_aptr = $slice_adaptor->db->get_AssemblyMapperAdaptor();

  my $cs = $to_slice->coord_system();
  my $slice_cs = $self->coord_system();
  my $to_slice_id = $to_slice->get_seq_region_id;

  my @projection;
  my $current_start = 1;

  # decompose this slice into its symlinked components.
  # this allows us to handle haplotypes and PARs
  my $normal_slice_proj =
    $slice_adaptor->fetch_normalized_slice_projection($self);
  foreach my $segment (@$normal_slice_proj) {
    my $normal_slice = $segment->[2];

    $slice_cs = $normal_slice->coord_system();

    my $asma = $self->adaptor->db->get_AssemblyMapperAdaptor();
    my $asm_mapper = $asma->fetch_by_CoordSystems($slice_cs, $cs);

    # perform the mapping between this slice and the requested system
    my @coords;

    if( defined $asm_mapper ) {
     @coords = $asm_mapper->map($normal_slice->seq_region_name(),
				 $normal_slice->start(),
				 $normal_slice->end(),
				 $normal_slice->strand(),
				 $slice_cs, undef, $to_slice);
    } else {
      $coords[0] = Bio::EnsEMBL::Mapper::Gap->new( $normal_slice->start(),
						   $normal_slice->end());
    }

    my $last_rank =0;
    #construct a projection from the mapping results and return it
    foreach my $coord (@coords) {
      my $coord_start  = $coord->start();
      my $coord_end    = $coord->end();
      my $length       = $coord_end - $coord_start + 1;


      if( $last_rank != $coord->rank){
	$current_start = 1;
      }
      $last_rank = $coord->rank;

      #skip gaps
      if($coord->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
	if($coord->id != $to_slice_id){ # for multiple mappings only get the correct one
	  $current_start += $length;
	  next;
	}
        my $coord_cs     = $coord->coord_system();

        # If the normalised projection just ended up mapping to the
        # same coordinate system we were already in then we should just
        # return the original region.  This can happen for example, if we
        # were on a PAR region on Y which refered to X and a projection to
        # 'toplevel' was requested.
#        if($coord_cs->equals($slice_cs)) {
#          # trim off regions which are not defined
#          return $self->_constrain_to_region();
#        }

        #create slices for the mapped-to coord system
        my $slice = $slice_adaptor->fetch_by_seq_region_id(
                                                    $coord->id(),
                                                    $coord_start,
                                                    $coord_end,
                                                    $coord->strand());

	my $current_end = $current_start + $length - 1;

        push @projection, bless([$current_start, $current_end, $slice],
                                "Bio::EnsEMBL::ProjectionSegment");
      }

      $current_start += $length;
    }
  }


  # delete the cache as we may want to map to different set next time and old
  # results will be cached.

  $mapper_aptr->delete_cache;

  return \@projection;
}


=head2 get_all_synonyms

  Args [1]   : String external_db_name The name of the database to retrieve 
               the synonym for
  Args [2]   : (optional) Integer external_db_version Optionally restrict
               results from external_db_name to a specific version of
               the the specified external database
  Example    : my @alternative_names = @{$slice->get_all_synonyms()};
               @alternative_names = @{$slice->get_all_synonyms('EMBL')};
  Description: get a list of alternative names for this slice
  Returntype : reference to list of SeqRegionSynonym objects.
  Exception  : none
  Caller     : general
  Status     : At Risk

=cut

sub get_all_synonyms {
  my ($self, $external_db_name, $external_db_version) = @_;

  if ( !defined( $self->{'synonym'} ) ) {
    my $adap = $self->adaptor->db->get_SeqRegionSynonymAdaptor();
    $self->{'synonym'} =
      $adap->get_synonyms( $self->get_seq_region_id($self) );
  }
  
  if(! $external_db_name) {
    return $self->{'synonym'};
  }
  my @args =  ($external_db_version) ? 
              ($external_db_name, $external_db_version) : 
              ($external_db_name, undef, 1);
  my $external_db_id = $self->adaptor->db()->get_DBEntryAdaptor()->get_external_db_id(@args);
  if(!$external_db_id) {
    my $extra = ($external_db_version) ? "and version $external_db_version " : q{};
    throw "The external database $external_db_name ${extra}did not result in a valid identifier";
  }

  return [ grep { $_->external_db_id() == $external_db_id } @{$self->{synonym}} ];
}

=head2 add_synonym

  Args[0]    : synonym.
  Example    : $slice->add_synonym("alt_name");
  Description: add an alternative name for this slice
  Returntype : none
  Exception  : none
  Caller     : general
  Status     : At Risk

=cut

sub add_synonym{
  my $self = shift;
  my $syn = shift;
  my $external_db_id = shift;
  
  my $adap = $self->adaptor->db->get_SeqRegionSynonymAdaptor();
  if ( !defined( $self->{'synonym'} ) ) {
    $self->{'synonym'} = $self->get_all_synonyms();
  }
  my $new_syn = Bio::EnsEMBL::SeqRegionSynonym->new( #-adaptor => $adap,
                                                     -synonym => $syn,
                                                     -external_db_id => $external_db_id, 
                                                     -seq_region_id => $self->get_seq_region_id($self));

  push (@{$self->{'synonym'}}, $new_syn);

  return;
}

=head2 summary_as_hash

  Example       : $slice_summary = $slice->summary_as_hash();
  Description   : Retrieves a textual summary of this slice.
  Returns       : hashref of descriptive strings
=cut

sub summary_as_hash {
  my $self = shift;
  my %summary;
  my @aliases = map { $_->name } @{$self->slice->get_all_synonyms()};

  $summary{'seq_region_name'} = $self->seq_region_name;
  $summary{'id'} = $self->seq_region_name;
  $summary{'start'} = $self->start;
  $summary{'end'} = $self->end;
  $summary{'strand'} = 0;
  $summary{'source'} = $self->source || $self->coord_system->version;
  $summary{'Alias'} = \@aliases if scalar(@aliases);
  $summary{'Is_circular'} = $self->is_circular ? "true" : undef;
  return \%summary;
}


sub slice {
  my $self = shift;
  return $self;
}

#
# Bioperl Bio::PrimarySeqI methods:
#

=head2 id

  Description: Included for Bio::PrimarySeqI interface compliance (0.7)

=cut

sub id { name(@_); }


=head2 display_id

  Description: Included for Bio::PrimarySeqI interface compliance (1.2)

=cut

sub display_id { name(@_); }


=head2 primary_id

  Description: Included for Bio::PrimarySeqI interface compliance (1.2)

=cut

sub primary_id { name(@_); }


=head2 desc

Description: Included for Bio::PrimarySeqI interface compliance (1.2)

=cut

sub desc{ return $_[0]->coord_system->name().' '.$_[0]->seq_region_name(); }


=head2 moltype

Description: Included for Bio::PrimarySeqI interface compliance (0.7)

=cut

sub moltype { return 'dna'; }

=head2 alphabet

  Description: Included for Bio::PrimarySeqI interface compliance (1.2)

=cut

sub alphabet { return 'dna'; }


=head2 accession_number

  Description: Included for Bio::PrimarySeqI interface compliance (1.2)

=cut

sub accession_number { name(@_); }

1;
