=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::CircularSlice - Arbitary Slice of a genome

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

A Slice object represents a region of a genome.  It can be used to
retrieve sequence or features from an area of interest.

=head1 METHODS

=cut

package Bio::EnsEMBL::CircularSlice;
use vars qw(@ISA);
use strict;

use Bio::PrimarySeqI;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception
  qw(throw deprecate warning stack_trace_dump);
use Bio::EnsEMBL::RepeatMaskedSlice;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Utils::Scalar qw( assert_ref );
use Bio::EnsEMBL::ProjectionSegment;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::MergedAdaptor;

use Bio::EnsEMBL::StrainSlice;
#use Bio::EnsEMBL::IndividualSlice;
#use Bio::EnsEMBL::IndividualSliceFactory;
use Bio::EnsEMBL::Mapper::RangeRegistry;
use Bio::EnsEMBL::Slice;
use Data::Dumper;
use Scalar::Util qw(weaken isweak);

my $reg = "Bio::EnsEMBL::Registry";

@ISA = qw(Bio::EnsEMBL::Slice);

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
  Example    :
  
    $slice =
      Bio::EnsEMBL::CircularSlice->new( -coord_system      => $cs,
                                        -start             => 1,
                                        -end               => 10000,
                                        -strand            => 1,
                                        -seq_region_name   => 'X',
                                        -seq_region_length => 12e6,
                                        -adaptor => $slice_adaptor );

  Description: Creates a new slice object.  A slice represents a
               region of sequence in a particular coordinate system.
               Slices can be used to retrieve sequence and features
               from an area of interest in a genome.

               Coordinates start at 1 and are inclusive.  Negative
               coordinates or coordinates exceeding the length of
               the seq_region are permitted.  Start must be less
               than or equal. to end regardless of the strand.

               Slice objects are immutable. Once instantiated their
               attributes (with the exception of the adaptor) may
               not be altered.  To change the attributes a new slice
               must be created.

  Returntype : Bio::EnsEMBL::CircularSlice
  Exceptions : throws if start, end, coordsystem or seq_region_name not
               specified or not of the correct type
  Caller     : general, Bio::EnsEMBL::SliceAdaptor
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  #new can be called as a class or object method
  my $class = ref($caller) || $caller;

  my ( $seq, $coord_system, $seq_region_name, $seq_region_length,
       $start, $end, $strand, $adaptor, $empty )
    = rearrange( [
      qw(SEQ COORD_SYSTEM SEQ_REGION_NAME SEQ_REGION_LENGTH
        START END STRAND ADAPTOR EMPTY) ],
    @_ );

  #empty is only for backwards compatibility
  if ($empty) {
    deprecate(   "Creation of empty slices is no longer needed "
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

  if ( !defined($seq_region_length) ) { $seq_region_length = $end }

  if ( $seq_region_length <= 0 ) {
    throw('SEQ_REGION_LENGTH must be > 0');
  }

  if ( defined($coord_system) ) {
    assert_ref( $coord_system, 'Bio::EnsEMBL::CoordSystem' );

    if ( $coord_system->is_top_level() ) {
      throw('Cannot create circular slice on toplevel CoordSystem.');
    }
  } else {
    warning("CircularSlice without coordinate system");
  }

  $strand ||= 1;

  if ( $strand != 1 && $strand != -1 ) {
    throw('STRAND argument must be -1 or 1');
  }

  if ( defined($adaptor) ) {
    assert_ref( $adaptor, 'Bio::EnsEMBL::DBSQL::SliceAdaptor' );
  }

  my $seq1 = { 'coord_system'      => $coord_system,
               'seq'               => $seq,
               'seq_region_name'   => $seq_region_name,
               'seq_region_length' => $seq_region_length,
               'start'             => int($start),
               'end'               => int($end),
               'strand'            => $strand };

  bless $seq1, $class;
  $seq1->adaptor($adaptor);
  return $seq1;
} ## end sub new


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

  my ( $s, $e, $length ) =
    ( $self->{'start'}, $self->{'end'}, $self->{'seq_region_length'} );

  if ( $s < $e ) {
    return ( $s + $e )/2;
  }

  my $r1 = $length - $s;
  my $r2 = $e;
  my $r  = ( $r1 + $r2 )/2;
  my $m  = $s + $r;

  if ( $m > $length ) {
    $m = $m - $length;
  }

  return $m;
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

  if ( $self->{'start'} < $self->{'end'} ) {
    return $self->{'end'} - $self->{'start'} + 1;
  }

  my $r1 = $self->{'seq_region_length'} - $self->{'start'};
  my $r2 = $self->{'end'};
  my $ln = $r1 + $r2 + 1;

  return $ln;
}

sub _split {
  my $self = shift;

  my $sl1 =
    Bio::EnsEMBL::CircularSlice->new(
                     -COORD_SYSTEM      => $self->{'coord_system'},
                     -SEQ_REGION_NAME   => $self->{'seq_region_name'},
                     -SEQ_REGION_LENGTH => $self->{'seq_region_length'},
                     -START             => $self->{'start'},
                     -END               => $self->{'seq_region_length'},
                     -STRAND            => $self->{'strand'},
                     -ADAPTOR           => $self->adaptor() );

  my $sl2 =
    Bio::EnsEMBL::CircularSlice->new(
                     -COORD_SYSTEM      => $self->{'coord_system'},
                     -SEQ_REGION_NAME   => $self->{'seq_region_name'},
                     -SEQ_REGION_LENGTH => $self->{'seq_region_length'},
                     -START             => 1,
                     -END               => $self->{'end'},
                     -STRAND            => $self->{'strand'},
                     -ADAPTOR           => $self->adaptor() );

  return ($sl1, $sl2);
}

=head2 seq

  Arg [1]    : none
  Example    : print "SEQUENCE = ", $slice->seq();
  Description: Returns the sequence of the region represented by this
               slice formatted as a string.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub seq {
  my $self = shift;

  # special case for in-between (insert) coordinates
  return '' if ( $self->start() == $self->end() + 1 );
  return $self->{'seq'} if ( $self->{'seq'} );

  if ( $self->adaptor() ) {

    my $seqAdaptor = $self->adaptor()->db()->get_SequenceAdaptor();
    if ( $self->{'start'} > $self->{'end'} ) {
      my $length = $self->{'seq_region_length'};

      my ($sl1, $sl2) = $self->_split;

      my $seq1 = ${
        $seqAdaptor->fetch_by_Slice_start_end_strand( $sl1, 1, undef,
                                                      1 ) };
      my $seq2 = ${
        $seqAdaptor->fetch_by_Slice_start_end_strand( $sl2, 1, undef,
                                                      1 ) };
      return $seq1 . $seq2;

    } else {
      my $seq1 = ${
        $seqAdaptor->fetch_by_Slice_start_end_strand( $self, 1, undef,
                                                      1 ) };
      return $seq1;
    }
  } ## end if ( $self->adaptor() )

  # no attached sequence, and no db, so just return Ns
  return 'N' x $self->length();
} ## end sub seq

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

  # handle 'between' case for insertions
  return '' if ( $start == $end + 1 );

  $strand = 1 unless ( defined $strand );

  if ( $strand != -1 && $strand != 1 ) {
    throw("Invalid strand [$strand] in call to Slice::subseq.");
  }
  my $subseq;
  my $length = $self->{'seq_region_length'};

  if ( $self->adaptor ) {

    my $seqAdaptor = $self->adaptor->db->get_SequenceAdaptor();
    if ( $end < $start ) {
      my $subseq1 = ${
        $seqAdaptor->fetch_by_Slice_start_end_strand( $self, $start,
                                                      $length, $strand )
        };
      my $subseq2 = ${
        $seqAdaptor->fetch_by_Slice_start_end_strand( $self, 1, $end,
                                                      $strand ) };
      $subseq = $subseq1 . $subseq2;

    } else {
      $subseq = ${
        $seqAdaptor->fetch_by_Slice_start_end_strand( $self, $start,
                                                      $end, $strand ) };
    }
  } else {
    ## check for gap at the beginning and pad it with Ns
    if ( $start < 1 ) {
      $subseq = "N" x ( 1 - $start );
      $start = 1;
    }
    $subseq .= substr( $self->seq(), $start - 1, $end - $start + 1 );
    ## check for gap at the end and pad it with Ns
    if ( $end > $self->length() ) {
      $subseq .= "N" x ( $end - $self->length() );
    }
    reverse_comp( \$subseq ) if ( $strand == -1 );
  }
  return $subseq;
} ## end sub subseq

=head2 project

  Arg [1]    : string $name
               The name of the coordinate system to project this slice onto
  Arg [2]    : string $version
               The version of the coordinate system (such as 'NCBI34') to
               project this slice onto
  Example    :
    my $clone_projection = $slice->project('clone');

    foreach my $seg (@$clone_projection) {
      my $clone = $segment->to_Slice();
      print $slice->seq_region_name(), ':', $seg->from_start(), '-',
        $seg->from_end(), ' -> ',
        $clone->seq_region_name(), ':', $clone->start(), '-',
        $clone->end(),
        $clone->strand(), "\n";
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
  my $self       = shift;
  my $cs_name    = shift;
  my $cs_version = shift;

  throw('Coord_system name argument is required') if ( !$cs_name );

  my $slice_adaptor = $self->adaptor();

  if ( !$slice_adaptor ) {
    warning("Cannot project without attached adaptor.");
    return [];
  }

  if ( !$self->coord_system() ) {
    warning("Cannot project without attached coord system.");
    return [];
  }

  my $db  = $slice_adaptor->db();
  my $csa = $db->get_CoordSystemAdaptor();
  my $cs  = $csa->fetch_by_name( $cs_name, $cs_version );

  my ($sl01, $sl02) = $self->_split;
  my @projection;
  my $current_start = 1;

  foreach my $sl2 ( $sl01, $sl02 ) {
    my $slice_cs = $sl2->coord_system();

    if ( !$cs ) {
      throw(   "Cannot project to unknown coordinate system "
             . "[$cs_name $cs_version]" );
    }

# no mapping is needed if the requested coord system is the one we are in
# but we do need to check if some of the slice is outside of defined regions
    if ( $slice_cs->equals($cs) ) {
      return $self->_constrain_to_region();
    }

    # decompose this slice into its symlinked components.
    # this allows us to handle haplotypes and PARs
    my $normal_slice_proj =
      $slice_adaptor->fetch_normalized_slice_projection($sl2);

    foreach my $segment (@$normal_slice_proj) {
      my $normal_slice = $segment->[2];

      $slice_cs = $normal_slice->coord_system();

      my $asma = $db->get_AssemblyMapperAdaptor();
      my $asm_mapper = $asma->fetch_by_CoordSystems( $slice_cs, $cs );

      # perform the mapping between this slice and the requested system
      my @coords;

      if ( defined $asm_mapper ) {
        @coords = $asm_mapper->map( $normal_slice->seq_region_name(),
                                    $normal_slice->start(),
                                    $normal_slice->end(),
                                    $normal_slice->strand(),
                                    $slice_cs );

      } else {
        $coords[0] =
          Bio::EnsEMBL::Mapper::Gap->new( $normal_slice->start(),
                                          $normal_slice->end() );
      }

      #construct a projection from the mapping results and return it
      foreach my $coord (@coords) {
        my $coord_start = $coord->start();
        my $coord_end   = $coord->end();
        my $length      = $coord_end - $coord_start + 1;

        #skip gaps
        if ( $coord->isa('Bio::EnsEMBL::Mapper::Coordinate') ) {
          my $coord_cs = $coord->coord_system();

      # If the normalised projection just ended up mapping to the
      # same coordinate system we were already in then we should just
      # return the original region.  This can happen for example, if we
      # were on a PAR region on Y which refered to X and a projection to
      # 'toplevel' was requested.

          if ( $coord_cs->equals($slice_cs) ) {
            # trim off regions which are not defined
            return $self->_constrain_to_region();
          }
          #create slices for the mapped-to coord system

          my $slice =
            $slice_adaptor->fetch_by_seq_region_id(
                                          $coord->id(), $coord_start,
                                          $coord_end,   $coord->strand()
            );

          my $current_end = $current_start + $length - 1;
          push @projection,
            bless( [ $current_start, $current_end, $slice ],
                   "Bio::EnsEMBL::ProjectionSegment" );
        } ## end if ( $coord->isa('Bio::EnsEMBL::Mapper::Coordinate'...))

        $current_start += $length;
      } ## end foreach my $coord (@coords)
    } ## end foreach my $segment (@$normal_slice_proj)
  }    #foreach

  return \@projection;
} ## end sub project

sub project_org {
  my $self       = shift;
  my $cs_name    = shift;
  my $cs_version = shift;

  throw('Coord_system name argument is required') if ( !$cs_name );

  my $slice_adaptor = $self->adaptor();

  if ( !$slice_adaptor ) {
    warning("Cannot project without attached adaptor.");
    return [];
  }

  if ( !$self->coord_system() ) {
    warning("Cannot project without attached coord system.");
    return [];
  }

  my $db       = $slice_adaptor->db();
  my $csa      = $db->get_CoordSystemAdaptor();
  my $cs       = $csa->fetch_by_name( $cs_name, $cs_version );
  my $slice_cs = $self->coord_system();

  if ( !$cs ) {
    throw(   "Cannot project to unknown coordinate system "
           . "[$cs_name $cs_version]" );
  }

  # No mapping is needed if the requested coord system is the one we
  # are in.  But we do need to check if some of the slice is outside of
  # defined regions.
  if ( $slice_cs->equals($cs) ) {
    return $self->_constrain_to_region();
  }

  my @projection;
  my $current_start = 1;

  # Decompose this slice into its symlinked components.  This allows us
  # to handle haplotypes and PARs.
  my $normal_slice_proj =
    $slice_adaptor->fetch_normalized_slice_projection($self);
  foreach my $segment (@$normal_slice_proj) {
    my $normal_slice = $segment->[2];

    $slice_cs = $normal_slice->coord_system();

    my $asma = $db->get_AssemblyMapperAdaptor();
    my $asm_mapper = $asma->fetch_by_CoordSystems( $slice_cs, $cs );

    # perform the mapping between this slice and the requested system
    my @coords;

    if ( defined $asm_mapper ) {
      @coords = $asm_mapper->map( $normal_slice->seq_region_name(),
                                  $normal_slice->start(),
                                  $normal_slice->end(),
                                  $normal_slice->strand(),
                                  $slice_cs );

    } else {
      $coords[0] =
        Bio::EnsEMBL::Mapper::Gap->new( $normal_slice->start(),
                                        $normal_slice->end() );
    }

    #construct a projection from the mapping results and return it
    foreach my $coord (@coords) {
      my $coord_start = $coord->start();
      my $coord_end   = $coord->end();
      my $length      = $coord_end - $coord_start + 1;

      #skip gaps
      if ( $coord->isa('Bio::EnsEMBL::Mapper::Coordinate') ) {
        my $coord_cs = $coord->coord_system();

        # If the normalised projection just ended up mapping to the
        # same coordinate system we were already in then we should just
        # return the original region.  This can happen for example,
        # if we were on a PAR region on Y which refered to X and a
        # projection to 'toplevel' was requested.

        if ( $coord_cs->equals($slice_cs) ) {
          # trim off regions which are not defined
          return $self->_constrain_to_region();
        }
        #create slices for the mapped-to coord system

        my $slice =
          $slice_adaptor->fetch_by_seq_region_id( $coord->id(),
                           $coord_start, $coord_end, $coord->strand() );

        my $current_end = $current_start + $length - 1;

        push @projection,
          bless( [ $current_start, $current_end, $slice ],
                 "Bio::EnsEMBL::ProjectionSegment" );
      } ## end if ( $coord->isa('Bio::EnsEMBL::Mapper::Coordinate'...))

      $current_start += $length;
    } ## end foreach my $coord (@coords)
  } ## end foreach my $segment (@$normal_slice_proj)

  return \@projection;
} ## end sub project_org

sub _constrain_to_region {
  my $self = shift;

  my $entire_len = $self->seq_region_length();

  # If the slice has negative coordinates or coordinates exceeding the
  # exceeding length of the sequence region we want to shrink the slice
  # to the defined region.

  if ( $self->{'start'} > $entire_len || $self->{'end'} < 1 ) {
    #none of this slice is in a defined region
    return [];
  }

  my $right_contract = 0;
  my $left_contract  = 0;
  if ( $self->{'end'} > $entire_len ) {
    $right_contract = $entire_len - $self->{'end'};
  }
  if ( $self->{'start'} < 1 ) {
    $left_contract = $self->{'start'} - 1;
  }

  my $new_slice;
  if ( $left_contract || $right_contract ) {
    $new_slice = $self->expand( $left_contract, $right_contract );
  } else {
    $new_slice = $self;
  }

  return [ bless [ 1 - $left_contract,
                   $self->length() + $right_contract,
                   $new_slice ],
           "Bio::EnsEMBL::ProjectionSegment" ];
} ## end sub _constrain_to_region

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

  my $new_start;
  my $new_end;
  my $sshift = $five_prime_shift;
  my $eshift = $three_prime_shift;

  if ( $self->{'strand'} != 1 ) {
    $eshift = $five_prime_shift;
    $sshift = $three_prime_shift;
  }

  $new_start = $self->{'start'} - $sshift;
  $new_end   = $self->{'end'} + $eshift;

#  if($new_start > $new_end) {
#      if ($force_expand) {    # Apply max possible shift, if force_expand is set
#	  if ($sshift < 0) {     # if we are contracting the slice from the start - move the start just before the end
#	      $new_start = $new_end - 1;
#	      $sshift = $self->{start} - $new_start;
#	  }

#	  if($new_start > $new_end) { # if the slice still has a negative length - try to move the end
#	      if ($eshift < 0) {
#		  $new_end = $new_start + 1;
#		  $eshift = $new_end - $self->{end};
#	      }
#	  }
# return the values by which the primes were actually shifted
#	  $$tpref = $self->{strand} == 1 ? $eshift : $sshift;
#	  $$fpref = $self->{strand} == 1 ? $sshift : $eshift;
#      }
#      if($new_start > $new_end) {
#	  throw('Slice start cannot be greater than slice end');
#      }
#  }

  #fastest way to copy a slice is to do a shallow hash copy
  my %new_slice = %$self;
  $new_slice{'start'} = int($new_start);
  $new_slice{'end'}   = int($new_end);

  return bless \%new_slice, ref($self);
} ## end sub expand



=head2 get_all_VariationFeatures

    Args       : $filter [optional]
    Description:returns all variation features on this slice. This function will only work 
                correctly if the variation database has been attached to the core database.
		        If $filter is "genotyped" return genotyped Snps only... (nice likkle hack);
    ReturnType : listref of Bio::EnsEMBL::Variation::VariationFeature
    Exceptions : none
    Caller     : contigview, snpview
    Status     : At Risk
               : Variation database is under development.

=cut

sub get_all_VariationFeatures {
  my $self   = shift;
  my $filter = shift;

  $filter ||= '';
  if ( !$self->adaptor() ) {
    warning('Cannot get variation features without attached adaptor');
    return [];
  }

  my $vf_adaptor =
    Bio::EnsEMBL::DBSQL::MergedAdaptor->new(
                            -species => $self->adaptor()->db()->species,
                            -type    => "VariationFeature" );
  if ($vf_adaptor) {
    if ( $filter eq 'genotyped' ) {
      return $vf_adaptor->fetch_all_genotyped_by_Slice($self);
    } else {
      return $vf_adaptor->fetch_all_by_Slice($self);
    }
  } else {
    warning( "Variation database must be attached to core database to "
             . "retrieve variation information" );
    return [];
  }
}




=head2 get_all_genotyped_VariationFeatures

    Args       : none
    Description: returns all variation features on this slice that have been genotyped.
                 This function will only work correctly if the variation database has 
                 been attached to the core database.
    ReturnType : listref of Bio::EnsEMBL::Variation::VariationFeature
    Exceptions : none
    Caller     : contigview, snpview, ldview
    Status     : At Risk
               : Variation database is under development.

=cut

sub get_all_genotyped_VariationFeatures {
  my $self = shift;
  my $vfa;
  if ( !$self->adaptor() ) {
    warning('Cannot get variation features without attached adaptor');
    return [];
  }

  my $vf_adaptor =
    Bio::EnsEMBL::DBSQL::MergedAdaptor->new(
                            -species => $self->adaptor()->db()->species,
                            -type    => "VariationFeature" );

  if ($vf_adaptor) {
    return $vf_adaptor->fetch_all_genotyped_by_Slice($self);
  } else {
    warning( "Variation database must be attached to core database to "
             . "retrieve variation information" );
    return [];
  }
}


=head2 get_all_QtlFeatures

  Args       : none
  Description: returns overlapping QtlFeatures
  Returntype : listref Bio::EnsEMBL::Map::QtlFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_QtlFeatures {
  my $self = shift;

  if ( !$self->adaptor() ) {
    warning('Cannot get QtlFeatures without attached adaptor');
    return [];
  }

  my $qfAdaptor;
  if ( $self->adaptor() ) {
    $qfAdaptor = $self->adaptor()->db()->get_QtlFeatureAdaptor();
  } else {
    return [];
  }

  ## circular BOF
  my ($sl1, $sl2) = $self->_split;
  my ( @arr, @arr1, @arr2 );
  @arr1 = @{ $qfAdaptor->fetch_all_by_Slice_constraint($sl1) };
  @arr2 = @{ $qfAdaptor->fetch_all_by_Slice_constraint($sl2) };
  push @arr, @arr1, @arr2;
  return \@arr;
  ## circular EOF

} ## end sub get_all_QtlFeatures

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

  if ( !$self->adaptor() ) {
    warning('Cannot get KaryotypeBands without attached adaptor');
    return [];
  }

  my $kadp = $self->adaptor->db->get_KaryotypeBandAdaptor();

  ## circular BOF
  my ($sl1, $sl2) = $self->_split;
  my ( @arr, @arr1, @arr2 );
  @arr1 = @{ $kadp->fetch_all_by_Slice($sl1) };
  @arr2 = @{ $kadp->fetch_all_by_Slice($sl2) };
  push @arr, @arr1, @arr2;
  return \@arr;
  ## circular EOF

} ## end sub get_all_KaryotypeBands


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
  my ( $self, $qy_species, $qy_assembly, $alignment_type, $compara_db )
    = @_;

  if ( !$self->adaptor() ) {
    warning(
           "Cannot retrieve DnaAlignFeatures without attached adaptor");
    return [];
  }

  unless ( $qy_species && $alignment_type    # && $qy_assembly
    )
  {
    throw(
"Query species and assembly and alignmemt type arguments are required"
    );
  }

  if ( !defined($compara_db) ) {
    $compara_db =
      Bio::EnsEMBL::Registry->get_DBAdaptor( "compara", "compara" );
  }
  unless ($compara_db) {
    warning(
         "Compara database must be attached to core database or passed "
           . "as an argument to "
           . "retrieve compara information" );
    return [];
  }

  my $dafa = $compara_db->get_DnaAlignFeatureAdaptor;

  ## circular BOF
  my ($sl1, $sl2) = $self->_split;
  my ( @arr, @arr1, @arr2 );
  @arr1 = @{
    $dafa->fetch_all_by_Slice( $sl1, $qy_species, $qy_assembly,
                               $alignment_type ) };
  @arr2 = @{
    $dafa->fetch_all_by_Slice( $sl2, $qy_species, $qy_assembly,
                               $alignment_type ) };
  push @arr, @arr1, @arr2;
  return \@arr;
  ## circular EOF

} ## end sub get_all_compara_DnaAlignFeatures

=head2 get_all_compara_Syntenies

  Arg [1]    : string $query_species e.g. "Mus_musculus" or "Mus musculus"
  Arg [2]    : string $method_link_type, default is "SYNTENY"
  Arg [3]    : <optional> compara dbadaptor to use.
  Description: gets all the compara syntenyies for a specfic species
  Returns    : arrayref of Bio::EnsEMBL::Compara::SyntenyRegion
  Status     : Stable

=cut

sub get_all_compara_Syntenies {
  my ( $self, $qy_species, $method_link_type, $compara_db ) = @_;

  if ( !$self->adaptor() ) {
    warning("Cannot retrieve features without attached adaptor");
    return [];
  }

  unless ($qy_species) {
    throw("Query species and assembly arguments are required");
  }

  unless ( defined $method_link_type ) {
    $method_link_type = "SYNTENY";
  }

  if ( !defined($compara_db) ) {
    $compara_db =
      Bio::EnsEMBL::Registry->get_DBAdaptor( "compara", "compara" );
  }
  unless ($compara_db) {
    warning(
         "Compara database must be attached to core database or passed "
           . "as an argument to "
           . "retrieve compara information" );
    return [];
  }
  my $gdba  = $compara_db->get_GenomeDBAdaptor();
  my $mlssa = $compara_db->get_MethodLinkSpeciesSetAdaptor();
  my $dfa   = $compara_db->get_DnaFragAdaptor();
  my $sra   = $compara_db->get_SyntenyRegionAdaptor();

  my $this_gdb =
    $gdba->fetch_by_core_DBAdaptor( $self->adaptor()->db() );
  my $query_gdb = $gdba->fetch_by_registry_name($qy_species);
  my $mlss =
    $mlssa->fetch_by_method_link_type_GenomeDBs( $method_link_type,
                                            [ $this_gdb, $query_gdb ] );

  my $cs = $self->coord_system()->name();
  my $sr = $self->seq_region_name();
  my ($dnafrag) =
    @{ $dfa->fetch_all_by_GenomeDB_region( $this_gdb, $cs, $sr ) };

  ## circular BOF
  my ($sl1, $sl2) = $self->_split;
  my ( @arr, @arr1, @arr2 );
  @arr1 = @{
    $sra->fetch_all_by_MethodLinkSpeciesSet_DnaFrag( $mlss, $dnafrag,
                                                $sl1->start, $sl1->end )
    };
  @arr2 = @{
    $sra->fetch_all_by_MethodLinkSpeciesSet_DnaFrag( $mlss, $dnafrag,
                                                $sl2->start, $sl2->end )
    };
  push @arr, @arr1, @arr2;
  return \@arr;
  ## circular EOF

} ## end sub get_all_compara_Syntenies

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
  my ( $self, $lite_flag ) = @_;

  if ( !$self->adaptor() ) {
    warning("Cannot retrieve features without attached adaptor");
    return [];
  }

  my $haplo_db = $self->adaptor->db->get_db_adaptor('haplotype');

  unless ($haplo_db) {
    warning( "Haplotype database must be attached to core database to "
             . "retrieve haplotype information" );
    return [];
  }

  my $haplo_adaptor = $haplo_db->get_HaplotypeAdaptor;

  ## circular BOF
  my ($sl1, $sl2) = $self->_split;
  my ( @arr, @arr1, @arr2 );
  @arr1 = @{ $haplo_adaptor->fetch_all_by_Slice( $sl1, $lite_flag ) };
  @arr2 = @{ $haplo_adaptor->fetch_all_by_Slice( $sl2, $lite_flag ) };
  push @arr, @arr1, @arr2;
  return \@arr;
  ## circular EOF

} ## end sub get_all_Haplotypes



=head2 get_all_DASFeatures

  Arg [1]    : none
  Example    : $features = $slice->get_all_DASFeatures;
  Description: Retreives a hash reference to a hash of DAS feature
               sets, keyed by the DNS, NOTE the values of this hash
               are an anonymous array containing:
                (1) a pointer to an array of features;
                (2) a pointer to the DAS stylesheet
  Returntype : hashref of Bio::SeqFeatures
  Exceptions : ?
  Caller     : webcode
  Status     : Stable

=cut
sub get_all_DASFeatures {
  my ( $self, $source_type ) = @_;

  if ( !$self->adaptor() ) {
    warning("Cannot retrieve features without attached adaptor");
    return [];
  }

  my %genomic_features = map {
    ( $_->adaptor->dsn =>
      [ $_->fetch_all_Features( $self, $source_type ) ] )
  } $self->adaptor()->db()->_each_DASFeatureFactory;
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
  my ( $self, $track_name ) = @_;
  if ( !$self->adaptor() ) {
    warning("Cannot retrieve features without attached adaptor");
    return [];
  }
  my $features    = [];
  my $xfa_hash    = $self->adaptor->db->get_ExternalFeatureAdaptors;
  my @xf_adaptors = ();
  if ($track_name) {
    #use a specific adaptor
    if ( exists $xfa_hash->{$track_name} ) {
      push @xf_adaptors, $xfa_hash->{$track_name};
    }
  } else {
    #use all of the adaptors
    push @xf_adaptors, values %$xfa_hash;
  }

  ## circular BOF
  my ($sl1, $sl2) = $self->_split;
  foreach my $xfa (@xf_adaptors) {
    push @$features, @{ $xfa->fetch_all_by_Slice($sl1) };
    push @$features, @{ $xfa->fetch_all_by_Slice($sl2) };
  }
  return $features;
  ## circular EOF

} ## end sub get_all_ExternalFeatures

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

  my ( $self, @names ) = @_;

  if ( !$self->adaptor() ) {
    warning('Cannot retrieve features without attached adaptor');
    return [];
  }

  my $db = $self->adaptor()->db();

  my %features = ();    # this will hold the results

  # get the adaptors for each feature
  my %adaptors = %{ $db->get_GenericFeatureAdaptors(@names) };

  foreach my $adaptor_name ( keys(%adaptors) ) {

    my $adaptor_obj = $adaptors{$adaptor_name};
    # get the features and add them to the hash
    ## circular BOF
    my ($sl1, $sl2) = $self->_split;
    my ( @arr1, @arr2 );
    my $features_ref;
    @arr1 = @{ $adaptor_obj->fetch_all_by_Slice($sl1) };
    @arr2 = @{ $adaptor_obj->fetch_all_by_Slice($sl2) };
    push @{$features_ref}, @arr1, @arr2;
    ## circular EOF

    # add each feature to the hash to be returned
    foreach my $feature (@$features_ref) {
      $features{$adaptor_name} = $feature;
    }
  } ## end foreach my $adaptor_name ( ...)

  return \%features;

} ## end sub get_generic_features

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
  my $self     = shift;
  my $to_slice = shift;

  throw('Slice argument is required') if ( !$to_slice );

  my $slice_adaptor = $self->adaptor();

  if ( !$slice_adaptor ) {
    warning("Cannot project without attached adaptor.");
    return [];
  }

  my $mapper_aptr = $slice_adaptor->db->get_AssemblyMapperAdaptor();

  my $cs       = $to_slice->coord_system();
  my $slice_cs = $self->coord_system();

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
    my $asm_mapper = $asma->fetch_by_CoordSystems( $slice_cs, $cs );

    # perform the mapping between this slice and the requested system
    my @coords;

    if ( defined $asm_mapper ) {
      @coords = $asm_mapper->map( $normal_slice->seq_region_name(),
                                  $normal_slice->start(),
                                  $normal_slice->end(),
                                  $normal_slice->strand(),
                                  $slice_cs,
                                  undef,
                                  $to_slice );
    } else {
      $coords[0] =
        Bio::EnsEMBL::Mapper::Gap->new( $normal_slice->start(),
                                        $normal_slice->end() );
    }

    #construct a projection from the mapping results and return it
    foreach my $coord (@coords) {
      my $coord_start = $coord->start();
      my $coord_end   = $coord->end();
      my $length      = $coord_end - $coord_start + 1;

      #skip gaps
      if ( $coord->isa('Bio::EnsEMBL::Mapper::Coordinate') ) {
        my $coord_cs = $coord->coord_system();

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
        my $slice =
          $slice_adaptor->fetch_by_seq_region_id( $coord->id(),
                           $coord_start, $coord_end, $coord->strand() );

        my $current_end = $current_start + $length - 1;

        push @projection,
          bless( [ $current_start, $current_end, $slice ],
                 "Bio::EnsEMBL::ProjectionSegment" );
      }

      $current_start += $length;
    } ## end foreach my $coord (@coords)
  } ## end foreach my $segment (@$normal_slice_proj)

# delete the cache as we may want to map to different set next time and old
# results will be cached.

  $mapper_aptr->delete_cache;

  return \@projection;
} ## end sub project_to_slice

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

sub desc {
  return $_[0]->coord_system->name() . ' ' . $_[0]->seq_region_name();
}

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

=head2 is_circular

  Description: Included for Bio::PrimarySeqI interface compliance (1.2)

=cut

sub is_circular {
  my ($self) = @_;

  if ( !defined( $self->{'circular'} ) ) {
    my @attrs =
      grep { $_ } @{ $self->get_all_Attributes('circular_seq') };
    $self->{'circular'} = @attrs ? 1 : 0;
  }

  return $self->{'circular'};
}

# sub DEPRECATED METHODS #
###############################################################################

=head1 DEPRECATED METHODS

=cut

=head2 chr_name

  DEPRECATED use seq_region_name() instead

=cut

sub chr_name {
  deprecate("Use seq_region_name() instead");
  seq_region_name(@_);
}

=head2 chr_start

  DEPRECATED use start() instead

=cut

sub chr_start {
  deprecate('Use start() instead');
  start(@_);
}

=head2 chr_end

  DEPRECATED use end() instead

  Returntype : int
  Exceptions : none
  Caller     : SliceAdaptor, general

=cut

sub chr_end {
  deprecate('Use end() instead');
  end(@_);
}

1;
