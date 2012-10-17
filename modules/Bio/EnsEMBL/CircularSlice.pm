=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
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
  $s{'strand'} = $self->{'strand'}*-1;

  # reverse compliment any attached sequence
  reverse_comp( \$s{'seq'} ) if ( $s{'seq'} );

  # bless and return the copy
  return bless \%s, ref $self;
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

      my $sl1 =
        Bio::EnsEMBL::Slice->new(
                     -COORD_SYSTEM      => $self->{'coord_system'},
                     -SEQ_REGION_NAME   => $self->{'seq_region_name'},
                     -SEQ_REGION_LENGTH => $self->{'seq_region_length'},
                     -START             => $self->{'start'},
                     -END               => $self->{'seq_region_length'},
                     -STRAND            => $self->{'strand'},
                     -ADAPTOR           => $self->adaptor() );

      my $sl2 =
        Bio::EnsEMBL::Slice->new(
                     -COORD_SYSTEM      => $self->{'coord_system'},
                     -SEQ_REGION_NAME   => $self->{'seq_region_name'},
                     -SEQ_REGION_LENGTH => $self->{'seq_region_length'},
                     -START             => 1,
                     -END               => $self->{'end'},
                     -STRAND            => $self->{'strand'},
                     -ADAPTOR           => $self->adaptor() );

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

  my $sl01 =
    Bio::EnsEMBL::CircularSlice->new(
                     -COORD_SYSTEM      => $self->{'coord_system'},
                     -SEQ_REGION_NAME   => $self->{'seq_region_name'},
                     -SEQ_REGION_LENGTH => $self->{'seq_region_length'},
                     -START             => $self->{'start'},
                     -END               => $self->{'seq_region_length'},
                     -STRAND            => $self->{'strand'},
                     -ADAPTOR           => $self->adaptor() );

  my $sl02 =
    Bio::EnsEMBL::CircularSlice->new(
                     -COORD_SYSTEM      => $self->{'coord_system'},
                     -SEQ_REGION_NAME   => $self->{'seq_region_name'},
                     -SEQ_REGION_LENGTH => $self->{'seq_region_length'},
                     -START             => 1,
                     -END               => $self->{'end'},
                     -STRAND            => $self->{'strand'},
                     -ADAPTOR           => $self->adaptor() );

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

=head2 sub_Slice_same

  Arg   1    : int $start
  Arg   2    : int $end
  Arge [3]   : int $strand
  Description: Makes another Slice that covers only part of this slice
               If a slice is requested which lies outside of the boundaries
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

sub sub_Slice_same {
  my ( $self, $start, $end, $strand ) = @_;

  if ( $start < 1 || $start > $self->{'end'} ) {
    # throw( "start argument not valid" );
    return undef;
  }

  if ( $end < $start || $end > $self->{'end'} ) {
    # throw( "end argument not valid" )
    return undef;
  }

  my ( $new_start, $new_end, $new_strand, $new_seq );
  if ( !defined $strand ) {
    $strand = 1;
  }

  if ( $self->{'strand'} == 1 ) {
    $new_start  = $self->{'start'} + $start - 1;
    $new_end    = $self->{'start'} + $end - 1;
    $new_strand = $strand;
  } else {
    $new_start  = $self->{'end'} - $end + 1;
    $new_end    = $self->{'end'} - $start + 1;
    $new_strand = -$strand;
  }

  if ( defined $self->{'seq'} ) {
    $new_seq = $self->subseq( $start, $end, $strand );
  }

  #fastest way to copy a slice is to do a shallow hash copy
  my %new_slice = %$self;
  $new_slice{'start'}  = int($new_start);
  $new_slice{'end'}    = int($new_end);
  $new_slice{'strand'} = $new_strand;
  if ($new_seq) {
    $new_slice{'seq'} = $new_seq;
  }

  return bless \%new_slice, ref($self);
} ## end sub sub_Slice_same

=head2 seq_region_Slice_same

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

sub seq_region_Slice_same {
  my $self = shift;

  if ( $self->{'seq'} ) {
    warning(
          "Cannot get a seq_region_Slice of a slice which has manually "
            . "attached sequence " );
    return undef;
  }

  # quick shallow copy
  my $slice;
  %{$slice} = %{$self};
  bless $slice, ref($self);

  $slice->{'start'}  = 1;
  $slice->{'end'}    = $slice->{'seq_region_length'};
  $slice->{'strand'} = 1;

  return $slice;
}

=head2 get_seq_region_id_same

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

sub get_seq_region_id_same {
  my ($self) = @_;

  if ( $self->adaptor ) {
    return $self->adaptor->get_seq_region_id($self);
  } else {
    warning('Cannot retrieve seq_region_id without attached adaptor.');
    return undef;
  }
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
  my $self        = shift;
  my $attrib_code = shift;

  my $result;
  my @results;

  if ( !$self->adaptor() ) {
    warning('Cannot get attributes without an adaptor.');
    return [];
  }

  my $attribute_adaptor = $self->adaptor->db->get_AttributeAdaptor();

## circular BOF
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
  my $pta_ref = [];

  if ( defined $attrib_code ) {
    my @res1 =
      grep { uc( $_->code() ) eq uc($attrib_code) }
      @{ $attribute_adaptor->fetch_all_by_Slice($sl1) };
    my @res2 =
      grep { uc( $_->code() ) eq uc($attrib_code) }
      @{ $attribute_adaptor->fetch_all_by_Slice($sl2) };
    my @res;
    push @res, @res1, @res2;
    $result = \@res;
  } else {
    my @res1 = @{ $attribute_adaptor->fetch_all_by_Slice($sl1) };
    my @res2 = @{ $attribute_adaptor->fetch_all_by_Slice($sl2) };
    my @res;
    push @res, @res1, @res2;
    $result = \@res;
  }

## circular EOF

  #  if( defined $attrib_code ) {
  #    @results = grep { uc($_->code()) eq uc($attrib_code) }
  #      @{$attribute_adaptor->fetch_all_by_Slice( $self )};
  #    $result = \@results;
  #  } else {
  #    $result = $attribute_adaptor->fetch_all_by_Slice( $self );
  #  }

  return $result;
} ## end sub get_all_Attributes

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
  my ( $self, $logic_name, $load_exons ) = @_;

  if ( !$self->adaptor() ) {
    warning(
           'Cannot get PredictionTranscripts without attached adaptor');
    return [];
  }
  my $pta = $self->adaptor()->db()->get_PredictionTranscriptAdaptor();

## circular BOF
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

  my ( @arr, @arr1, @arr2 );
  @arr1 =
    @{ $pta->fetch_all_by_Slice( $sl1, $logic_name, $load_exons ) };
  @arr2 =
    @{ $pta->fetch_all_by_Slice( $sl2, $logic_name, $load_exons ) };
  push @arr, @arr1, @arr2;
  return \@arr;
## circular EOF

  #return $pta->fetch_all_by_Slice($sl1, $logic_name, $load_exons);
} ## end sub get_all_PredictionTranscripts

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
               retrieved.
  Returntype : listref of Bio::EnsEMBL::DnaDnaAlignFeatures
  Exceptions : warning if slice does not have attached adaptor
  Caller     : general
  Status     : Stable

=cut

sub get_all_DnaAlignFeatures {
  my ( $self, $logic_name, $score, $dbtype, $hcoverage ) = @_;

  if ( !$self->adaptor() ) {
    warning('Cannot get DnaAlignFeatures without attached adaptor');
    return [];
  }

  my $db;

  if ($dbtype) {
    $db = $self->adaptor->db->get_db_adaptor($dbtype);
    if ( !$db ) {
      warning("Don't have db $dbtype returning empty list\n");
      return [];
    }
  } else {
    $db = $self->adaptor->db;
  }

  my $dafa = $db->get_DnaAlignFeatureAdaptor();

  if ( defined($score) and defined($hcoverage) ) {
    warning "cannot specify score and hcoverage. Using score only";
  }

## circular BOF
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

  my ( @arr, @arr1, @arr2 );
  my %union;
  if ( defined($score) ) {
    @arr1 = @{ $dafa->fetch_all_by_Slice_and_score( $sl1, $score,
                                                    $logic_name ) };
    @arr2 = @{ $dafa->fetch_all_by_Slice_and_score( $sl2, $score,
                                                    $logic_name ) };
    push @arr, @arr1, @arr2;
    return \@arr;
  }
  @arr1 = @{
    $dafa->fetch_all_by_Slice_and_hcoverage( $sl1, $hcoverage,
                                             $logic_name ) };
  @arr2 = @{
    $dafa->fetch_all_by_Slice_and_hcoverage( $sl2, $hcoverage,
                                             $logic_name ) };
  push @arr, @arr1, @arr2;
  return \@arr;

## circular EOF

# if(defined($score)){
#   return $dafa->fetch_all_by_Slice_and_score($self,$score, $logic_name);
# }
# return $dafa->fetch_all_by_Slice_and_hcoverage($self,$hcoverage, $logic_name);
} ## end sub get_all_DnaAlignFeatures

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
  my ( $self, $logic_name, $score, $dbtype, $hcoverage ) = @_;
  if ( !$self->adaptor() ) {
    warning('Cannot get ProteinAlignFeatures without attached adaptor');
    return [];
  }
  my $db;
  if ($dbtype) {
    $db = $self->adaptor->db->get_db_adaptor($dbtype);
    if ( !$db ) {
      warning("Don't have db $dbtype returning empty list\n");
      return [];
    }
  } else {
    $db = $self->adaptor->db;
  }

  my $pafa = $db->get_ProteinAlignFeatureAdaptor();

  if ( defined($score) and defined($hcoverage) ) {
    warning "cannot specify score and hcoverage. Using score only";
  }

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
  my ( @arr, @arr1, @arr2 );
  if ( defined($score) ) {
    @arr1 = @{ $pafa->fetch_all_by_Slice_and_score( $sl1, $score,
                                                    $logic_name ) };
    @arr2 = @{ $pafa->fetch_all_by_Slice_and_score( $sl2, $score,
                                                    $logic_name ) };

    push @arr, @arr1, @arr2;
    return \@arr;
  }

  @arr1 = @{
    $pafa->fetch_all_by_Slice_and_hcoverage( $sl1, $hcoverage,
                                             $logic_name ) };
  @arr2 = @{
    $pafa->fetch_all_by_Slice_and_hcoverage( $sl2, $hcoverage,
                                             $logic_name ) };

  push @arr, @arr1, @arr2;
  return \@arr;

#  if(defined($score)){
#    return $pafa->fetch_all_by_Slice_and_score($self,$score, $logic_name);
#  }
#  return $pafa->fetch_all_by_Slice_and_hcoverage($self,$hcoverage, $logic_name);

} ## end sub get_all_ProteinAlignFeatures

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
  my ( $self, $logic_name, $score ) = @_;

  my @out = ();

  push @out,
    @{ $self->get_all_ProteinAlignFeatures( $logic_name, $score ) };
  push @out,
    @{ $self->get_all_DnaAlignFeatures( $logic_name, $score ) };

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
  my ( $self, $logic_name, $score, $dbtype ) = @_;

  if ( !$self->adaptor() ) {
    warning('Cannot get SimpleFeatures without attached adaptor');
    return [];
  }

  my $db;
  if ($dbtype) {
    $db = $self->adaptor->db->get_db_adaptor($dbtype);
    if ( !$db ) {
      warning("Don't have db $dbtype returning empty list\n");
      return [];
    }
  } else {
    $db = $self->adaptor->db;
  }

  my $sfa = $db->get_SimpleFeatureAdaptor();

  return $sfa->fetch_all_by_Slice_and_score( $self, $score,
                                             $logic_name );
}

#### STOP !!! ########################################

=head2 get_all_RepeatFeatures

  Arg [1]    : (optional) string $logic_name
               The name of the analysis performed on the repeat features
               to obtain.
  Arg [2]    : (optional) string $repeat_type
               Limits features returned to those of the specified repeat_type
  Arg [3]    : (optional) string $db
               Key for database e.g. core/vega/cdna/....
  Example    : @repeat_feats = @{$slice->get_all_RepeatFeatures(undef,'LTR')};
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
  my ( $self, $logic_name, $repeat_type, $dbtype ) = @_;

  if ( !$self->adaptor() ) {
    warning('Cannot get RepeatFeatures without attached adaptor');
    return [];
  }

  my $db;
  if ($dbtype) {
    $db = $self->adaptor->db->get_db_adaptor($dbtype);
    if ( !$db ) {
      warning("Don't have db $dbtype returning empty list\n");
      return [];
    }
  } else {
    $db = $self->adaptor->db;
  }

  my $rpfa = $db->get_RepeatFeatureAdaptor();

  return $rpfa->fetch_all_by_Slice( $self, $logic_name, $repeat_type );
}

=head2 get_all_LD_values

    Arg [1]     : (optional) Bio::EnsEMBL::Variation::Population $population
    Description : returns all LD values on this slice. This function will only work correctly if the variation
                  database has been attached to the core database. If the argument is passed, will return the LD information
                  in that population
    ReturnType  : Bio::EnsEMBL::Variation::LDFeatureContainer
    Exceptions  : none
    Caller      : contigview, snpview
     Status     : At Risk
                : Variation database is under development.

=cut

sub get_all_LD_values {
  my $self       = shift;
  my $population = shift;

  if ( !$self->adaptor() ) {
    warning('Cannot get LDFeatureContainer without attached adaptor');
    return [];
  }

  my $variation_db = $self->adaptor->db->get_db_adaptor('variation');

  unless ($variation_db) {
    warning( "Variation database must be attached to core database to "
             . "retrieve variation information" );
    return [];
  }

  my $ld_adaptor = $variation_db->get_LDFeatureContainerAdaptor;

  if ($ld_adaptor) {
    return $ld_adaptor->fetch_by_Slice( $self, $population );
  } else {
    return [];

  }

#     my $ld_adaptor = Bio::EnsEMBL::DBSQL::MergedAdaptor->new(-species => $self->adaptor()->db()->species, -type => "LDFeatureContainer");

#   if( $ld_adaptor ) {
#       my $ld_values = $ld_adaptor->fetch_by_Slice($self,$population);
#       if (@{$ld_values} > 1){
# 	  warning("More than 1 variation database attached. Trying to merge LD results");
# 	  my $ld_value_merged = shift @{$ld_values};
# 	  #with more than 1 variation database attached, will try to merge in one single LDContainer object.
# 	  foreach my $ld (@{$ld_values}){
# 	      #copy the ld values to the result hash
# 	      foreach my $key (keys %{$ld->{'ldContainer'}}){
# 		  $ld_value_merged->{'ldContainer'}->{$key} = $ld->{'ldContainer'}->{$key};
# 	      }
# 	      #and copy the variationFeatures as well
# 	      foreach my $key (keys %{$ld->{'variationFeatures'}}){
# 		  $ld_value_merged->{'variationFeatures'}->{$key} = $ld->{'variationFeatures'}->{$key};
# 	      }

# 	  }
# 	  return $ld_value_merged;
#       }
#       else{
# 	  return shift @{$ld_values};
#       }
# } else {
#     warning("Variation database must be attached to core database to " .
# 		"retrieve variation information" );
#     return [];
# }
} ## end sub get_all_LD_values

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
} ## end sub get_all_VariationFeatures

=head2 get_all_IndividualSlice

    Args        : none
    Example     : my $individualSlice = $slice->get_by_Population($population);
    Description : Gets the specific Slice for all the individuls in the population
    ReturnType  : listref of Bio::EnsEMB::IndividualSlice
    Exceptions  : none
    Caller      : general

=cut

sub get_all_IndividualSlice {
  my $self = shift;

  my $individualSliceFactory =
    Bio::EnsEMBL::IndividualSliceFactory->new(
                     -START             => $self->{'start'},
                     -END               => $self->{'end'},
                     -STRAND            => $self->{'strand'},
                     -ADAPTOR           => $self->adaptor(),
                     -SEQ_REGION_NAME   => $self->{'seq_region_name'},
                     -SEQ_REGION_LENGTH => $self->{'seq_region_length'},
                     -COORD_SYSTEM      => $self->{'coord_system'}, );
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

sub get_by_Individual {
  my $self       = shift;
  my $individual = shift;

  return Bio::EnsEMBL::IndividualSlice->new(
    -START   => $self->{'start'},
    -END     => $self->{'end'},
    -STRAND  => $self->{'strand'},
    -ADAPTOR => $self->adaptor(),
    #					  -SEQ     => $self->{'seq'},
    -SEQ_REGION_NAME   => $self->{'seq_region_name'},
    -SEQ_REGION_LENGTH => $self->{'seq_region_length'},
    -COORD_SYSTEM      => $self->{'coord_system'},
    -INDIVIDUAL        => $individual );

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
  my $self        = shift;
  my $strain_name = shift;

  return
    Bio::EnsEMBL::StrainSlice->new(
                     -START             => $self->{'start'},
                     -END               => $self->{'end'},
                     -STRAND            => $self->{'strand'},
                     -ADAPTOR           => $self->adaptor(),
                     -SEQ               => $self->{'seq'},
                     -SEQ_REGION_NAME   => $self->{'seq_region_name'},
                     -SEQ_REGION_LENGTH => $self->{'seq_region_length'},
                     -COORD_SYSTEM      => $self->{'coord_system'},
                     -STRAIN_NAME       => $strain_name );

}

sub calculate_theta {
  my $self    = shift;
  my $strains = shift;
  my $feature = shift
    ; #optional parameter. Name of the feature in the Slice you want to calculate

  if ( !$self->adaptor() ) {
    warning('Cannot get variation features without attached adaptor');
    return 0;
  }
  my $variation_db = $self->adaptor->db->get_db_adaptor('variation');

  unless ($variation_db) {
    warning( "Variation database must be attached to core database to "
             . "retrieve variation information" );
    return 0;
  }

  #need to get coverage regions for the slice in the different strains
  my $coverage_adaptor = $variation_db->get_ReadCoverageAdaptor;
  my $strain;
  my $differences = [];
  my $slices      = [];
  if ($coverage_adaptor) {
    my $num_strains = scalar( @{$strains} ) + 1;
    if ( !defined $feature ) {
      #we want to calculate for the whole slice
      push @{$slices},
        $self;  #add the slice as the slice to calculate the theta value
    } else {
      #we have features, get the slices for the different features
      my $features = $self->get_all_Exons();
      map { push @{$slices}, $_->feature_Slice }
        @{$features};    #add the slices of the features
    }
    my $length_regions = 0;
    my $snps           = 0;
    my $theta          = 0;
    my $last_position  = 0;
    #get all the differences in the slice coordinates
    foreach my $strain_name ( @{$strains} ) {
      my $strain = $self->get_by_strain($strain_name)
        ;                #get the strainSlice for the strain

      my $results = $strain->get_all_differences_Slice;
      push @{$differences}, @{$results} if ( defined $results );
    }
#when we finish, we have, in max_level, the regions covered by all the sample
#sort the differences by the genomic position
    my @differences_sorted =
      sort { $a->start <=> $b->start } @{$differences};
    foreach my $slice ( @{$slices} ) {
      my $regions_covered =
        $coverage_adaptor->fetch_all_regions_covered( $slice,
                                                      $strains );
      if ( defined $regions_covered ) {
        foreach my $range ( @{$regions_covered} ) {
          $length_regions +=
            ( $range->[1] - $range->[0] ) +
            1;    #add the length of the genomic region
          for ( my $i = $last_position; $i < @differences_sorted; $i++ )
          {
            if (    $differences_sorted[$i]->start >= $range->[0]
                 && $differences_sorted[$i]->end <= $range->[1] )
            {
              $snps++;    #count differences in the region
            } elsif ( $differences_sorted[$i]->end > $range->[1] ) {
              $last_position = $i;
              last;
            }
          }
        }
        #when all the ranges have been iterated, calculate rho
        #this is an intermediate variable called a in the formula
        #  a = sum i=2..strains 1/i-1
      }
    } ## end foreach my $slice ( @{$slices...})
    my $a = _calculate_a($num_strains);
    $theta = $snps/( $a*$length_regions );
    return $theta;
  } else {
    return 0;
  }
} ## end sub calculate_theta

sub _calculate_a {
  my $max_level = shift;

  my $a = 0;
  for ( my $i = 2; $i <= $max_level + 1; $i++ ) {
    $a += 1/( $i - 1 );
  }
  return $a;
}

sub calculate_pi {
  my $self    = shift;
  my $strains = shift;
  my $feature = shift;

  if ( !$self->adaptor() ) {
    warning('Cannot get variation features without attached adaptor');
    return 0;
  }
  my $variation_db = $self->adaptor->db->get_db_adaptor('variation');

  unless ($variation_db) {
    warning( "Variation database must be attached to core database to "
             . "retrieve variation information" );
    return 0;
  }

  #need to get coverage regions for the slice in the different strains
  my $coverage_adaptor = $variation_db->get_ReadCoverageAdaptor;
  my $differences      = [];
  my $slices           = [];
  if ($coverage_adaptor) {
    my $num_strains = scalar( @{$strains} ) + 1;
    if ( !defined $feature ) {
      #we want to calculate for the whole slice
      push @{$slices},
        $self;  #add the slice as the slice to calculate the theta value
    } else {
      #we have features, get the slices for the different features
      my $features = $self->get_all_Exons();
      map { push @{$slices}, $_->feature_Slice }
        @{$features};    #add the slices of the features
    }
    my @range_differences = ();
    my $pi                = 0;
    my $regions           = 0;
    my $last_position =
      0;    #last position visited in the sorted list of differences
    my $triallelic    = 0;
    my $is_triallelic = 0;
    foreach my $slice ( @{$slices} ) {

      foreach my $strain_name ( @{$strains} ) {
        my $strain = $slice->get_by_strain($strain_name)
          ;    #get the strainSlice for the strain
        my $results = $strain->get_all_differences_Slice;
        push @{$differences}, @{$results} if ( defined $results );
      }
      my @differences_sorted =
        sort { $a->start <=> $b->start } @{$differences};

      my $regions_covered =
        $coverage_adaptor->fetch_all_regions_covered( $slice,
                                                      $strains );
#when we finish, we have, in max_level, the regions covered by all the sample
#sort the differences
      if ( defined $regions_covered ) {
        foreach my $range ( @{$regions_covered} ) {
          for ( my $i = $last_position; $i < @differences_sorted; $i++ )
          {
            if (    $differences_sorted[$i]->start >= $range->[0]
                 && $differences_sorted[$i]->end <= $range->[1] )
            {
              #check wether it is the same region or different
              if ( !defined $range_differences[0]
                   || ( $differences_sorted[$i]->start ==
                        $range_differences[0]->start ) )
              {
                if ( defined $range_differences[0]
                     && ( $differences_sorted[$i]->allele_string ne
                          $range_differences[0]->allele_string ) )
                {
                  $is_triallelic = 1;
                }
                push @range_differences, $differences_sorted[$i];
              } else {
                #new site, calc pi for the previous one
                $pi +=
                  2*
                  ( @range_differences/($num_strains) )*
                  ( 1 - ( @range_differences/$num_strains ) );
                if ($is_triallelic) {
                  $triallelic++;
                  $is_triallelic = 0;
                }
                $regions++;
                @range_differences = ();
                #and start a new range
                push @range_differences, $differences_sorted[$i];
              }
            } elsif ( $differences_sorted[$i]->end > $range->[1] ) {
              $last_position = $i;
              last;
            }
          } ## end for ( my $i = $last_position...)
          #calculate pi for last site, if any
          if ( defined $range_differences[0] ) {
            $pi +=
              2*
              ( @range_differences/$num_strains )*
              ( 1 - ( @range_differences/$num_strains ) );
            $regions++;
          }
        } ## end foreach my $range ( @{$regions_covered...})
      } ## end if ( defined $regions_covered)
      $pi = $pi/$regions;    #calculate average pi
      print
"Regions with variations in region $regions and triallelic $triallelic\n\n";
    } ## end foreach my $slice ( @{$slices...})
    return $pi;
  } else {
    return 0;
  }

} ## end sub calculate_pi

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

=head2 get_all_SNPs

 Description: DEPRECATED. Use get_all_VariationFeatures instead

=cut

sub get_all_SNPs {
  my $self = shift;

  deprecate('Use get_all_VariationFeatures() instead.');

  my $snps;
  my $vf = $self->get_all_genotyped_VariationFeatures();
  if ( $vf->[0] ) {
    #necessary to convert the VariationFeatures into SNP objects
    foreach my $variation_feature ( @{$vf} ) {
      push @{$snps}, $variation_feature->convert_to_SNP();
    }
    return $snps;
  } else {
    return [];
  }
}

=head2 get_all_genotyped_SNPs

  Description   : DEPRECATED. Use get_all_genotyped_VariationFeatures instead

=cut

sub get_all_genotyped_SNPs {
  my $self = shift;

  deprecate("Use get_all_genotyped_VariationFeatures instead");
  my $vf = $self->get_all_genotyped_VariationFeatures;
  my $snps;
  if ( $vf->[0] ) {
    foreach my $variation_feature ( @{$vf} ) {
      push @{$snps}, $variation_feature->convert_to_SNP();
    }
    return $snps;
  } else {
    return [];
  }
}

sub get_all_SNPs_transcripts {
  my $self = shift;

  deprecate("DEPRECATED");

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
  Description: Retrieves all genes that overlap this slice.
  Returntype : listref of Bio::EnsEMBL::Genes
  Exceptions : none
  Caller     : none
  Status     : Stable

=cut

sub get_all_Genes {
  my ( $self, $logic_name, $dbtype, $load_transcripts, $source,
       $biotype )
    = @_;

  if ( !$self->adaptor() ) {
    warning('Cannot get Genes without attached adaptor');
    return [];
  }

  my $ga;
  if ($dbtype) {
    my $db = $reg->get_db( $self->adaptor()->db(), $dbtype );
    if ( defined($db) ) {
      $ga = $reg->get_adaptor( $db->species(), $db->group(), "Gene" );
    } else {
      $ga = $reg->get_adaptor( $self->adaptor()->db()->species(),
                               $dbtype, "Gene" );
    }
    if ( !defined $ga ) {
      warning("$dbtype genes not available");
      return [];
    }
  } else {
    $ga = $self->adaptor->db->get_GeneAdaptor();
  }

## circular BOF
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

  my ( @arr, @arr1, @arr2 );
  @arr1 = @{
    $ga->fetch_all_by_Slice( $sl1, $logic_name, $load_transcripts,
                             $source, $biotype ) };
  @arr2 = @{
    $ga->fetch_all_by_Slice( $sl2, $logic_name, $load_transcripts,
                             $source, $biotype ) };
  push @arr, @arr1, @arr2;
  return \@arr;

  ## circular EOF

} ## end sub get_all_Genes

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
               If not provided all genes will be retrieved instead.

  Returntype : listref of Bio::EnsEMBL::Genes
  Exceptions : none
  Caller     : genebuilder, general
  Status     : Stable

=cut

sub get_all_Genes_by_type {
  my ( $self, $type, $logic_name, $load_transcripts ) = @_;

  if ( !$self->adaptor() ) {
    warning('Cannot get Genes without attached adaptor');
    return [];
  }

  return
    $self->get_all_Genes( $logic_name, undef, $load_transcripts, undef,
                          $type );
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

  Returntype : listref of Bio::EnsEMBL::Genes
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Genes_by_source {
  my ( $self, $source, $load_transcripts ) = @_;

  if ( !$self->adaptor() ) {
    warning('Cannot get Genes without attached adaptor');
    return [];
  }

  return $self->get_all_Genes( undef, undef, $load_transcripts,
                               $source );
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
               through the transcripts of each gene.
  Returntype : reference to a list of Bio::EnsEMBL::Transcripts
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Transcripts {
  my $self       = shift;
  my $load_exons = shift;
  my $logic_name = shift;
  my $dbtype     = shift;
  if ( !$self->adaptor() ) {
    warning('Cannot get Transcripts without attached adaptor');
    return [];
  }

  my $ta;
  if ($dbtype) {
    my $db = $reg->get_db( $self->adaptor()->db(), $dbtype );
    if ( defined($db) ) {
      $ta =
        $reg->get_adaptor( $db->species(), $db->group(), "Transcript" );
    } else {
      $ta = $reg->get_adaptor( $self->adaptor()->db()->species(),
                               $dbtype, "Transcript" );
    }
    if ( !defined $ta ) {
      warning("$dbtype genes not available");
      return [];
    }
  } else {
    $ta = $self->adaptor->db->get_TranscriptAdaptor();
  }

## circular BOF
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

  my ( @arr, @arr1, @arr2 );
  @arr1 =
    @{ $ta->fetch_all_by_Slice( $sl1, $load_exons, $logic_name ) };
  @arr2 =
    @{ $ta->fetch_all_by_Slice( $sl2, $load_exons, $logic_name ) };
  push @arr, @arr1, @arr2;
  return \@arr;
## circular EOF
  #return $ta->fetch_all_by_Slice($self, $load_exons, $logic_name);
} ## end sub get_all_Transcripts

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

  if ( !$self->adaptor() ) {
    warning('Cannot get Exons without attached adaptor');
    return [];
  }

  ## circular BOF
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

  my ( @arr, @arr1, @arr2 );
  @arr1 =
    @{ $sl1->adaptor->db->get_ExonAdaptor->fetch_all_by_Slice($sl1) };
  @arr2 =
    @{ $sl2->adaptor->db->get_ExonAdaptor->fetch_all_by_Slice($sl2) };
  push @arr, @arr1, @arr2;
  return \@arr;
  ## circular EOF

  #rturn $self->adaptor->db->get_ExonAdaptor->fetch_all_by_Slice($self);
} ## end sub get_all_Exons

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

  my ( @arr, @arr1, @arr2 );
  @arr1 = @{ $qfAdaptor->fetch_all_by_Slice_constraint($sl1) };
  @arr2 = @{ $qfAdaptor->fetch_all_by_Slice_constraint($sl2) };
  push @arr, @arr1, @arr2;
  return \@arr;
  ## circular EOF

  #return $qfAdaptor->fetch_all_by_Slice_constraint( $self );
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

  my ( @arr, @arr1, @arr2 );
  @arr1 = @{ $kadp->fetch_all_by_Slice($sl1) };
  @arr2 = @{ $kadp->fetch_all_by_Slice($sl2) };
  push @arr, @arr1, @arr2;
  return \@arr;
  ## circular EOF

  #return $kadp->fetch_all_by_Slice($self);
} ## end sub get_all_KaryotypeBands

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
  my ( $self, $logic_names, $soft_mask, $not_default_masking_cases ) =
    @_;

  return
    Bio::EnsEMBL::RepeatMaskedSlice->new(
                -START             => $self->{'start'},
                -END               => $self->{'end'},
                -STRAND            => $self->{'strand'},
                -ADAPTOR           => $self->adaptor(),
                -SEQ               => $self->{'seq'},
                -SEQ_REGION_NAME   => $self->{'seq_region_name'},
                -SEQ_REGION_LENGTH => $self->{'seq_region_length'},
                -COORD_SYSTEM      => $self->{'coord_system'},
                -REPEAT_MASK       => $logic_names,
                -SOFT_MASK         => $soft_mask,
                -NOT_DEFAULT_MASKING_CASES => $not_default_masking_cases
    );
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
  my ( $self, $dnaref, $repeats, $soft_mask,
       $not_default_masking_cases ) = @_;

  $soft_mask = 0 unless ( defined $soft_mask );
  $not_default_masking_cases = {}
    unless ( defined $not_default_masking_cases );

  # explicit CORE::length call, to avoid any confusion with the Slice
  # length method
  my $dnalen = CORE::length($$dnaref);

REP: foreach my $old_f ( @{$repeats} ) {
    my $f      = $old_f->transfer($self);
    my $start  = $f->start;
    my $end    = $f->end;
    my $length = ( $end - $start ) + 1;

    # check if we get repeat completely outside of expected slice range
    if ( $end < 1 || $start > $dnalen ) {
  # warning("Unexpected: Repeat completely outside slice coordinates.");
      next REP;
    }

    # repeat partly outside slice range, so correct
    # the repeat start and length to the slice size if needed
    if ( $start < 1 ) {
      $start  = 1;
      $length = ( $end - $start ) + 1;
    }

    # repeat partly outside slice range, so correct
    # the repeat end and length to the slice size if needed
    if ( $end > $dnalen ) {
      $end    = $dnalen;
      $length = ( $end - $start ) + 1;
    }

    $start--;

    my $padstr;
# if we decide to define masking on the base of the repeat_type, we'll need
# to add the following, and the other commented line few lines below.
# my $rc_type = "repeat_type_" . $f->repeat_consensus->repeat_type;
    my $rc_class = "repeat_class_" . $f->repeat_consensus->repeat_class;
    my $rc_name  = "repeat_name_" . $f->repeat_consensus->name;

    my $masking_type;
# $masking_type = $not_default_masking_cases->{$rc_type} if (defined $not_default_masking_cases->{$rc_type});
    $masking_type = $not_default_masking_cases->{$rc_class}
      if ( defined $not_default_masking_cases->{$rc_class} );
    $masking_type = $not_default_masking_cases->{$rc_name}
      if ( defined $not_default_masking_cases->{$rc_name} );

    $masking_type = $soft_mask unless ( defined $masking_type );

    if ($masking_type) {
      $padstr = lc substr( $$dnaref, $start, $length );
    } else {
      $padstr = 'N' x $length;
    }
    substr( $$dnaref, $start, $length ) = $padstr;
  } ## end foreach my $old_f ( @{$repeats...})
} ## end sub _mask_features

=head2 get_all_SearchFeatures

  Arg [1]    : scalar $ticket_ids
  Example    : $slice->get_all_SearchFeatures('BLA_KpUwwWi5gY');
  Description: Retreives all search features for stored blast
               results for the ticket that overlap this slice
  Returntype : listref of Bio::EnsEMBL::SeqFeatures
  Exceptions : none
  Caller     : general (webby!)
  Status     : Stable

=cut

sub get_all_SearchFeatures {
  my $self   = shift;
  my $ticket = shift;
  local $_;
  unless ($ticket) {
    throw("ticket argument is required");
  }

  if ( !$self->adaptor() ) {
    warning("Cannot get SearchFeatures without an attached adaptor");
    return [];
  }

  my $sfa = $self->adaptor()->db()->get_db_adaptor('blast');

  my $offset = $self->start - 1;

  my $features = $sfa
    ? $sfa->get_all_SearchFeatures( $ticket, $self->seq_region_name,
                                    $self->start, $self->end )
    : [];

  foreach (@$features) {
    $_->start( $_->start - $offset );
    $_->end( $_->end - $offset );
  }
  return $features;

} ## end sub get_all_SearchFeatures

=head2 get_all_AssemblyExceptionFeatures

  Arg [1]    : string $set (optional)
  Example    : $slice->get_all_AssemblyExceptionFeatures();
  Description: Retreives all misc features which overlap this slice. If
               a set code is provided only features which are members of
               the requested set are returned.
  Returntype : listref of Bio::EnsEMBL::AssemblyExceptionFeatures
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_AssemblyExceptionFeatures {
  my $self     = shift;
  my $misc_set = shift;

  my $adaptor = $self->adaptor();

  if ( !$adaptor ) {
    warning('Cannot retrieve features without attached adaptor.');
    return [];
  }

  my $aefa = $adaptor->db->get_AssemblyExceptionFeatureAdaptor();

  ## circular BOF
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

  my ( @arr, @arr1, @arr2 );
  @arr1 = @{ $aefa->fetch_all_by_Slice($sl1) };
  @arr2 = @{ $aefa->fetch_all_by_Slice($sl2) };
  push @arr, @arr1, @arr2;
  return \@arr;
  ## circular EOF
  #return $aefa->fetch_all_by_Slice($self);
} ## end sub get_all_AssemblyExceptionFeatures

=head2 get_all_MiscFeatures

  Arg [1]    : string $set (optional)
  Arg [2]    : string $database (optional)
  Example    : $slice->get_all_MiscFeatures('cloneset');
  Description: Retreives all misc features which overlap this slice. If
               a set code is provided only features which are members of
               the requested set are returned.
  Returntype : listref of Bio::EnsEMBL::MiscFeatures
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_MiscFeatures {
  my $self     = shift;
  my $misc_set = shift;
  my $dbtype   = shift;
  my $msa;

  my $adaptor = $self->adaptor();
  if ( !$adaptor ) {
    warning('Cannot retrieve features without attached adaptor.');
    return [];
  }

  my $mfa;
  if ($dbtype) {
    my $db = $reg->get_db( $adaptor->db(), $dbtype );
    if ( defined($db) ) {
      $mfa = $reg->get_adaptor( lc( $db->species() ),
                                $db->group(), "miscfeature" );
    } else {
      $mfa = $reg->get_adaptor( $adaptor->db()->species(),
                                $dbtype, "miscfeature" );
    }
    if ( !defined $mfa ) {
      warning("$dbtype misc features not available");
      return [];
    }
  } else {
    $mfa = $adaptor->db->get_MiscFeatureAdaptor();
  }

  ## circular BOF
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

  my ( @arr, @arr1, @arr2 );

  if ($misc_set) {
    @arr1 =
      @{ $mfa->fetch_all_by_Slice_and_set_code( $sl1, $misc_set ) };
    @arr2 =
      @{ $mfa->fetch_all_by_Slice_and_set_code( $sl2, $misc_set ) };
    push @arr, @arr1, @arr2;
    return \@arr;
  }
  @arr1 = @{ $mfa->fetch_all_by_Slice($sl1) };
  @arr2 = @{ $mfa->fetch_all_by_Slice($sl2) };
  push @arr, @arr1, @arr2;
  return \@arr;
  ## circular EOF

  #  if($misc_set) {
  #    return $mfa->fetch_all_by_Slice_and_set_code($self,$misc_set);
  #  }

  #  return $mfa->fetch_all_by_Slice($self);
} ## end sub get_all_MiscFeatures

=head2 get_all_AffyFeatures

  Args       : (optional) list of strings - array names
  Example    : $slice->get_all_AffyFeatures();
  Description: Retrieves all AffyFeature objects which overlap this slice.
  Returntype : listref of Bio::EnsEMBL::AffyFeature objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_AffyFeatures {
  my $self       = shift;
  my @arraynames = @_;

  my $sa = $self->adaptor();
  if ( !$sa ) {
    warning("Cannot retrieve features without attached adaptor.");
  }
  my $fa = $sa->db()->get_AffyFeatureAdaptor();
  my $features;

  ## circular BOF
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

  my ( @arr, @arr1, @arr2 );

  if (@arraynames) {
    @arr1 = @{ $fa->fetch_all_by_Slice_arrayname( $sl1, @arraynames ) };
    @arr2 = @{ $fa->fetch_all_by_Slice_arrayname( $sl2, @arraynames ) };
  } else {
    @arr1 = @{ $fa->fetch_all_by_Slice($sl1) };
    @arr2 = @{ $fa->fetch_all_by_Slice($sl2) };
  }
  push @arr, @arr1, @arr2;
  return \@arr;
  ## circular EOF

#    if ( @arraynames ) {
#        $features = $fa->fetch_all_by_Slice_arrayname( $self, @arraynames );
#    } else {
#        $features = $fa->fetch_all_by_Slice( $self );
#    }
#    return $features;
} ## end sub get_all_AffyFeatures

=head2 get_all_OligoFeatures

  Args       : (optional) list of strings - array names
  Example    : $slice->get_all_OligoFeatures();
  Description: Retrieves all OligoFeature objects which overlap this slice.
               Optionally just retrieve OligoFeature objects generated by
               probes from named arrays.
  Returntype : listref of Bio::EnsEMBL::OligoFeature objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_OligoFeatures {
  my $self       = shift;
  my @arraynames = @_;

  my $sa = $self->adaptor();
  if ( !$sa ) {
    warning("Cannot retrieve features without attached adaptor.");
  }
  my $fa = $sa->db()->get_OligoFeatureAdaptor();
  my $features;

  ## circular BOF
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

  my ( @arr, @arr1, @arr2 );

  if (@arraynames) {
    @arr1 = @{ $fa->fetch_all_by_Slice_arrayname( $sl1, @arraynames ) };
    @arr2 = @{ $fa->fetch_all_by_Slice_arrayname( $sl2, @arraynames ) };
  } else {
    @arr1 = @{ $fa->fetch_all_by_Slice($sl1) };
    @arr2 = @{ $fa->fetch_all_by_Slice($sl2) };
  }
  push @arr, @arr1, @arr2;
  return \@arr;
  ## circular EOF

#    if ( @arraynames ) {
#        $features = $fa->fetch_all_by_Slice_arrayname( $self, @arraynames );
#    } else {
#        $features = $fa->fetch_all_by_Slice( $self );
#    }
#    return $features;
} ## end sub get_all_OligoFeatures

=head2 get_all_OligoFeatures_by_type

  Arg [1]    : string - type of array (e.g. AFFY or OLIGO)
  Arg [2]    : (optional) string - logic name
  Example    : $slice->get_all_OligoFeatures_by_type('OLIGO');
  Description: Retrieves all OligoFeature objects which overlap this slice and
               were created by probes from the specified type of array.
  Returntype : listref of Bio::EnsEMBL::OligoFeature objects
  Exceptions : throws if no type
  Caller     : general
  Status     : Stable

=cut

sub get_all_OligoFeatures_by_type {
  my ( $self, $type, $logic_name ) = @_;

  throw('Need type as parameter') if !$type;

  my $sa = $self->adaptor();
  if ( !$sa ) {
    warning("Cannot retrieve features without attached adaptor.");
  }
  my $fa = $sa->db()->get_OligoFeatureAdaptor();

  ## circular BOF
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

  my ( @arr, @arr1, @arr2 );
  @arr1 = @{ $fa->fetch_all_by_Slice_type( $sl1, $type, $logic_name ) };
  @arr2 = @{ $fa->fetch_all_by_Slice_type( $sl2, $type, $logic_name ) };
  push @arr, @arr1, @arr2;
  return \@arr;
  ## circular EOF

#my $features = $fa->fetch_all_by_Slice_type( $self, $type, $logic_name );
#return $features;
} ## end sub get_all_OligoFeatures_by_type

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
  my ( $self, $logic_name, $priority, $map_weight ) = @_;

  if ( !$self->adaptor() ) {
    warning('Cannot retrieve MarkerFeatures without attached adaptor.');
    return [];
  }

  my $ma = $self->adaptor->db->get_MarkerFeatureAdaptor;

  ## circular BOF
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

  my ( @arr, @arr1, @arr2 );
  @arr1 = @{
    $ma->fetch_all_by_Slice_and_priority( $sl1, $priority, $map_weight,
                                          $logic_name ) };
  @arr2 = @{
    $ma->fetch_all_by_Slice_and_priority( $sl2, $priority, $map_weight,
                                          $logic_name ) };
  push @arr, @arr1, @arr2;
  return \@arr;
  ## circular EOF

#  my $feats = $ma->fetch_all_by_Slice_and_priority($self, $priority, $map_weight, $logic_name);
#  return $feats;
} ## end sub get_all_MarkerFeatures

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
#return $dafa->fetch_all_by_Slice($self, $qy_species, $qy_assembly, $alignment_type);
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
#return $sra->fetch_all_by_MethodLinkSpeciesSet_DnaFrag($mlss, $dnafrag, $self->start, $self->end);
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

  my ( @arr, @arr1, @arr2 );
  @arr1 = @{ $haplo_adaptor->fetch_all_by_Slice( $sl1, $lite_flag ) };
  @arr2 = @{ $haplo_adaptor->fetch_all_by_Slice( $sl2, $lite_flag ) };
  push @arr, @arr1, @arr2;
  return \@arr;
  ## circular EOF

#my $haplotypes = $haplo_adaptor->fetch_all_by_Slice($self, $lite_flag);
#return $haplotypes;
} ## end sub get_all_Haplotypes

sub get_all_DASFactories {
  my $self = shift;
  return [ $self->adaptor()->db()->_each_DASFeatureFactory ];
}


sub get_all_DASFeatures_dsn {
  my ( $self, $source_type, $dsn ) = @_;

  if ( !$self->adaptor() ) {
    warning("Cannot retrieve features without attached adaptor");
    return [];
  }
  my @X =
    grep { $_->adaptor->dsn eq $dsn }
    $self->adaptor()->db()->_each_DASFeatureFactory;

  return [ $X[0]->fetch_all_Features( $self, $source_type ) ];
}

sub get_all_DAS_Features {
  my ($self) = @_;

  $self->{_das_features} ||= {};    # Cache
  $self->{_das_styles}   ||= {};    # Cache
  $self->{_das_segments} ||= {};    # Cache
  my %das_features;
  my %das_styles;
  my %das_segments;
  my $slice = $self;

  foreach my $dasfact ( @{ $self->get_all_DASFactories } ) {
    my $dsn  = $dasfact->adaptor->dsn;
    my $name = $dasfact->adaptor->name;
    #   my $type = $dasfact->adaptor->type;
    my $url = $dasfact->adaptor->url;

    my ($type) = $dasfact->adaptor->mapping;
    if ( ref $type eq 'ARRAY' ) {
      $type = shift @$type;
    }
    $type ||= $dasfact->adaptor->type;
 # Construct a cache key : SOURCE_URL/TYPE
 # Need the type to handle sources that serve multiple types of features

    my $key = join( '/', $name, $type );
    if ( $self->{_das_features}->{$key} ) {    # Use cached
      $das_features{$name} = $self->{_das_features}->{$key};
      $das_styles{$name}   = $self->{_das_styles}->{$key};
      $das_segments{$name} = $self->{_das_segments}->{$key};
    } else {                                   # Get fresh data
      my ( $featref, $styleref, $segref ) =
        $dasfact->fetch_all_Features( $slice, $type );
      $self->{_das_features}->{$key} = $featref;
      $self->{_das_styles}->{$key}   = $styleref;
      $self->{_das_segments}->{$key} = $segref;
      $das_features{$name}           = $featref;
      $das_styles{$name}             = $styleref;
      $das_segments{$name}           = $segref;
    }
  } ## end foreach my $dasfact ( @{ $self...})

  return ( \%das_features, \%das_styles, \%das_segments );
} ## end sub get_all_DAS_Features

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

sub old_get_all_DASFeatures {
  my ( $self, @args ) = @_;

  if ( !$self->adaptor() ) {
    warning("Cannot retrieve features without attached adaptor");
    return [];
  }

  my %genomic_features =
    map { ( $_->adaptor->dsn => [ $_->fetch_all_by_Slice($self) ] ) }
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

  foreach my $xfa (@xf_adaptors) {
    push @$features, @{ $xfa->fetch_all_by_Slice($sl1) };
    push @$features, @{ $xfa->fetch_all_by_Slice($sl2) };
  }
  return $features;
  ## circular EOF
} ## end sub get_all_ExternalFeatures

=head2 get_all_DitagFeatures

  Arg [1]    : (optional) string ditag type
  Arg [1]    : (optional) string logic_name
  Example    : @dna_dna_align_feats = @{$slice->get_all_DitagFeatures};
  Description: Retrieves the DitagFeatures of a specific type which overlap
               this slice with. If type is not defined, all features are
               retrieved.
  Returntype : listref of Bio::EnsEMBL::DitagFeatures
  Exceptions : warning if slice does not have attached adaptor
  Caller     : general
  Status     : Stable

=cut

sub get_all_DitagFeatures {
  my ( $self, $type, $logic_name ) = @_;

  if ( !$self->adaptor() ) {
    warning('Cannot get DitagFeatures without attached adaptor');
    return [];
  }

  my $dfa = $self->adaptor->db->get_DitagFeatureAdaptor();

  ## circular BOF
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

  my ( @arr, @arr1, @arr2 );
  @arr1 = @{ $dfa->fetch_all_by_Slice( $sl1, $type, $logic_name ) };
  @arr2 = @{ $dfa->fetch_all_by_Slice( $sl2, $type, $logic_name ) };
  push @arr, @arr1, @arr2;
  return \@arr;
  ## circular EOF
  #return $dfa->fetch_all_by_Slice($self, $type, $logic_name);
} ## end sub get_all_DitagFeatures

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

    my ( @arr1, @arr2 );
    my $features_ref;
    @arr1 = @{ $adaptor_obj->fetch_all_by_Slice($sl1) };
    @arr2 = @{ $adaptor_obj->fetch_all_by_Slice($sl2) };
    push @{$features_ref}, @arr1, @arr2;
    ## circular EOF
    #my $features_ref = $adaptor_obj->fetch_all_by_Slice($self);

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

=head2 get_all_supercontig_Slices

  DEPRECATED use get_tiling_path("NTcontig") instead

=cut

sub get_all_supercontig_Slices {
  my $self = shift;

  deprecate("Use get_tiling_path('NTcontig') instead");

  my $result = [];

  if ( $self->adaptor() ) {
    my $superctg_names =
      $self->adaptor()->list_overlapping_supercontigs($self);

    for my $name (@$superctg_names) {
      my $slice;
      $slice = $self->adaptor()->fetch_by_supercontig_name($name);
      $slice->name($name);
      push( @$result, $slice );
    }
  } else {
    warning(
        "Slice needs to be attached to a database to get supercontigs");
  }

  return $result;
}

=head2 get_Chromosome

  DEPRECATED use this instead:
  $slice_adp->fetch_by_region('chromosome',$slice->seq_region_name)

=cut

sub get_Chromosome {
  my $self = shift @_;

  deprecate(   "Use SliceAdaptor::fetch_by_region('chromosome',"
             . '$slice->seq_region_name) instead' );

  my $csa = $self->adaptor->db->get_CoordSystemAdaptor();
  my ($top_cs) = @{ $csa->fetch_all() };

  return
    $self->adaptor->fetch_by_region( $top_cs->name(),
                              $self->seq_region_name(),
                              undef, undef, undef, $top_cs->version() );
}

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

=head2 assembly_type

  DEPRECATED use version instead

=cut

sub assembly_type {
  my $self = shift;
  deprecate('Use $slice->coord_system()->version() instead.');
  return $self->coord_system->version();
}

=head2 get_tiling_path

  DEPRECATED use project instead

=cut

sub get_tiling_path {
  my $self = shift;
  deprecate('Use $slice->project("seqlevel") instead.');
  return [];
}

=head2 dbID

  Description: DEPRECATED use SliceAdaptor::get_seq_region_id instead

=cut

sub dbID {
  my $self = shift;
  deprecate('Use SliceAdaptor::get_seq_region_id instead.');
  if ( !$self->adaptor ) {
    warning('Cannot retrieve seq_region_id without attached adaptor.');
    return 0;
  }
  return $self->adaptor->get_seq_region_id($self);
}

=head2 get_all_MapFrags

  DEPRECATED use get_all_MiscFeatures instead

=cut

sub get_all_MapFrags {
  my $self = shift;
  deprecate('Use get_all_MiscFeatures instead');
  return $self->get_all_MiscFeatures(@_);
}

=head2 has_MapSet

  DEPRECATED use get_all_MiscFeatures instead

=cut

sub has_MapSet {
  my ( $self, $mapset_name ) = @_;
  deprecate('Use get_all_MiscFeatures instead');
  my $mfs = $self->get_all_MiscFeatures($mapset_name);
  return ( @$mfs > 0 );
}

1;
