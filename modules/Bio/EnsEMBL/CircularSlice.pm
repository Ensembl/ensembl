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
        $seqAdaptor->fetch_by_Slice_start_end_strand( $sl1, 1, $sl1->end - $sl1->start + 1, $sl1->strand) };
                                                      
      my $seq2 = ${
        $seqAdaptor->fetch_by_Slice_start_end_strand( $sl2, 1, $sl2->end - $sl2->start + 1, $sl2->strand) };
                                                      
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


=head2 get_all_DASFeatures

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


1;
