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

Bio::EnsEMBL::AssemblyMapper - 
Handles mapping between two coordinate systems using the information
stored in the assembly table.

=head1 SYNOPSIS

    $db   = Bio::EnsEMBL::DBSQL::DBAdaptor->new(...);
    $asma = $db->get_AssemblyMapperAdaptor();
    $csa  = $db->get_CoordSystemAdaptor();

    my $chr_cs = $cs_adaptor->fetch_by_name( 'chromosome', 'NCBI33' );
    my $ctg_cs = $cs_adaptor->fetch_by_name('contig');

    $asm_mapper = $map_adaptor->fetch_by_CoordSystems( $cs1, $cs2 );

    # Map to contig coordinate system from chromosomal.
    @ctg_coords =
      $asm_mapper->map( 'X', 1_000_000, 2_000_000, 1, $chr_cs );

    # Map to chromosome coordinate system from contig.
    @chr_coords =
      $asm_mapper->map( 'AL30421.1.200.92341', 100, 10000, -1,
      $ctg_cs );

    # List contig names for a region of chromsome.
    @ctg_ids = $asm_mapper->list_ids( '13', 1_000_000, 1, $chr_cs );

    # List chromosome names for a contig region.
    @chr_ids =
      $asm_mapper->list_ids( 'AL30421.1.200.92341', 1, 1000, -1,
      $ctg_cs );

=head1 DESCRIPTION

The AssemblyMapper is a database aware mapper which faciliates
conversion of coordinates between any two coordinate systems with an
relationship explicitly defined in the assembly table.  In the future
it may be possible to perform multiple step (implicit) mapping between
coordinate systems.

It is implemented using the Bio::EnsEMBL::Mapper object, which is a
generic mapper object between disjoint coordinate systems.

=head1 METHODS

=cut


package Bio::EnsEMBL::AssemblyMapper;

use strict;
use warnings;

use Bio::EnsEMBL::Mapper;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate);
use Scalar::Util qw(weaken);
use Bio::EnsEMBL::Utils::Scalar qw( check_ref);

my $ASSEMBLED = 'assembled';
my $COMPONENT = 'component';

my $DEFAULT_MAX_PAIR_COUNT = 1000;


=head2 new

  Arg [1]    : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
  Arg [2]    : Bio::EnsEMBL::CoordSystem $asm_cs
  Arg [3]    : Bio::EnsEMBL::CoordSystem $cmp_cs
  Example    : Should use AssemblyMapperAdaptor->fetch_by_CoordSystems()
  Description: Creates a new AssemblyMapper
  Returntype : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
  Exceptions : Throws if multiple coord_systems are provided
  Caller     : AssemblyMapperAdaptor
  Status     : Stable

=cut

sub new {
  my ( $proto, $adaptor, @coord_systems ) = @_;

  my $class = ref($proto) || $proto;

  my $self = bless( {}, $class );

  $self->adaptor($adaptor);

  $adaptor->cache_seq_ids_with_mult_assemblys();

  if ( @coord_systems != 2 ) {
    throw(   'Can only map between two coordinate systems. '
           . scalar(@coord_systems)
           . ' were provided' );
  }

  # Set the component and assembled coordinate systems
  $self->{'asm_cs'} = $coord_systems[0];
  $self->{'cmp_cs'} = $coord_systems[1];

  # We load the mapper calling the 'ASSEMBLED' the 'from' coord system
  # and the 'COMPONENT' the 'to' coord system.

  $self->{'mapper'} = Bio::EnsEMBL::Mapper->new( $ASSEMBLED, $COMPONENT,
                                 $coord_systems[0], $coord_systems[1] );

  $self->{'max_pair_count'} = $DEFAULT_MAX_PAIR_COUNT;

  return $self;
} ## end sub new

=head2 max_pair_count

  Arg [1]    : (optional) int $max_pair_count
  Example    : $mapper->max_pair_count(100000)
  Description: Getter/Setter for the number of mapping pairs allowed
               in the internal cache.  This can be used to override
               the default value (1000) to tune the performance and
               memory usage for certain scenarios.  Higher value
               means bigger cache, more memory used.
  Return type: int
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub max_pair_count {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'max_pair_count'} = $value;
  }

  return $self->{'max_pair_count'};
}

=head2 register_all

  Arg [1]    : None
  Example    : $mapper->max_pair_count(10e6);
               $mapper->register_all();
  Description: Pre-registers all assembly information in this
               mapper.  The cache size should be set to a
               sufficiently large value so that all of the
               information can be stored.  This method is useful
               when *a lot* of mapping will be done in regions
               which are distributed around the genome.  After
               registration the mapper will consume a lot of memory
               but will not have to perform any SQL and will be
               faster.
  Return type: None
  Exceptions : None
  Caller     : Specialised programs doing a lot of mapping.
  Status     : Stable

=cut

sub register_all {
  my ($self) = @_;

  $self->adaptor()->register_all($self);
}

=head2 map

  Arg [1]    : string $frm_seq_region
               The name of the sequence region to transform FROM.
  Arg [2]    : int $frm_start
               The start of the region to transform FROM.
  Arg [3]    : int $frm_end
               The end of the region to transform FROM.
  Arg [4]    : int $strand
               The strand of the region to transform FROM.
  Arg [5]    : Bio::EnsEMBL::CoordSystem
               The coordinate system to transform FROM
  Example    : @coords =
                $asm_mapper->map( 'X', 1_000_000, 2_000_000, 1,
                                  $chr_cs );
  Description: Transforms coordinates from one coordinate system to
               another.
  Return type: List of Bio::EnsEMBL::Mapper::Coordinate and/or
               Bio::EnsEMBL::Mapper:Gap objects.
  Exceptions : Throws if if the specified TO coordinat system is not
               one of the coordinate systems associated with this
               assembly mapper.
  Caller     : General
  Status     : Stable

=cut

sub map {
  throw('Incorrect number of arguments.') if (!( @_ >= 6));

  my ( $self, $frm_seq_region_name, $frm_start, $frm_end, $frm_strand,
       $frm_cs, $to_slice )
    = @_;

  my $mapper  = $self->{'mapper'};
  my $asm_cs  = $self->{'asm_cs'};
  my $cmp_cs  = $self->{'cmp_cs'};
  my $adaptor = $self->{'adaptor'};
  my $frm;


  my $seq_region_id =
    $self->adaptor()
    ->seq_regions_to_ids( $frm_cs, [$frm_seq_region_name] )->[0];

  # Speed critical section:
  # Try to do simple pointer equality comparisons of the coord system
  # objects first since this is likely to work most of the time and is
  # much faster than a function call.

  if ( $frm_cs == $cmp_cs
       || ( $frm_cs != $asm_cs && $frm_cs->equals($cmp_cs) ) )
  {
    if ( !$self->{'cmp_register'}->{$seq_region_id} ) {
      $adaptor->register_component( $self, $seq_region_id );
    }
    $frm = $COMPONENT;

  } elsif ( $frm_cs == $asm_cs || $frm_cs->equals($asm_cs) ) {

    # This can be probably be sped up some by only calling registered
    # assembled if needed.
    $adaptor->register_assembled( $self, $seq_region_id, $frm_start,
                                  $frm_end );
    $frm = $ASSEMBLED;

  } else {

    throw(
           sprintf( "Coordinate system %s %s is neither the assembled "
                      . "nor the component coordinate system "
                      . "of this AssemblyMapper",
                    $frm_cs->name(), $frm_cs->version() ) );

  }

  my @coords = 
    $mapper->map_coordinates( $seq_region_id, $frm_start, $frm_end,
                              $frm_strand, $frm );
  
  # decorate (org,)mapped coordinates with their corresponding region names
  map {
    check_ref($_, 'Bio::EnsEMBL::Mapper::Coordinate') && # exclude gap
      $_->name($adaptor->seq_ids_to_regions([$_->id])->[0])
    } @coords;

  return @coords;
  
} ## end sub map


=head2 flush

  Args       : None
  Example    : None
  Description: Remove all cached items from this AssemblyMapper.
  Return type: None
  Exceptions : None
  Caller     : AssemblyMapperAdaptor
  Status     : Stable

=cut

sub flush {
  my ($self) = @_;

  $self->{'mapper'}->flush();
  $self->{'cmp_register'} = {};
  $self->{'asm_register'} = {};
}

=head2 size

  Args       : None
  Example    : $num_of_pairs = $mapper->size();
  Description: Returns the number of pairs currently stored.
  Return type: int
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub size {
  my ($self) = @_;

  return $self->{'mapper'}->{'pair_count'};
}

=head2 fastmap

  Arg [1]    : string $frm_seq_region
               The name of the sequence region to transform FROM.
  Arg [2]    : int $frm_start
               The start of the region to transform FROM.
  Arg [3]    : int $frm_end
               The end of the region to transform FROM.
  Arg [4]    : int $strand
               The strand of the region to transform FROM.
  Arg [5]    : Bio::EnsEMBL::CoordSystem
               The coordinate system to transform FROM.
  Example    : @coords =
                $asm_mapper->map( 'X', 1_000_000, 2_000_000, 1,
                                  $chr_cs );
  Description: Transforms coordinates from one coordinate system to
               another.
  Return type: List of Bio::EnsEMBL::Mapper::Coordinate and/or
               Bio::EnsEMBL::Mapper:Gap objects.
  Exceptions : Throws if the specified TO coordinat system is not
               one of the coordinate systems associated with this
               assembly mapper.
  Caller     : General
  Status     : Stable

=cut

sub fastmap {
  if ( @_ != 6 ) {
    throw('Incorrect number of arguments.');
  }

  my ( $self, $frm_seq_region_name, $frm_start, $frm_end, $frm_strand,
       $frm_cs )
    = @_;

  my $mapper  = $self->{'mapper'};
  my $asm_cs  = $self->{'asm_cs'};
  my $cmp_cs  = $self->{'cmp_cs'};
  my $adaptor = $self->adaptor();
  my $frm;

  my @tmp;
  push @tmp, $frm_seq_region_name;

  my $seq_region_id =
    $self->adaptor()->seq_regions_to_ids( $frm_cs, \@tmp )->[0];

  # Speed critical section:
  # Try to do simple pointer equality comparisons of the coord system
  # objects first since this is likely to work most of the time and is
  # much faster than a function call.

  if ( $frm_cs == $cmp_cs
       || ( $frm_cs != $asm_cs && $frm_cs->equals($cmp_cs) ) )
  {

    if ( !$self->{'cmp_register'}->{$seq_region_id} ) {
      $adaptor->register_component( $self, $seq_region_id );
    }
    $frm = $COMPONENT;

  } elsif ( $frm_cs == $asm_cs || $frm_cs->equals($asm_cs) ) {

    # This can be probably be sped up some by only calling registered
    # assembled if needed
    $adaptor->register_assembled( $self, $seq_region_id, $frm_start,
                                  $frm_end );
    $frm = $ASSEMBLED;

  } else {

    throw(
           sprintf( "Coordinate system %s %s is neither the assembled "
                      . "nor the component coordinate system "
                      . "of this AssemblyMapper",
                    $frm_cs->name(), $frm_cs->version() ) );

  }

  return
    $mapper->fastmap( $seq_region_id, $frm_start, $frm_end, $frm_strand,
                      $frm );
} ## end sub fastmap

=head2 list_ids

  Arg [1]    : string $frm_seq_region
               The name of the sequence region of interest.
  Arg [2]    : int $frm_start
               The start of the region of interest.
  Arg [3]    : int $frm_end
               The end of the region to transform of interest.
  Arg [5]    : Bio::EnsEMBL::CoordSystem $frm_cs
               The coordinate system to obtain overlapping IDs of.
  Example    : foreach my $id (
                        $asm_mapper->list_ids( 'X', 1, 1000, $ctg_cs ) )
                { ... }
  Description: Retrieves a list of overlapping seq_region names of
               another coordinate system.  This is the same as the
               list_ids method but uses seq_region names rather
               internal IDs.
  Return type: List of strings.
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub list_ids {
  if ( @_ != 5 ) {
    throw('Incorrect number of arguments.');
  }

  my ( $self, $frm_seq_region_name, $frm_start, $frm_end, $frm_cs ) =
    @_;

  my @tmp = ($frm_seq_region_name);

  my $seq_region_id =
    $self->adaptor()->seq_regions_to_ids( $frm_cs, \@tmp )->[0];

  if ( $frm_cs->equals( $self->component_CoordSystem() ) ) {

    if ( !$self->have_registered_component($seq_region_id) ) {
      $self->adaptor->register_component( $self, $seq_region_id );
    }

    # Pull out the 'from' identifiers of the mapper pairs.  The we
    # loaded the assembled side as the 'from' side in the constructor.

    return
      map ( { $_->from()->id() }
            $self->mapper()->list_pairs(
                        $seq_region_id, $frm_start, $frm_end, $COMPONENT
            ) );

  } elsif ( $frm_cs->equals( $self->assembled_CoordSystem() ) ) {

    $self->adaptor->register_assembled( $self, $seq_region_id,
                                        $frm_start, $frm_end );

    # Pull out the 'to' identifiers of the mapper pairs we loaded the
    # component side as the 'to' coord system in the constructor.

    return
      map ( { $_->to->id() }
            $self->mapper()->list_pairs(
                        $seq_region_id, $frm_start, $frm_end, $ASSEMBLED
            ) );

  } else {

    throw(
           sprintf( "Coordinate system %s %s is neither the assembled "
                      . "nor the component coordinate system "
                      . "of this AssemblyMapper",
                    $frm_cs->name(), $frm_cs->version() ) );

  }
} ## end sub list_ids

=head2 list_seq_regions

  Arg [1]    : string $frm_seq_region
               The name of the sequence region of interest.
  Arg [2]    : int $frm_start
               The start of the region of interest.
  Arg [3]    : int $frm_end
               The end of the region to transform of interest.
  Arg [5]    : Bio::EnsEMBL::CoordSystem $frm_cs
               The coordinate system to obtain overlapping IDs of.
  Example    : foreach my $id (
                                 $asm_mapper->list_seq_regions(
                                                   'X', 1, 1000, $chr_cs
                                 ) ) { ... }
  Description: Retrieves a list of overlapping seq_region internal
               identifiers of another coordinate system.  This is
               the same as the list_seq_regions method but uses
               internal identfiers rather than seq_region strings.
  Return type: List of ints.
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub list_seq_regions {
  if ( @_ != 5 ) {
    throw('Incorrect number of arguments.');
  }

  my ( $self, $frm_seq_region, $frm_start, $frm_end, $frm_cs ) = @_;

  # Retrieve the seq_region names.

  my @seq_ids =
    $self->list_ids( $frm_seq_region, $frm_start, $frm_end, $frm_cs );

  # The seq_regions are from the 'to' coordinate system not the from
  # coordinate system we used to obtain them.

  my $to_cs;
  if ( $frm_cs->equals( $self->assembled_CoordSystem() ) ) {
    $to_cs = $self->component_CoordSystem();
  } else {
    $to_cs = $self->assembled_CoordSystem();
  }

  # Convert them to IDs.
  return @{ $self->adaptor()->seq_ids_to_regions( \@seq_ids ) };
}


=head2 have_registered_component

  Arg [1]    : string $cmp_seq_region
               The name of the sequence region to check for
               registration.
  Example    : if ( $asm_mapper->have_registered_component('AL240214.1') ) {}
  Description: Returns true if a given component region has
               been registered with this assembly mapper.  This
               should only be called by this class or the
               AssemblyMapperAdaptor.  In other words, do not use
               this method unless you really know what you are
               doing.
  Return type: Boolean (0 or 1)
  Exceptions : Throws on incorrect arguments.
  Caller     : Internal, AssemblyMapperAdaptor
  Status     : Stable

=cut

sub have_registered_component {
  my ( $self, $cmp_seq_region ) = @_;

  if ( !defined($cmp_seq_region) ) {
    throw('cmp_seq_region argument is required');
  }

  if ( exists( $self->{'cmp_register'}->{$cmp_seq_region} ) ) {
    return 1;
  }

  return 0;
}

=head2 have_registered_assembled

  Arg [1]    : string $asm_seq_region
               The name of the sequence region to check for
               registration.
  Arg [2]    : int $chunk_id
               The chunk number of the provided seq_region to check
               for registration.
  Example    : if ( $asm_mapper->have_registered_component( 'X', 9 ) ) { }
  Description: Returns true if a given assembled region chunk
               has been registered with this assembly mapper.
               This should only be called by this class or the
               AssemblyMapperAdaptor.  In other words, do not use
               this method unless you really know what you are
               doing.
  Return type: Boolean (0 or 1)
  Exceptions : Throws on incorrect arguments
  Caller     : Internal, AssemblyMapperAdaptor
  Status     : Stable

=cut

sub have_registered_assembled {
  my ( $self, $asm_seq_region, $chunk_id ) = @_;

  if ( !defined($asm_seq_region) ) {
    throw('asm_seq_region argument is required');
  }
  if ( !defined($chunk_id) ) {
    throw('chunk_id is required');
  }

  if (
     exists( $self->{'asm_register'}->{$asm_seq_region}->{$chunk_id} ) )
  {
    return 1;
  }

  return 0;
}


=head2 register_component

  Arg [1]    : integer $cmp_seq_region
               The dbID of the component sequence region to
               register.
  Example    : $asm_mapper->register_component('AL312341.1');
  Description: Flags a given component sequence region as registered
               in this assembly mapper.  This should only be called
               by this class or the AssemblyMapperAdaptor.
  Return type: None
  Exceptions : Throws on incorrect arguments
  Caller     : Internal, AssemblyMapperAdaptor
  Status     : Stable

=cut

sub register_component {
  my ( $self, $cmp_seq_region ) = @_;

  if ( !defined($cmp_seq_region) ) {
    throw('cmp_seq_region argument is required');
  }

  $self->{'cmp_register'}->{$cmp_seq_region} = 1;
}

=head2 register_assembled

  Arg [1]    : integer $asm_seq_region
               The dbID of the sequence region to register.
  Arg [2]    : int $chunk_id
               The chunk number of the provided seq_region to register.
  Example    : $asm_mapper->register_assembled( 'X', 4 );
  Description: Flags a given assembled region as registered in this
               assembly mapper.  This should only be called by this
               class or the AssemblyMapperAdaptor.  Do not call this
               method unless you really know what you are doing.
  Return type: None
  Exceptions : Throws on incorrect arguments
  Caller     : Internal, AssemblyMapperAdaptor
  Status     : Stable

=cut

sub register_assembled {
  my ( $self, $asm_seq_region, $chunk_id ) = @_;

  if ( !defined($asm_seq_region) ) {
    throw('asm_seq_region argument is required');
  }
  if ( !defined($chunk_id) ) {
    throw('chunk_id srgument is required');
  }

  $self->{'asm_register'}->{$asm_seq_region}->{$chunk_id} = 1;
}

=head2 mapper

  Arg [1]    : None
  Example    : $mapper = $asm_mapper->mapper();
  Description: Retrieves the internal mapper used by this Assembly
               Mapper.  This is unlikely to be useful unless you
               _really_ know what you are doing.
  Return type: Bio::EnsEMBL::Mapper
  Exceptions : None
  Caller     : Internal, AssemblyMapperAdaptor
  Status     : Stable

=cut

sub mapper {
  my ($self) = @_;

  return $self->{'mapper'};
}

=head2 assembled_CoordSystem

  Arg [1]    : None
  Example    : $cs = $asm_mapper->assembled_CoordSystem();
  Description: Retrieves the assembled CoordSystem from this
               assembly mapper.
  Return type: Bio::EnsEMBL::CoordSystem
  Exceptions : None
  Caller     : Internal, AssemblyMapperAdaptor
  Status     : Stable

=cut

sub assembled_CoordSystem {
  my ($self) = @_;

  return $self->{'asm_cs'};
}

=head2 component_CoordSystem

  Arg [1]    : None
  Example    : $cs = $asm_mapper->component_CoordSystem();
  Description: Retrieves the component CoordSystem from this
               assembly mapper.
  Return type: Bio::EnsEMBL::CoordSystem
  Exceptions : None
  Caller     : Internal, AssemblyMapperAdaptor
  Status     : Stable

=cut

sub component_CoordSystem {
  my ($self) = @_;

  return $self->{'cmp_cs'};
}

=head2 adaptor

  Arg [1]    : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor $adaptor
  Description: Getter/set terfor this object's database adaptor.
  Returntype : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub adaptor {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    weaken($self->{'adaptor'} = $value);
  }

  return $self->{'adaptor'};
}

1;
