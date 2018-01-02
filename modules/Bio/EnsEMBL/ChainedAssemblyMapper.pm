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

Bio::EnsEMBL::ChainedAssemblyMapper -
Handles mapping between two coordinate systems using the information
stored in the assembly table

=head1 SYNOPSIS

    $db   = Bio::EnsEMBL::DBSQL::DBAdaptor->new(...);
    $asma = $db->get_AssemblyMapperAdaptor();
    $csa  = $db->get_CoordSystemAdaptor();

    my $chr_cs = $cs_adaptor->fetch_by_name( 'chromosome', 'NCBI33' );
    my $cln_cs = $cs_adaptor->fetch_by_name('clone');

    $asm_mapper = $map_adaptor->fetch_by_CoordSystems( $cs1, $cs2 );

    # Map to contig coordinate system from chromosomal
    @cln_coords =
      $asm_mapper->map( 'X', 1_000_000, 2_000_000, 1, $chr_cs );

    # Map to chromosome coordinate system from contig
    @chr_coords =
      $asm_mapper->map( 'AL30421.1', 100, 10000, -1, $cln_cs );

    # List contig names for a region of chromsome
    @cln_ids = $asm_mapper->list_ids( '13', 1_000_000, 1, $chr_cs );

    # List chromosome names for a contig region
    @chr_ids =
      $asm_mapper->list_ids( 'AL30421.1', 1, 1000, -1, $cln_cs );

=head1 DESCRIPTION

The ChainedAssemblyMapper is an extension of the regular AssemblyMapper
that allows for mappings between coordinate systems that require
multi-step mapping.  For example if explicit mappings are defined
between the following coordinate systems,

  chromosome <-> contig
  contig     <-> clone

the ChainedAssemblyMapper would be able to perform implicit mapping
between the chromosome and clone coordinate systems.  This should be
transparent to the user of this module, and users should not even
realise that they are using a chained assembly mapper as opposed to a
normal assembly mapper.

=head1 METHODS

=cut

package Bio::EnsEMBL::ChainedAssemblyMapper;

use strict;
use warnings;
use integer; #use proper arithmetic bitshifts

use Bio::EnsEMBL::Mapper;
use Bio::EnsEMBL::Mapper::RangeRegistry;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate);
use Scalar::Util qw(weaken);
use Bio::EnsEMBL::Utils::Scalar qw( check_ref);

my $FIRST = 'first';
my $MIDDLE = 'middle';
my $LAST  = 'last';

#2^20 = approx 10^6
my $CHUNKFACTOR = 20;

# max size of the pair cache in the mappers
my $DEFAULT_MAX_PAIR_COUNT = 6000;

=head2 new

  Arg [1]    : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
  Arg [2]    : Bio::EnsEMBL::CoordSystem $src_cs
  Arg [3]    : Bio::EnsEMBL::CoordSystem $int_cs
  Arg [4]    : Bio::EnsEMBL::CoordSystem $dst_cs
  Example    : Should use AssemblyMapperAdaptor->fetch_by_CoordSystems
  Description: Creates a new AssemblyMapper
  Returntype : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
  Exceptions : thrown if wrong number of coord_systems are provided
  Caller     : AssemblyMapperAdaptor
  Status     : Stable

=cut

sub new {
  my ($caller,$adaptor,@coord_systems) = @_;

  my $class = ref($caller) || $caller;

  my $self = {};
  bless $self, $class;

  $self->adaptor($adaptor);

  if(@coord_systems != 3) {
    throw('ChainedMapper can only map between 3 coordinate systems. ' .
          scalar(@coord_systems) . ' were provided');
  }

  $adaptor->cache_seq_ids_with_mult_assemblys();

  # Set the component, intermediate and assembled coordinate systems
  $self->{'first_cs'}   = $coord_systems[0];
  $self->{'mid_cs'}   = $coord_systems[1];
  $self->{'last_cs'}   = $coord_systems[2];

  #maps between first and intermediate coord systems
  $self->{'first_mid_mapper'} = Bio::EnsEMBL::Mapper->new($FIRST, $MIDDLE);

  #maps between last and intermediate
  $self->{'last_mid_mapper'} = Bio::EnsEMBL::Mapper->new($LAST, $MIDDLE);

  #mapper that is actually used and is loaded by the mappings generated
  #by the other two mappers
  $self->{'first_last_mapper'} = Bio::EnsEMBL::Mapper->new($FIRST, $LAST,
                                                           $coord_systems[0],
                                                           $coord_systems[2]);

  #need registries to keep track of what regions are registered in source
  #and destination coordinate systems
  $self->{'first_registry'} = Bio::EnsEMBL::Mapper::RangeRegistry->new();
  $self->{'last_registry'} = Bio::EnsEMBL::Mapper::RangeRegistry->new();

  $self->{'max_pair_count'} = $DEFAULT_MAX_PAIR_COUNT;

  return $self;
}


=head2 max_pair_count

  Arg [1]    : (optional) int $max_pair_count
  Example    : $mapper->max_pair_count(100000)
  Description: Getter/Setter for the number of mapping pairs allowed in the
               internal cache. This can be used to override the default value
               (6000) to tune the performance and memory usage for certain
               scenarios. Higher value = bigger cache, more memory used
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub max_pair_count {
  my $self = shift;
  $self->{'max_pair_count'} = shift if(@_);
  return $self->{'max_pair_count'};
}




=head2 register_all

  Arg [1]    : none
  Example    : $mapper->max_pair_count(10e6);
               $mapper->register_all();
  Description: Pre-registers all assembly information in this mapper.  The
               cache size should be set to a sufficiently large value
               so that all of the information can be stored.  This method
               is useful when *a lot* of mapping will be done in regions
               which are distributed around the genome.   After registration
               the mapper will consume a lot of memory but will not have to
               perform any SQL and will be faster.
  Returntype : none
  Exceptions : none
  Caller     : specialised programs doing a lot of mapping
  Status     : Stable

=cut

sub register_all {
  my $self = shift;
  $self->adaptor->register_all_chained($self);
  return;
}




sub flush {
  my $self = shift;
  $self->{'first_registry'}->flush();
  $self->{'last_registry'}->flush();

  $self->{'first_mid_mapper'}->flush();
  $self->{'last_mid_mapper'}->flush();
  $self->{'first_last_mapper'}->flush();
}

=head2 size

  Args       : none
  Example    : $num_of_pairs = $mapper->size();
  Description: return the number of pairs currently stored.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub size {
  my $self = shift;
  return ( $self->{'first_last_mapper'}->{'pair_count'} +
           $self->{'last_mid_mapper'}->{'pair_count'} +
           $self->{'first_mid_mapper'}->{'pair_count'} );
}



=head2 map

  Arg [1]    : string $frm_seq_region
               The name of the sequence region to transform FROM
  Arg [2]    : int $frm_start
               The start of the region to transform FROM
  Arg [3]    : int $frm_end
               The end of the region to transform FROM
  Arg [4]    : int $strand
               The strand of the region to transform FROM
  Arg [5]    : Bio::EnsEMBL::CoordSystem
               The coordinate system to transform FROM
  Arg [6]    : (optional) fastmap
  Arg [7]    : (optional) Bio::Ensembl::Slice
               The slice to transform TO
  Example    : @coords = $asm_mapper->map('X', 1_000_000, 2_000_000,
                                            1, $chr_cs);
  Description: Transforms coordinates from one coordinate system
               to another.
  Returntype : List of Bio::EnsEMBL::Mapper::Coordinate and/or
               Bio::EnsEMBL::Mapper:Gap objects
  Exceptions : thrown if the specified TO coordinat system is not one
               of the coordinate systems associated with this assembly mapper
  Caller     : general
  Status     : Stable

=cut

sub map {
  throw('Incorrect number of arguments.') if(@_ < 6);

  my ($self, $frm_seq_region_name, $frm_start,
      $frm_end, $frm_strand, $frm_cs, $fastmap, $to_slice) = @_;

  my $mapper  = $self->{'first_last_mapper'};
  my $first_cs  = $self->{'first_cs'};
  my $last_cs  = $self->{'last_cs'};

  my $is_insert = ($frm_end + 1 == $frm_start);

  my $frm;
  my $registry;




  my @tmp;
  push @tmp, $frm_seq_region_name;
  my $seq_region_id = @{$self->adaptor()->seq_regions_to_ids($frm_cs, \@tmp)}[0];

  #speed critical section:
  #try to do simple pointer equality comparisons of the coord system objects
  #first since this is likely to work most of the time and is much faster
  #than a function call

  if($frm_cs == $first_cs ||
     ($frm_cs != $last_cs && $frm_cs->equals($first_cs))) {
    $frm = $FIRST;
    $registry = $self->{'first_registry'};
  } elsif($frm_cs == $last_cs || $frm_cs->equals($last_cs)) {
    $frm = $LAST;
    $registry = $self->{'last_registry'};
  } else {
    throw("Coordinate system " . $frm_cs->name . " " . $frm_cs->version .
          " is neither the first nor the last coordinate system " .
          " of this ChainedAssemblyMapper");
  }

  #the minimum area we want to register if registration is necessary is
  #about 1MB. Break requested ranges into chunks of 1MB and then register
  #this larger region if we have a registry miss.

  #use bitwise shift for fast and easy integer multiplication and division
  my ($min_start, $min_end);

  if($is_insert) {
    $min_start = (($frm_end >> $CHUNKFACTOR) << $CHUNKFACTOR);
    $min_end   = ((($frm_start >> $CHUNKFACTOR) + 1) << $CHUNKFACTOR) - 1 ;
  } else {
    $min_start = (($frm_start >> $CHUNKFACTOR) << $CHUNKFACTOR);
    $min_end   = ((($frm_end >> $CHUNKFACTOR) + 1) << $CHUNKFACTOR) - 1 ;
  }

  #get a list of ranges in the requested region that have not been registered,
  #and register them at the same

  my $ranges;

  if($is_insert) {
    $ranges = $registry->check_and_register($seq_region_id, $frm_end,
                                            $frm_start, $min_start, $min_end);
  } else {
    $ranges = $registry->check_and_register($seq_region_id, $frm_start,
                                            $frm_end, $min_start, $min_end);
  }

  if(defined($ranges)) {
    if( $self->size() > $self->{'max_pair_count'} ) {
      $self->flush();

      if($is_insert) {
        $ranges = $registry->check_and_register
          ($seq_region_id, $frm_end, $frm_start, $min_start, $min_end);
      } else {
        $ranges = $registry->check_and_register
          ($seq_region_id, $frm_start, $frm_end, $min_start, $min_end);
      }
    }
    $self->adaptor->register_chained($self,$frm,$seq_region_id,$ranges,$to_slice);
  }

  if($fastmap) {
    return $mapper->fastmap($seq_region_id, $frm_start, $frm_end,
                            $frm_strand, $frm);
  }

  my @coords = $mapper->map_coordinates($seq_region_id, $frm_start, $frm_end,
					$frm_strand, $frm);

  # decorate (org,)mapped coordinates with their corresponding region names
  map {
    check_ref($_, 'Bio::EnsEMBL::Mapper::Coordinate') && # exclude gap
      $_->name($self->adaptor->seq_ids_to_regions([$_->id])->[0])
    } @coords;

  return @coords;
}


sub fastmap {
  my $self = shift;
  return $self->map(@_,1);
}


=head2 list_ids

  Arg [1]    : string $frm_seq_region
               The name of the sequence region of interest
  Arg [2]    : int $frm_start
               The start of the region of interest
  Arg [3]    : int $frm_end
               The end of the region to transform of interest
  Arg [5]    : Bio::EnsEMBL::CoordSystem $frm_cs
               The coordinate system to obtain overlapping ids of
  Example    : foreach $id ($asm_mapper->list_ids('X',1,1000,$chr_cs)) {...}
  Description: Retrieves a list of overlapping seq_region internal identifiers
               of another coordinate system.  This is the same as the
               list_seq_regions method but uses internal identfiers rather
               than seq_region strings
  Returntype : List of ints
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub list_ids {
  throw('Incorrect number of arguments.') if(@_ != 5);
  my($self, $frm_seq_region_name, $frm_start, $frm_end, $frm_cs) = @_;

  my $is_insert = ($frm_start == $frm_end + 1);

  #the minimum area we want to register if registration is necessary is
  #about 1MB. Break requested ranges into chunks of 1MB and then register
  #this larger region if we have a registry miss.

  #use bitwise shift for fast and easy integer multiplication and division
  my ($min_start, $min_end);

  if($is_insert) {
    $min_start = (($frm_end >> $CHUNKFACTOR) << $CHUNKFACTOR);
    $min_end   = ((($frm_start >> $CHUNKFACTOR) + 1) << $CHUNKFACTOR) - 1;
  } else {
    $min_start = (($frm_start >> $CHUNKFACTOR) << $CHUNKFACTOR);
    $min_end   = ((($frm_end >> $CHUNKFACTOR) + 1) << $CHUNKFACTOR) - 1;
  }

  my @tmp;
  push @tmp, $frm_seq_region_name;
  my $seq_region_id = @{$self->adaptor()->seq_regions_to_ids($frm_cs, \@tmp)}[0];

  if($frm_cs->equals($self->{'first_cs'})) {
    my $registry = $self->{'first_registry'};

    my $ranges;


    if($is_insert) {
      $ranges = $registry->check_and_register
        ($seq_region_id, $frm_end, $frm_start, $min_start, $min_end);
    } else {
      $ranges = $registry->check_and_register
        ($seq_region_id, $frm_start, $frm_end, $min_start, $min_end);
    }

    if(defined($ranges)) {
      $self->adaptor->register_chained($self,$FIRST,$seq_region_id,$ranges);
    }

    return map {$_->to()->id()}
      $self->first_last_mapper()->list_pairs($seq_region_id, $frm_start,
				  $frm_end, $FIRST);

  } elsif($frm_cs->equals($self->{'last_cs'})) {
    my $registry = $self->{'last_registry'};

    my $ranges;
    if($is_insert) {
      $ranges = $registry->check_and_register
        ($seq_region_id, $frm_end, $frm_start, $min_start, $min_end);
    } else {
      $ranges = $registry->check_and_register
        ($seq_region_id, $frm_start, $frm_end, $min_start, $min_end);
    }

    if(defined($ranges)) {
      $self->adaptor->register_chained($self,$LAST,$seq_region_id,$ranges);
    }

    return map {$_->from()->id()}
      $self->first_last_mapper()->list_pairs($seq_region_id, $frm_start,
				  $frm_end, $LAST);
  } else {
    throw("Coordinate system " . $frm_cs->name . " " . $frm_cs->version .
          " is neither the first nor the last coordinate system " .
          " of this ChainedAssemblyMapper");
  }
}


=head2 list_seq_regions

  Arg [1]    : string $frm_seq_region
               The name of the sequence region of interest
  Arg [2]    : int $frm_start
               The start of the region of interest
  Arg [3]    : int $frm_end
               The end of the region to transform of interest
  Arg [5]    : Bio::EnsEMBL::CoordSystem $frm_cs
               The coordinate system to obtain overlapping ids of
  Example    : foreach $id ($asm_mapper->list_ids('X',1,1000,$ctg_cs)) {...}
  Description: Retrieves a list of overlapping seq_region internal identifiers
               of another coordinate system.  This is the same as the
               list_ids method but uses seq_region names rather internal ids
  Returntype : List of strings
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub list_seq_regions {
  throw('Incorrect number of arguments.') if(@_ != 5);
  my($self, $frm_seq_region, $frm_start, $frm_end, $frm_cs) = @_;

  #retrieve the seq_region names
  my @seq_regs =
    $self->list_ids($frm_seq_region,$frm_start,$frm_end,$frm_cs);

  #The seq_regions are from the 'to' coordinate system not the
  #from coordinate system we used to obtain them
  my $to_cs;
  if($frm_cs->equals($self->first_CoordSystem())) {
    $to_cs = $self->last_CoordSystem();
  } else {
    $to_cs = $self->first_CoordSystem();
  }

  #convert them to names
  return @{$self->adaptor()->seq_ids_to_regions(\@seq_regs)};
}






=head2 first_last_mapper

  Args       : none
  Example    : $mapper = $cam->first_last_mapper();
  Description: return the mapper.
  Returntype : Bio::EnsEMBL::Mapper
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub first_last_mapper {
  my $self = shift;
  return $self->{'first_last_mapper'};
}

=head2 first_middle_mapper

  Args       : none
  Example    : $mapper = $cam->first_middle_mapper();
  Description: return the mapper.
  Returntype : Bio::EnsEMBL::Mapper
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut


sub first_middle_mapper {
  my $self = shift;
  return $self->{'first_mid_mapper'};
}

=head2 last_middle_mapper

  Args       : none
  Example    : $mapper = $cam->last_middle_mapper();
  Description: return the mapper.
  Returntype : Bio::EnsEMBL::Mapper
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub last_middle_mapper {
  my $self = shift;
  return $self->{'last_mid_mapper'};
}


=head2 first_CoordSystem

  Args       : none
  Example    : $coordsys = $cam->first_CoordSystem();
  Description: return the CoordSystem.
  Returntype : Bio::EnsEMBL::CoordSystem
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub first_CoordSystem {
  my $self = shift;
  return $self->{'first_cs'};
}


=head2 middle_CoordSystem

  Args       : none
  Example    : $coordsys = $cam->middle_CoordSystem();
  Description: return the CoordSystem.
  Returntype : Bio::EnsEMBL::CoordSystem
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub middle_CoordSystem {
  my $self = shift;
  return $self->{'mid_cs'};
}

=head2 last_CoordSystem

  Args       : none
  Example    : $coordsys = $cam->last_CoordSystem();
  Description: return the CoordSystem.
  Returntype : Bio::EnsEMBL::CoordSystem
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub last_CoordSystem {
  my $self = shift;
  return $self->{'last_cs'};
}

=head2 first_registry

  Args       : none
  Example    : $rr = $cam->first_registry();
  Description: return the Registry.
  Returntype : Bio::EnsEMBL::Mapper::RangeRegistry
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub first_registry {
  my $self = shift;
  return $self->{'first_registry'};
}

=head2 last_registry

  Args       : none
  Example    : $rr = $cam->last_registry();
  Description: return the Registry.
  Returntype : Bio::EnsEMBL::Mapper::RangeRegistry
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub last_registry {
  my $self = shift;
  return $self->{'last_registry'};
}


#
# Methods supplied to maintain polymorphism with AssemblyMapper there
# is no real assembled or component in the chained mapper, since the
# ordering is arbitrary and both ends might actually be assembled, but
# these methods provide convenient synonyms
#

=head2 mapper

  Args       : none
  Example    : $mapper = $cam->mapper();
  Description: return the first_last_mapper.
  Returntype : Bio::EnsEMBL::Mapper
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub mapper {
  my $self = shift;
  return $self->first_last_mapper();
}

=head2 assembled_CoordSystem

  Args       : none
  Example    : $coordsys = $cam->assembled_CoordSystem();
  Description: return the first CoordSystem.
  Returntype : Bio::EnsEMBL::CoordSystem
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut


sub assembled_CoordSystem {
  my $self = shift;
  return $self->{'first_cs'};
}

=head2 component_CoordSystem

  Args       : none
  Example    : $coordsys = $cam->component_CoordSystem();
  Description: return the last CoordSystem.
  Returntype : Bio::EnsEMBL::CoordSystem
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub component_CoordSystem {
  my $self = shift;
  return $self->{'last_cs'};
}


=head2 adaptor

  Arg [1]    : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor $adaptor
  Description: get/set for this objects database adaptor
  Returntype : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub adaptor {
  my $self = shift;
  weaken($self->{'adaptor'} = shift) if(@_);
  return $self->{'adaptor'};
}


1;
