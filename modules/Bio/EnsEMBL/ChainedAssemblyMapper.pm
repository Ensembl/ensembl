
#
# Ensembl module for Bio::EnsEMBL::ChainedAssemblyMapper
#
# Written by Graham McVicker
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::ChainedAssemblyMapper - 
Handles mapping between two coordinate systems using the information stored in
the assembly table

=head1 SYNOPSIS
    $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(...);
    $asma = $db->get_AssemblyMapperAdaptor();
    $csa  = $db->get_CoordSystemAdaptor();

    my $chr_cs = $cs_adaptor->fetch_by_name('chromosome', 'NCBI33');
    my $cln_cs   = $cs_adaptor->fetch_by_name('clone');

    $asm_mapper = $map_adaptor->fetch_by_CoordSystems($cs1, $cs2);

    #map to contig coordinate system from chromosomal
    @cln_coords = $asm_mapper->map('X', 1_000_000, 2_000_000, 1, $chr_cs);

    #map to chromosome coordinate system from contig
    @chr_coords = $asm_mapper->map('AL30421.1',100,10000,-1,$cln_cs);

    #list contig names for a region of chromsome
    @cln_ids = $asm_mapper->list_ids('13', 1_000_000, 1, $chr_cs);

    #list chromosome names for a contig region
    @chr_ids = $asm_mapper->list_ids('AL30421.1',1,1000,-1,$cln_cs);

=head1 DESCRIPTION

The ChainedAssemblyMapper is an extension of the regular AssemblyMapper that
allows for mappings between coordinate systems that require multi-step mapping.
For example if explicit mappings are defined between the following 
coordinate systems,
  chromosome <-> contig
  contig     <-> clone
the ChainedAssemblyMapper would be able to perform implicit mapping between
the chromosome and clone coordinate systems.  This should be transparent to
the user of this module, and users should not even realise that they are using
a chained assembly mapper as opposed to a normal assembly mapper.

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 METHODS

=cut


my $FIRST = 'first';
my $MIDDLE = 'middle';
my $LAST  = 'last';

package Bio::EnsEMBL::ChainedAssemblyMapper;

use strict;
use warnings;
use integer; #use proper arithmetic bitshifts

use Bio::EnsEMBL::Mapper;
use Bio::EnsEMBL::Mapper::RangeRegistry;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate);

my $CHUNKFACTOR = 20; #2^20 = approx 10^6
my $MAX_PAIR_COUNT = 6000; # max size of the pair cache in the mappers

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

  return $self;
}


sub flush {
  my $self = shift;
  $self->{'first_registry'}->flush();
  $self->{'last_registry'}->flush();

  $self->{'first_mid_mapper'}->flush();
  $self->{'last_mid_mapper'}->flush();
  $self->{'first_last_mapper'}->flush();
}

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
  Example    : @coords = $asm_mapper->map('X', 1_000_000, 2_000_000,
                                            1, $chr_cs);
  Description: Transforms coordinates from one coordinate system
               to another.
  Returntype : List of Bio::EnsEMBL::Mapper::Coordinate and/or
               Bio::EnsEMBL::Mapper:Gap objects
  Exceptions : thrown if the specified TO coordinat system is not one
               of the coordinate systems associated with this assembly mapper
  Caller     : general

=cut

sub map {
  throw('Incorrect number of arguments.') if(@_ < 6);

  my ($self, $frm_seq_region, $frm_start,
      $frm_end, $frm_strand, $frm_cs, $fastmap) = @_;

  my $mapper  = $self->{'first_last_mapper'};
  my $first_cs  = $self->{'first_cs'};
  my $last_cs  = $self->{'last_cs'};

  my $frm;
  my $registry;

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

  $min_start = (($frm_start >> $CHUNKFACTOR) << $CHUNKFACTOR);
  $min_end   = ((($frm_end >> $CHUNKFACTOR) + 1) << $CHUNKFACTOR) - 1 ;

  #get a list of ranges in the requested region that have not been registered,
  #and register them at the same

  #print STDERR "frm_start=$frm_start frm_end=$frm_end" .
  #              "min_start=$min_start min_end=$min_end\n";

  my $ranges =
    $registry->check_and_register($frm_seq_region, $frm_start, $frm_end,
				  $min_start, $min_end);

  if(defined($ranges)) {
    if( $self->size() > $MAX_PAIR_COUNT ) {
      $self->flush();
      $ranges =
	$registry->check_and_register($frm_seq_region, $frm_start, $frm_end,
				      $min_start, $min_end);
    }
    $self->adaptor->register_chained($self,$frm,$frm_seq_region,$ranges);
  }

  if($fastmap) {
    return $mapper->fastmap($frm_seq_region, $frm_start, $frm_end,
			    $frm_strand, $frm);
  }

  return $mapper->map_coordinates($frm_seq_region, $frm_start, $frm_end,
                                  $frm_strand, $frm);
}


sub fastmap {
  my $self = shift;
  return $self->map(@_,1);
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

=cut


sub list_seq_regions {
  throw('Incorrect number of arguments.') if(@_ != 5);
  my($self, $frm_seq_region, $frm_start, $frm_end, $frm_cs) = @_;

  #the minimum area we want to register if registration is necessary is
  #about 1MB. Break requested ranges into chunks of 1MB and then register
  #this larger region if we have a registry miss.

  #use bitwise shift for fast and easy integer multiplication and division
  my $min_start;
  my $min_end;

  $min_start = (($frm_start >> $CHUNKFACTOR) << $CHUNKFACTOR);
  $min_end   = ((($frm_end >> $CHUNKFACTOR) + 1) << $CHUNKFACTOR) - 1;

  if($frm_cs->equals($self->{'first_cs'})) {
    my $registry = $self->{'first_registry'};
    my $ranges =
      $registry->check_and_register($frm_seq_region, $frm_start, $frm_end,
				   $min_start, $min_end);

    if(defined($ranges)) {
      $self->adaptor->register_chained($self,$FIRST,$frm_seq_region,$ranges);
    }

    return map {$_->to()->id()}
      $self->first_last_mapper()->list_pairs($frm_seq_region, $frm_start,
				  $frm_end, $FIRST);

  } elsif($frm_cs->equals($self->{'last_cs'})) {
    my $registry = $self->{'last_registry'};
    my $ranges =
      $registry->check_and_register($frm_seq_region, $frm_start, $frm_end,
				    $min_start, $min_end);

    if(defined($ranges)) {
      $self->adaptor->register_chained($self,$LAST,$frm_seq_region,$ranges);
    }

    return map {$_->from()->id()}
      $self->first_last_mapper()->list_pairs($frm_seq_region, $frm_start,
				  $frm_end, $LAST);
  } else {
    throw("Coordinate system " . $frm_cs->name . " " . $frm_cs->version .
          " is neither the first nor the last coordinate system " .
          " of this ChainedAssemblyMapper");
  }
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

=cut

sub list_ids {
  throw('Incorrect number of arguments.') if(@_ != 5);
  my($self, $frm_seq_region, $frm_start, $frm_end, $frm_cs) = @_;

  #retrieve the seq_region names
  my @seq_regs =
    $self->list_seq_regions($frm_seq_region,$frm_start,$frm_end,$frm_cs);

  #The seq_regions are from the 'to' coordinate system not the
  #from coordinate system we used to obtain them
  my $to_cs;
  if($frm_cs->equals($self->first_CoordSystem())) {
    $to_cs = $self->last_CoordSystem();
  } else {
    $to_cs = $self->first_CoordSystem();
  }

  #convert them to ids
  return @{$self->adaptor()->seq_regions_to_ids($to_cs, \@seq_regs)};
}






sub first_last_mapper {
  my $self = shift;
  return $self->{'first_last_mapper'};
}

sub first_middle_mapper {
  my $self = shift;
  return $self->{'first_mid_mapper'};
}

sub last_middle_mapper {
  my $self = shift;
  return $self->{'last_mid_mapper'};
}


sub first_CoordSystem {
  my $self = shift;
  return $self->{'first_cs'};
}


sub middle_CoordSystem {
  my $self = shift;
  return $self->{'mid_cs'};
}

sub last_CoordSystem {
  my $self = shift;
  return $self->{'last_cs'};
}

sub first_registry {
  my $self = shift;
  return $self->{'first_registry'};
}

sub last_registry {
  my $self = shift;
  return $self->{'last_registry'};
}





#
# Methods supplied to maintain polymorphism with AssemblyMapper
# there is no real assembled or component in the chained mapper, since the
# ordering is arbitrary and both ends might actually be assembled, but these
# methods provide convenient synonyms
#
sub mapper {
  my $self = shift;
  return $self->first_last_mapper();
}
sub assembled_CoordSystem {
  my $self = shift;
  return $self->{'first_cs'};
}
sub component_CoordSystem {
  my $self = shift;
  return $self->{'last_cs'};
}


=head2 adaptor

  Arg [1]    : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor $adaptor
  Example    : none
  Description: get/set for this objects database adaptor
  Returntype : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
  Exceptions : none
  Caller     : general

=cut

sub adaptor {
  my $self = shift;
  $self->{'adaptor'} = shift if(@_);
  return $self->{'adaptor'};
}


=head2 in_assembly

  Description: Deprecated. Use map() or list_ids() instead

=cut

sub in_assembly {
  my ($self, $object) = @_;

  deprecate('Use map() or list_ids() instead.');

  my $csa = $self->db->get_CoordSystemAdaptor();

  my $top_level = $csa->fetch_top_level();

  my $asma = $self->adaptor->fetch_by_CoordSystems($object->coord_system(),
                                                   $top_level);

  my @list = $asma->list_ids($object->seq_region(), $object->start(),
                             $object->end(), $object->coord_system());

  return (@list > 0);
}


=head2 map_coordinates_to_assembly

  Description: DEPRECATED use map() instead

=cut

sub map_coordinates_to_assembly {
  my ($self, $contig_id, $start, $end, $strand) = @_;

  deprecate('Use map() instead.');

  #not sure if contig_id is seq_region_id or name...
  return $self->map($contig_id, $start, $end, $strand,
                   $self->contig_CoordSystem());

}


=head2 fast_to_assembly

  Description: DEPRECATED use map() instead

=cut

sub fast_to_assembly {
  my ($self, $contig_id, $start, $end, $strand) = @_;

  deprecate('Use map() instead.');

  #not sure if contig_id is seq_region_id or name...
  return $self->map($contig_id, $start, $end, $strand,
                    $self->contig_CoordSystem());
}


=head2 map_coordinates_to_rawcontig

  Description: DEPRECATED use map() instead

=cut

sub map_coordinates_to_rawcontig {
  my ($self, $chr_name, $start, $end, $strand) = @_;

  deprecate('Use map() instead.');

  return $self->map($chr_name, $start, $end, $strand,
                    $self->assembled_CoordSystem());

}

=head2 list_contig_ids
  Description: DEPRECATED Use list_ids instead

=cut

sub list_contig_ids {
  my ($self, $chr_name, $start, $end) = @_;

  deprecate('Use list_ids() instead.');

  return $self->list_ids($chr_name, $start, $end, 
                         $self->assembled_CoordSystem());
}



1;
