
#
# Ensembl module for Bio::EnsEMBL::AssemblyMapper
#
# Written by Arne Stabenau <stabenau@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::AssemblyMapper - 
Handles mapping between two coordinate systems using the information stored in
the assembly table

=head1 SYNOPSIS
    $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(...);
    $asma = $db->get_AssemblyMapperAdaptor();
    $csa  = $db->get_CoordSystemAdaptor();

    my $clone_cs = $cs_adaptor->fetch_by_name('chromosome', 'NCBI33');
    my $ctg_cs   = $cs_adaptor->fetch_by_name('contig');

    $asm_mapper = $map_adaptor->fetch_by_CoordSystems($cs1, $cs2);

    #map to contig coordinate system from chromosomal
    @ctg_coord_list = $asm_mapper->map(23142,1_000_000,2_000_000, 1, $ctg_cs);

    #map to chromosome coordinate system from contig
    @chr_coord_list = $asm_mapper->map(5431, 1000, 10_000, -1, $chr_cs);

    #list contig ids for a region of chromsome with seq_region_id 23142
    @ctg_ids = $asm_mapper->list_ids(23142, 1_000_000, 1, $chr_cs);

    #list chromosome ids for a contig
    @chr_ids = $asm_mapper->list_ids(5431, 10_000, -1, $ctg_cs);

=head1 DESCRIPTION

The AssemblyMapper is a database aware mapper which faciliates conversion
of coordinates between any two coordinate systems with an relationship
explicitly defined in the assembly table.  In the future it may be possible to
perform multiple step (implicit) mapping between coordinate systems.

It is implemented using the Bio::EnsEMBL::Mapper object, which is a generic
mapper object between disjoint coordinate systems.

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut



package Bio::EnsEMBL::AssemblyMapper;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Mapper;

@ISA = qw(Bio::EnsEMBL::Root);


my $CHUNKFACTOR = 20;
my $CHUNKSIZE   = 2**$CHUNKFACTOR;

=head2 new

  Arg [1]    : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
  Arg [2]    : Bio::EnsEMBL::CoordSystem $asm_cs
  Arg [3]    : Bio::EnsEMBL::CoordSystem $cmp_cs
  Example    : Should use AssemblyMapperAdaptor->fetch_by_CoordSystems
  Description: Creates a new AssemblyMapper
  Returntype : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
  Exceptions : thrown if multiple coord_systems are provided
  Caller     : 

=cut

sub new {
  my ($caller,$adaptor,@coord_systems) = @_;

  my $class = ref($caller) || $caller;

  my $self = {};
  bless $self, $class;

  $self->adaptor($adaptor);

  if(@coord_systems != 2) {
    throw('Can currently only map between 2 coordinate systems. ' .
          scalar(@coord_systems) . ' were provided');
  }

  $self->{'asm_cs'} = $coord_system[0];
  $self->{'cmp_cs'} = $coord_system[1];

  my $mapper = Bio::EnsEMBL::Mapper->new($coord_system[0]->dbID,
                                         $coord_system[1]->dbID);

  $self->{'mapper'} = $mapper;

  return $self;
}



sub map {
  my ($self, $frm_id, $frm_start, $frm_end, $frm_strand, $to_cs) = @_;

  #
  # Function calls are slow and this method is used in several
  # tight loops when mapping large numbers of features.
  # First just check pointer equality to see if it same coord system object
  # is used. This will work most of the time because the CoordSystemAdaptor
  # only creates a single CoordSystem object for each coord system.
  # If it doesn't work go with the slow approach of a function call
  #
  if($to_cs == $self->{'asm_cs'}) {
    if(!$self->{'cmp_register'}->{$id}) {
      $self->{'adaptor'}->register_component($self,$frm_id);
    }
    $self->{'mapper'}->map($frm_id, $frm_start,$frm_end,$frm_strand,
                           $to_cs->dbID());
  }


  } elsif($to_cs == $self->{'cmp_cs'} || $to_cs->equals($self->{'cmp_cs'})) {
    #
    # This can be probably be sped up some by only calling registered
    # assembled if needed
    #
    $self->{'adaptor'}->register_assembled($self,$frm_id,$frm_start,$frm_end);
    $self->{'mapper'}->map($frm_id,$frm_start,$frm_end,$frm_strand,
                           $to_cs->dbID());
  }

  #ptr test failed, do slow comparison
  elsif($to_cs->equals($self->{'asm_cs'})) {
    if(!$self->{'_cmp_register'}->{$id}) {
      $self->{'adaptor'}->register_component($self,$from_id);
    }
    $self->{'mapper'}->map($frm_id, $frm_start, $frm_end, $frm_strand,
                           $to_cs->dbID());

  } else {
    throw("Coordinate system " . $to_cs->name . " " . $to_cs->version .
          " is neither the assembled nor the component coordinate system " .
          " of this AssemblyMapper");
  }
}



sub list_ids {
  my($self, $frm_start, $frm_end, $frm_strand, $to_cs) = @_;

  if($to_cs->equals($self->assembled_CoordSystem())) {
    
  } elsif
}


=head2 map_coordinates_to_assembly

  Arg  1     : int $contig_id
               raw contig internal ID
  Arg  2     : int $start
               start position on raw contig
  Arg  3     : int $end
               end position on raw contig
  Arg  4     : int $strand
               raw contig orientation (+/- 1)
  Example    : none
  Description: takes RawContig coordinates and remaps it
               to Assembly coordinates. Areas that dont map produce 
               Gap objects.
  Returntype : list of Bio::EnsEMBL::Mapper::Coordinate,
               Bio::EnsEMBL::Mapper::Gap
  Exceptions : throws if args are not numeric
  Caller     : general

=cut


sub map_coordinates_to_assembly {
    my ($self, $contig_id, $start, $end, $strand) = @_;

    if( ! exists $self->{'_contig_register'}->{$contig_id} ) {
      $self->register_region_around_contig( $contig_id, 0, 0 );
    }

    return $self->{'_mapper'}->map_coordinates($contig_id, $start, 
					       $end, $strand, 'rawcontig');
}

=head2 fast_to_assembly

  Arg  1     : int $contig_id
               raw contig internal ID
  Arg  2     : int $start
               start position on raw contig
  Arg  3     : int $end
               end position on raw contig
  Arg  4     : int $strand
               raw contig orientation (+/- 1)
  Example    : none
  Description: takes RawContig coordinates and remaps it
               to Assembly coordinates. This is a fast simple version
               that only maps when there are no gaps and no splits.
               It will just return id, start, end, strand
  Returntype : list of results
  Exceptions : throws if args are not numeric
  Caller     : general

=cut


sub fast_to_assembly {
    my ($self, $contig_id, $start, $end, $strand) = @_;

    if( ! exists $self->{'_contig_register'}->{$contig_id} ) {
      $self->register_region_around_contig( $contig_id, 0, 0 );
    }
    return $self->{'_mapper'}->fastmap($contig_id, $start, $end, 
				       $strand, 'rawcontig');
}


=head2 map_coordinates_to_rawcontig

  Arg  1     : string $chr_name
  Arg  2     : int $chr_start
  Arg  3     : int $chr_end
  Arg  4     : int $chr_strand
               From p to q is 1 reverse is -1
  Example    : ( "X", 10000, 20000, -1 )
  Description: takes region in Assembly coordinates and
               remaps it to RawContig coordinates.
  Returntype : list of Bio::EnsEMBL::Mapper::Gap
               and/or   Bio::EnsEMBL::Mapper::Coordinate
	       The id method of each Coordinate object
	       is the numeric contig_id of the RawContig
	       it maps to
  Exceptions : argument type is checked where appropriat
  Caller     : general

=cut


sub map_coordinates_to_rawcontig {
  my ($self, $chr_name, $start, $end, $strand) = @_;
  
  $self->register_region($chr_name, $start, $end);
  
  return $self->_mapper->map_coordinates($chr_name, $start, $end, 
					 $strand, 'assembly');
}



=head2 list_contig_ids

  Arg  1     : string $chr_name
  Arg  2     : int $chr_start
  Arg  3     : int $chr_end
  Example    : ( "X", 1, 1000000 )
  Description: Returns a list of RawContig internal IDs
               which overlap the given chromosome region
  Returntype : list of int
  Exceptions : arguments are type checked
  Caller     : general, used for SQL query generation

=cut


sub list_contig_ids {
  my ($self, $chr_name, $start, $end) = @_;
  
  # may not have registered this region yet
  
  $self->register_region($chr_name, $start, $end);
  
  my @pairs = $self->_mapper->list_pairs($chr_name, $start, $end, 'assembly');
  
  my @ids;
  
  foreach my $pair ( @pairs ) {    
    push(@ids,$pair->from->id);
  }
  
  return @ids;
}



=head2 register_region

  Arg  1     : string $chr_name
  Arg  2     : int $chr_start
  Arg  3     : int $chr_end
  Example    : ( "X", 1, 1000000 )
  Description: Declares a chromosomal region to the AssemblyMapper.
               This extracts the relevant data from the assembly
               table and stores it in a Bio::EnsEMBL::Mapper.
               It therefore must be called before any mapping is
               attempted on that region. Otherwise only gaps will
               be returned! Is done automatically before mapping.
  Returntype : none
  Exceptions : arguments are checked
  Caller     : internal

=cut


sub register_region {
  my ($self, $chr_name, $start, $end) = @_;
  
  my $first_chunk = int( $start / $self->_chunksize() );
  my $last_chunk = int( $end / $self->_chunksize() );
  
  $self->_chunk_register_region( $chr_name, $first_chunk, $last_chunk );
}



=head2 register_region_around_contig

  Arg  1     : int $contig_id
  Arg  2     : int $upstream_bases
  Arg  3     : int $downstream_bases
  Example    : ( 2000, 0, 0 )
  Description: Declares a region around a given RawContig to this
               AssemblyMapper. Is done automatically before a mapping
               is attempted. The registering is cached, so registering 
               multiple times is cheap.
               With 0 as up and downstream tests if given contig is mappable.
               returns 1 if its is and 0 if it isnt. 
  Returntype : int 0,1
  Exceptions : all Args must be ints
  Caller     : internal

=cut


sub register_region_around_contig {
   my ($self, $contig_id, $left, $right) = @_;

   if( $self->_have_registered_contig( $contig_id ) 
       && $left == 0 && $right==0 ) {
     if( $self->_mapper->list_pairs( $contig_id, -1, -1, "rawcontig" )) {
       return 1;
     } else {
       return 0;
     }
   }
   
   my ( $chr_name, $chr_start, $chr_end ) = 
     $self->adaptor()->register_contig( $self, $self->_type, $contig_id );

   if( defined $chr_name ) {
     $self->register_region( $chr_name, $chr_start-$left, $chr_end+$right );
     return 1;
   } else {
     return 0;
   }
}


=head2 Internal functions

Internal functions

=cut



sub have_registered_component {
  my $self = shift;
  my $id = shift;

  if($self->{'cmp_register'}->{$id}) {
    return 1;
  }

  return 0;
}


sub have_registered_assembled {
  my $self = shift;
  my $id   = shift;

  if($self->{'_assembled_register'}->{$id}) {
    return 1;
  }

  return 0;
}



sub register_component {
  my $self = shift;
  my $id = shift;

  $self->{'cmp_register'}->{$id} = 1;
}


sub register_assembled {
  my $self = shift;
  my $id = shift;

  $self->{'cmp_register'}->{$id} = 1;
}



sub mapper {
  my $self = shift;
  return $self->{'mapper'};
}

sub assembled_CoordSystem {
  my $self = shift;
  return $self->{'asm_cs'};
}


sub component_CoordSystem {
  my $self = shift;
  return $self->{'cmp_cs'};
}


sub chunk_size {
  return $CHUNKSIZE;
}

sub chunk_factor {
  return $CHUNKFACTOR;
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



=head2 _chunk_register_region

  Arg  1     : string $chr_name
  Arg  2     : int $first_chunk
  Arg  3     : int $last_chunk
  Example    : none
  Description: registers the region on given chromosome from
               $first_chunk*chunksize until ($last_chunk+1)*chunksize-1
  Returntype : none
  Exceptions : none
  Caller     : register_region

=cut


sub _chunk_register_region {
  my ( $self, $chr_name, $first_chunk, $last_chunk ) = @_;

  for( my $i = $first_chunk; $i <= $last_chunk; $i++ ) {
    if( exists $self->{$chr_name}{$i} ) {
      next;
    } else {
      $self->{$chr_name}{$i} = 1;
      my $start = $i * $self->_chunksize();
      my $end = $start + $self->_chunksize() - 1;
      $self->adaptor->register_region
	( $self, $self->_type, $chr_name, $start, $end);
      
    }
  }
}


=head2 in_assembly

  Arg  1     : Bio::EnsEMBL::Clone or
               Bio::EnsEMBL::RawContig $object_in_assembly
                
  Example    : none
  Description: tests if the given Clone or RawContig object is in the 
               assembly
  Returntype : int 0,1
  Exceptions : argument type is checked
  Caller     : general

=cut


sub in_assembly {
  my ($self, $object) = @_;

  my @contigs;

  unless(ref $object) {
    $self->throw("$object is not an object reference");
  }

  if($object->isa("Bio::EnsEMBL::Clone")) {
    #get contigs from this clone
    @contigs = @{$object->get_all_Contigs()}; 
  } elsif ($object->isa("Bio::EnsEMBL::RawContig")) {
    #we already have the contig we need
    @contigs = ($object);
  } else {
    #object is not a clone or a raw contig
    $self->throw("$object is not a RawContig or Clone object");
  }

  #verify at least one of these contigs is mapped to the assembly
  foreach my $contig (@contigs) {
    if($self->register_region_around_contig( $contig->dbID(),
					     0, 0)) {
      return 1;
    }
  }

  #none of the contigs was in the assembly (golden path)
  return 0;
}


1;
