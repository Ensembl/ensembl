
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
Handles mapping from raw contigs to assembly coordinates

=head1 SYNOPSIS

    $map_adaptor = $dbadaptor->get_AssemblyMapperAdaptor();
    $mapper = $map_adaptor->fetch_by_type('UCSC');

    $mapper->register_region('chr1', 1, 100000);

    $mapper->register_region_around_contig(30001, 1000, 1000);

    my @chr_coordlist = $mapper->map_coordinates_to_assembly(627012, 2, 5, -1);

    my @raw_coordlist = $mapper->map_coordinates_to_rawcontig("1",10002,10020);

    my @cid_list = $mapper->list_contig_ids("1", 10002, 10020);

=head1 DESCRIPTION

The AssemblyMapper is a database aware mapper which handles the raw
contig to assembly mapping. It allows mappings both from raw contigs
to assembly coordinates and from assembly coordinates back to raw
contigs.

It extracts the mapping dynamically from the database

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

# Object preamble - inheriets from Bio::EnsEMBL::Root

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Mapper;

@ISA = qw(Bio::EnsEMBL::Root);


=head2 new

  Arg [1]    : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
  Arg [2]    : string $type
               The type of assembly mapper
  Example    : $mapper = new Bio::EnsEMBL::AssemblyMapper($adaptor,'assembly');
  Description: Creates a new AssemblyMapper object.
  Returntype : Bio::EnsEMBL::AssemblyMapper;
  Exceptions : thrown if arg count is wrong
  Caller     : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor

=cut

sub new {
  my($class,$adaptor,$type) = @_;

  my $self = {};
  bless $self,$class;

  if( !defined $type ) {
      $self->throw("Must have adaptor, type for the assembly mapper");
  }
  $self->{'_contig_register'} = {};
  $self->adaptor($adaptor);
  $self->_type($type);

  my $mapper = Bio::EnsEMBL::Mapper->new('rawcontig','assembly');
  $self->_mapper($mapper);
  return $self;
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
  
  return map
    { $_->from->id }
    $self->_mapper->list_pairs($chr_name, $start, $end, 'assembly');
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



=head2 _have_registered_contig

  Arg  1     : int $rawContig_dbID
  Example    : none
  Description: checks if given raw contig was registered before
  Returntype : int 0,1
  Exceptions : none
  Caller     : internal

=cut



sub _have_registered_contig {
   my ($self,$id) = @_;

   if( $self->{'_contig_register'}->{$id} ) {
       return 1;
   } else {
       return 0;
   }

}


=head2 _register_contig

  Arg  1     : int $rawContig_dbID
  Example    : none
  Description: marks given raw contig as registered
  Returntype : none
  Exceptions : none
  Caller     : internal

=cut



sub _register_contig {
   my ($self,$id) = @_;
  
   $self->{'_contig_register'}->{$id} = 1;

}



=head2 _type

  Arg [1]    : string $type
  Example    : none
  Description: get/set of attribute _type
  Returntype : string
  Exceptions : none
  Caller     : ?

=cut


sub _type {
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_type'} = $value;
    }
    return $self->{'_type'};

}


=head2 _mapper

  Arg [1]    : Bio::EnsEMBL::Mapper $mapper
               The mapper object which will be used to map between
               assembly and raw contigs
  Example    : none
  Description: get/set of the attribute _mapper 
  Returntype : Bio::EnsEMBL::Mapper
  Exceptions : none
  Caller     : ?

=cut


sub _mapper {
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_mapper'} = $value;
    }
    return $self->{'_mapper'};

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
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'adaptor'} = $value;
    }
    return $self->{'adaptor'};

}


# this function should be customized for the assemblies that are in 
# EnsEMBL. (Could be a hash of assembly_type ->size )


=head2 _chunksize

  Args       : none
  Example    : none
  Description: returns the size in which this mapper registers chromosomal
               regions. It should be optimized on the size of pieces 
               in assembly table, so that not too many and not too little
               rows are retrieved. 
  Returntype : int
  Exceptions : none
  Caller     : internal

=cut


sub _chunksize {
  return 1000000;
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

  foreach my $i ( $first_chunk..$last_chunk ) {
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
