
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

Bio::EnsEMBL::AssemblyMapper - Handles mapping from raw contigs to assembly coordinates

=head1 SYNOPSIS

    $map_adaptor = $dbadaptor->get_AssmeblyMapperAdaptor();
    $mapper = $map_adaptor->fetch_by_type('UCSC');

    $mapper->register_region('chr1', 1, 100000);

    $mapper->register_region_around_contig(30001, 1000, 1000);

    my @chr_coordlist = $mapper->map_coordinates_to_assembly(627012, 2, 5, -1);

    my @raw_coordlist = $mapper->map_coordinates_to_rawcontig("1", 10002, 10020);

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

# Object preamble - inheriets from Bio::Root::RootI

use Bio::Root::RootI;
use Bio::EnsEMBL::Mapper;

@ISA = qw(Bio::Root::RootI);


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

    Arg  1      int $contig_id
                raw contig internal ID
    Arg  2      int $start
                start position on raw contig
    Arg  3      int $end
                end position on raw contig
    Arg  4      int $strand
                raw contig orientation (+/- 1)
    Function    takes RawContig coordinates and remaps it
                to Assembly coordinates
    Returntype  Bio::EnsEMBL::Mapper::Coordinate
    Exceptions  Exception thrown if arguments not of valid type
    Caller      Bio::EnsEMBL::AssemblyMapper

=cut

sub map_coordinates_to_assembly {
    my ($self, $contig_id, $start, $end, $strand) = @_;

    unless ($contig_id =~ /^\d+$/) {
        $self->throw("Expecting numeric contig_id, but got '$contig_id'");
    }

    unless ($start =~ /^\d+$/) {
        $self->throw("Expecting integer for contig start coord, but got '$start'");
    }

    unless ($end =~ /^\d+$/) {
        $self->throw("Expecting integer for contig end coord, but got '$end'");
    }

    unless ($strand =~ /^[+-]?1$/ || $strand == 0) {
        $self->throw("Expecting +/- 1 for contig strand, but got '$strand'");
    }

    return $self->_mapper->map_coordinates($contig_id, $start, $end, $strand, 'rawcontig');
}


=head2 map_coordinates_to_rawcontig

    Arg  1      char $chr_name
                chromosome name (e.g. 'X')
    Arg  2      int $start
                start position on chromosome
    Arg  3      int $end
                end position on chromosome
    Arg  4      int $strand
                orientation of seq on assembly (+/- 1)
    Function    takes region in Assembly coordinates and
                remaps it to RawContig coordinates
    Returntype  array of Bio::EnsEMBL::Mapper::Gap
                and/or   Bio::EnsEMBL::Mapper::Coordinate
		The id method of each Coordinate object
		is the numeric contig_id of the RawContig
		it maps to
    Exceptions  Exception thrown if arguments not of valid type
    Caller      Bio::EnsEMBL::AssemblyMapper

=cut

sub map_coordinates_to_rawcontig {
    my ($self, $chr_name, $start, $end, $strand) = @_;

    unless ($chr_name =~ /^\S+$/) {
        $self->throw("Expecting sensible chromosome id, but got '$chr_name'");
    }

    unless ($start =~ /^\d+$/) {
        $self->throw("Expecting integer for chromosome start, but got '$start'");
    }

    unless ($end =~ /^\d+$/) {
        $self->throw("Expecting integer for chromosome end, but got '$end'");
    }

    unless ($strand =~ /^[+-]?1$/ || $strand == 0) {
        $self->throw("Expecting +/- 1 for chromosome strand, but got '$strand'");
    }

    $self->register_region($chr_name, $start, $end);

    return $self->_mapper->map_coordinates($chr_name, $start, $end, $strand, 'assembly');
}


=head2 list_contig_ids

    Arg  1      char $chr_name
                chromosome name (e.g. 'X')
    Arg  2      int $start
                start position on chromosome
    Arg  3      int $end
                end position on chromosome
    Function    Returns a list of RawContig internal IDs
                which overlap the given chromosome region
    Returntype  @int - list of contig internal IDs
    Exceptions  none
    Caller      Bio::EnsEMBL::AssemblyMapper

=cut

sub list_contig_ids {
   my ($self, $chr_name, $start, $end) = @_;
  
   unless ($chr_name =~ /^\S+$/) {
      $self->throw("Expecting sensible chromosome id, but got '$chr_name'");
   }

   unless ($start =~ /^\d+$/) {
      $self->throw("Expecting integer for chromosome start, but got '$start'");
   }

   unless ($end =~ /^\d+$/) {
      $self->throw("Expecting integer for chromosome end, but got '$end'");
   }

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

  Arg  1      char $chr_name
              chromosome name (e.g. '2', 'X')
  Arg  2      char $start
              chromosome start coordinate
  Arg  3      char $end
              chromosome end coordinate
  Function    Declares a chromosomal region to the AssemblyMapper.
              This extracts the relevant data from the assembly
              table and stores it in a Bio::EnsEMBL::Mapper.
              It therefore must be called before any mapping is
              attempted on that region. Otherwise only gaps will
              be returned!
  Returntype  none
  Exceptions  none
  Caller      Bio::EnsEMBL::AssemblyMapper

=cut

sub register_region {
   my ($self, $chr_name, $start, $end) = @_;
   
   unless ($chr_name =~ /^\S+$/) {
      $self->throw("Expecting sensible chromosome id, but got '$chr_name'");
   }

   unless ($start =~ /^\d+$/) {
      $self->throw("Expecting integer for chromosome start, but got '$start'");
   }

   unless ($end =~ /^\d+$/) {
      $self->throw("Expecting integer for chromosome end, but got '$end'");
   }

 
   $self->adaptor->register_region($self, $self->_type, $chr_name, $start, $end);
}


=head2 register_region_around_contig

  Arg  1      int $contig_id
              contig internal ID
  Arg  2      int $left
              5 prime (chromosomal) extension
  Arg [3]     int $right
              optional 3 prime extension
	      (same as 5 prime if not defined)
  Function    Declares a chromosomal region to the AssemblyMapper
	      based around a RawContig.
              This extracts the relevant data from the assembly
              table and stores it in a Bio::EnsEMBL::Mapper.
              It therefore must be called before any mapping is
              attempted on that region. Otherwise only gaps will
              be returned!
  Returntype  none
  Exceptions  none
  Caller      Bio::EnsEMBL::AssemblyMapper

=cut

sub register_region_around_contig {
   my ($self, $contig_id, $left, $right) = @_;

   unless ($contig_id =~ /^\d+$/) {
      $self->throw("Expecting integer for RawContig id, but got '$contig_id'");
   }

   unless ($left =~ /^\d+$/) {
      $self->throw("Expecting integer for 5 prime extension, but got '$left'");
   }

   unless ($right =~ /^\d+$/) {
      $self->throw("Expecting integer for 3 prime extension, but got '$right'");
   }

   $self->adaptor->register_region_around_contig($self, $self->_type, $contig_id, $left, $right);
}


=head2 Internal functions

Internal functions

=cut


=head2 _have_registered_contig

 Title   : _have_registered_contig
 Usage   :
 Function:
 Example :
 Returns :
 Args    :


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

 Title   : _register_contig
 Usage   :
 Function:
 Example :
 Returns :
 Args    :


=cut

sub _register_contig {
   my ($self,$id) = @_;
  
   $self->{'_contig_register'}->{$id} = 1;

}


=head2 _type

 Title   : _type
 Usage   : $obj->_type($newval)
 Function:
 Example :
 Returns : value of _type
 Args    : newvalue (optional)


=cut

sub _type {
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_type'} = $value;
    }
    return $self->{'_type'};

}


=head2 _mapper

 Title   : _mapper
 Usage   : $obj->_mapper($newval)
 Function:
 Example :
 Returns : value of _mapper
 Args    : newvalue (optional)


=cut

sub _mapper {
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_mapper'} = $value;
    }
    return $self->{'_mapper'};

}



=head2 adaptor

 Title   : adaptor
 Usage   : $obj->adaptor($newval)
 Function:
 Example :
 Returns : value of adaptor
 Args    : newvalue (optional)


=cut

sub adaptor {
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'adaptor'} = $value;
    }
    return $self->{'adaptor'};

}


1;
