
#
# Ensembl module for Bio::EnsEMBL::AssemblyMapper
#
# Cared for by Arne Stabenau <stabenau@ebi.ac.uk>
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
  
    $mapper->register_region('chr1',1,100000);

    my @chr_coordlist = $mapper->map_coordinates_to_assembly(2,5,-1,627012);
 
    my @raw_coordlist = $mapper->map_coordinates_to_rawcontig(10002,10020,1,"chr1");

    my @cid_list = $mapper->list_contig_ids(10002,10020,"chr1");

=head1 DESCRIPTION

The AssemblyMapper is a database aware mapper which handles the raw
contig to assembly mapping. It allows mappings both from raw contigs
to assembly coordinates and from assembly coordinates back to raw
contigs. 

It extracts the mapping dynamically from the database

It is implemented using the Bio::EnsEMBL::Mapper object, which is a generic
mapper object between disjoint coordinate systems.

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


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

    my @coord = $ass_mapr->map_coordinates_to_assembly(
        $raw_start, $raw_end, $contig_ori, $contig_id);
    
    # Example:
    my @coord = $ass_mapr->map_coordinates_to_assembly(
        12003, 12135, -1, 3002);

Takes a coordinate pair in RawContig space, and
remaps it to Assembly space.

Arguments are the start, end, and orientation of
the coordinate pair in RawContig coordinates, and
the contig_id (db_ID) of the RawContig.

The return values are a list of (usually one)
C<Bio::EnsEMBL::Mapper::Coordinate> objects.

=cut

sub map_coordinates_to_assembly {
    my( $self, $start, $end, $strand, $contig_id ) = @_;

    unless ($contig_id =~ /^\d+$/) {
        $self->throw("Expecting numeric contig_id, but got '$contig_id'");
    }

    return $self->_mapper->map_coordinates($start,$end,$strand,$contig_id,'rawcontig');
}

=head2 map_coordinates_to_rawcontig

    my @coord = $ass_mapr->map_coordinates_to_rawcontig(
        $ass_start, $ass_end, $ass_ori, $chr_name);
    
    # Example:
    my @coord = $ass_mapr->map_coordinates_to_rawcontig(
        13000444, 13000444, -1, 'chr2');

Takes a coordinate pair in Assembly space, and
remaps it to RawContig space.

Arguments are the start, end, and orientation of
the coordinate pair in genomic assembly
coordinates, and the name of the chromosome (or
piece of assembly).

The return values are a list of (usually one)
C<Bio::EnsEMBL::Mapper::Coordinate> objects. The
C<id> method of each Coordinate object is the
numeric contig_id (db_ID) of the RawContig it
maps to.

=cut

sub map_coordinates_to_rawcontig {
    my( $self, $start, $end, $strand, $chr_name ) = @_;

    return $self->_mapper->map_coordinates($start,$end,$strand,$chr_name,'assembly');
}



=head2 list_contig_ids

 Title   : list_contig_ids
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub list_contig_ids{
   my ($self,$start,$end,$chr) = @_;

   # may not have registered this region yet

   $self->register_region($start,$end,$chr);

   my @pairs = $self->_mapper->list_pairs($start,$end,$chr,'assembly');
   
   my @ids;

   foreach my $pair ( @pairs ) {
       push(@ids,$pair->from->id);
   }

   return @ids;
}


=head2 register_region

    $ass_mapr->register_region(12_000_000, 13_000_000, 'X');

Causes the assembly information (needed for
coordinate mapping) for a region of a chromosome
to be fetched from the database.

Arguments are chr_start, chr_end, chr_name.  The
example would load the assembly information for
the region of chromosome B<X> between 12Mbp and
13Mbp.

=cut

sub register_region{
   my ($self,$start,$end,$chr_name) = @_;

   $self->adaptor->register_region($self,$self->_type,$chr_name,$start,$end);
}





=head2 Internal functions

Internal functions

=cut


=head2 _have_loaded_contig

 Title   : _have_loaded_contig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _have_loaded_contig{
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

sub _register_contig{
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

sub _type{
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

sub _mapper{
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

sub adaptor{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'adaptor'} = $value;
    }
    return $self->{'adaptor'};

}



1;
