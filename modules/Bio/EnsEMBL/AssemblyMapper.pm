
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

    my @chr_coordlist = $mapper->map_coordinates(2,5,-1,627012,"rawcontig");
 
    my @raw_coordlist = $mapper->map_coordinates(10002,10020,1,"chr1","assembly");
    

=head1 DESCRIPTION

The AssemblyMapper is a database aware mapper which handles the raw
contig to assembly mapping. It allows mappings both from raw contigs
to assembly coordinates and from assembly coordinates back to raw
contigs. 

It extracts the mapping dynamically from the database

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



=head2 map_coordinates

 Title   : map_coordinates
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub map_coordinates{
   my ($self,$start,$end,$strand,$id,$maptype) = @_;

   if( !defined $maptype ) {
       $self->throw("map_coordinates start,end,strand,id,maptype");
   }

   if( $maptype ne "rawcontig" && $maptype ne "assembly" ) {
       $self->throw("maptype must be either rawcontig or assembly, not $maptype");
   }

   if( $maptype eq "rawcontig" && $id =~ /^\w/ ) {
       $self->throw("You have a rawcontig type, but the id is $id. AssemblyMappers expect internal_ids of the contigs in their mappers - probably you have given us a raw contig text id");
   }

   # ok map it 

   return $self->_mapper->map_coordinates($start,$end,$strand,$id,$maptype);
}


=head2 register_region

 Title   : register_region
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub register_region{
   my ($self,$id,$start,$end) = @_;

   $self->adaptor->register_region($self,$self->_type,$id,$start,$end);
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
