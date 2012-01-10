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

Bio::EnsEMBL::MappedSliceContainer - container for mapped slices

=head1 SYNOPSIS

  # get a reference slice
  my $slice =
    $slice_adaptor->fetch_by_region( 'chromosome', 14, 900000, 950000 );

  # create MappedSliceContainer based on the reference slice
  my $msc = Bio::EnsEMBL::MappedSliceContainer->new( -SLICE => $slice );

  # set the adaptor for fetching AssemblySlices
  my $asa = $slice->adaptor->db->get_AssemblySliceAdaptor;
  $msc->set_AssemblySliceAdaptor($asa);

  # add an AssemblySlice to your MappedSliceContainer
  $msc->attach_AssemblySlice('NCBIM36');

  foreach my $mapped_slice ( @{ $msc->get_all_MappedSlices } ) {
    print $mapped_slice->name, "\n";

    foreach my $sf ( @{ $mapped_slice->get_all_SimpleFeatures } ) {
      print "  ", &to_string($sf), "\n";
    }
  }

=head1 DESCRIPTION

NOTE: this code is under development and not fully functional nor tested
yet.  Use only for development.

A MappedSliceContainer holds a collection of one or more
Bio::EnsEMBL::MappedSlices. It is based on a real reference slice and
contains an artificial "container slice" which defines the common
coordinate system used by all attached MappedSlices. There is also a
mapper to convert coordinates between the reference and the container
slice.

Attaching MappedSlices to the container is delegated to adaptors
(which act more as object factories than as traditional Ensembl db
adaptors). The adaptors will also modify the container slice and
associated mapper if required. This design allows us to keep the
MappedSliceContainer generic and encapsulate the data source specific
code in the adaptor/factory module.

In the simplest use case, all required MappedSlices are attached to the
MappedSliceContainer at once (by a single call to the adaptor). This
object should also allow "hot-plugging" of MappedSlices (e.g. attach a
MappedSlice representing a strain to a container that already contains a
multi-species alignment). The methods for attaching new MappedSlice will
be responsable to perform the necessary adjustments to coordinates and
mapper on the existing MappedSlices.

=head1 METHODS

  new
  set_adaptor
  get_adaptor
  set_AssemblySliceAdaptor
  get_AssemblySliceAdaptor
  set_AlignSliceAdaptor (not implemented yet)
  get_AlignSliceAdaptor (not implemented yet)
  set_StrainSliceAdaptor (not implemented yet)
  get_StrainSliceAdaptor (not implemented yet)
  attach_AssemblySlice
  attach_AlignSlice (not implemented yet)
  attach_StrainSlice (not implemented yet)
  get_all_MappedSlices
  sub_MappedSliceContainer (not implemented yet)
  ref_slice
  container_slice
  mapper
  expanded

=head1 RELATED MODULES

  Bio::EnsEMBL::MappedSlice
  Bio::EnsEMBL::DBSQL::AssemblySliceAdaptor
  Bio::EnsEMBL::Compara::AlignSlice
  Bio::EnsEMBL::Compara::AlignSlice::Slice
  Bio::EnsEMBL::AlignStrainSlice
  Bio::EnsEMBL::StrainSlice

=cut

package Bio::EnsEMBL::MappedSliceContainer;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Mapper;


# define avalable adaptormajs to use with this container
my %adaptors = map { $_ => 1 } qw(assembly align strain);


=head2 new

  Arg [SLICE]     : Bio::EnsEMBL::Slice $slice - the reference slice for this
                    container
  Arg [EXPANDED]  : (optional) Boolean $expanded - set expanded mode (default:
                    collapsed)
  Example     : my $slice = $slice_adaptor->fetch_by_region('chromosome', 1,
                  9000000, 9500000);
                my $msc = Bio::EnsEMBL::MappedSliceContainer->new(
                    -SLICE    => $slice,
                    -EXPANDED => 1,
                );
  Description : Constructor. See the general documentation of this module for 
                details about this object. Note that the constructor creates an
                empty container, so you'll have to attach MappedSlices to it to
                be useful (this is usually done by an adaptor/factory).
  Return type : Bio::EnsEMBL::MappedSliceContainer
  Exceptions  : thrown on wrong or missing argument
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($ref_slice, $expanded) = rearrange([qw(SLICE EXPANDED)], @_);

  # argument check
  unless ($ref_slice and ref($ref_slice) and
          ($ref_slice->isa('Bio::EnsEMBL::Slice') or $ref_slice->isa('Bio::EnsEMBL::LRGSlice')) ) {
    throw("You must provide a reference slice.");
  }

  my $self = {};
  bless ($self, $class);

  # initialise object
  $self->{'ref_slice'} = $ref_slice;
  $self->{'expanded'} = $expanded || 0;

  $self->{'mapped_slices'} = [];

  # create the container slice
  $self->_create_container_slice($ref_slice);

  return $self;
}


#
# Create an artificial slice which represents the common coordinate system used
# for this MappedSliceContainer
#
sub _create_container_slice {
  my $self = shift;
  my $ref_slice = shift;

  # argument check
  unless ($ref_slice and ref($ref_slice) and
          ($ref_slice->isa('Bio::EnsEMBL::Slice') or $ref_slice->isa('Bio::EnsEMBL::LRGSlice')) ) {
    throw("You must provide a reference slice.");
  }

  # create an artificial coordinate system for the container slice
  my $cs = Bio::EnsEMBL::CoordSystem->new(
      -NAME => 'container',
      -RANK => 1,
  );

  # Create a new artificial slice spanning your container. Initially this will
  # simply span your reference slice
  my $container_slice = Bio::EnsEMBL::Slice->new(
      -COORD_SYSTEM     => $cs,
      -START            => 1,
      -END              => $ref_slice->length,
      -STRAND           => 1,
      -SEQ_REGION_NAME  => 'container',
  );

  $self->{'container_slice'} = $container_slice;

  # Create an Mapper to map to/from the reference slice to the container coord
  # system.
  my $mapper = Bio::EnsEMBL::Mapper->new('ref_slice', 'container');
  
  $mapper->add_map_coordinates(
      $ref_slice->seq_region_name,
      $ref_slice->start,
      $ref_slice->end,
      1,
      $container_slice->seq_region_name,
      $container_slice->start,
      $container_slice->end,
  );

  $self->{'mapper'} = $mapper;
}


=head2 set_adaptor

  Arg[1]      : String $type - the type of adaptor to set
  Arg[2]      : Adaptor $adaptor - the adaptor to set
  Example     : my $adaptor = Bio::EnsEMBL::DBSQL::AssemblySliceAdaptor->new;
                $msc->set_adaptor('assembly', $adaptor);
  Description : Parameterisable wrapper for all methods that set adaptors (see
                below).
  Return type : same as Arg 2
  Exceptions  : thrown on missing type
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub set_adaptor {
  my $self = shift;
  my $type = shift;
  my $adaptor = shift;

  # argument check
  unless ($type and $adaptors{$type}) {
    throw("Missing or unknown adaptor type.");
  }

  $type = ucfirst($type);
  my $method = "set_${type}SliceAdaptor";

  return $self->$method($adaptor);
}


=head2 get_adaptor

  Arg[1]      : String $type - the type of adaptor to get
  Example     : my $assembly_slice_adaptor = $msc->get_adaptor('assembly');
  Description : Parameterisable wrapper for all methods that get adaptors (see
                below).
  Return type : An adaptor for the requested type of MappedSlice.
  Exceptions  : thrown on missing type
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_adaptor {
  my $self = shift;
  my $type = shift;

  # argument check
  unless ($type and $adaptors{$type}) {
    throw("Missing or unknown adaptor type.");
  }

  $type = ucfirst($type);
  my $method = "get_${type}SliceAdaptor";

  return $self->$method;
}


=head2 set_AssemblySliceAdaptor

  Arg[1]      : Bio::EnsEMBL::DBSQL::AssemblySliceAdaptor - the adaptor to set
  Example     : my $adaptor = Bio::EnsEMBL::DBSQL::AssemblySliceAdaptor->new;
                $msc->set_AssemblySliceAdaptor($adaptor);
  Description : Sets an AssemblySliceAdaptor for this container. The adaptor can
                be used to attach MappedSlice for alternative assemblies.
  Return type : Bio::EnsEMBL::DBSQL::AssemblySliceAdaptor
  Exceptions  : thrown on wrong or missing argument
  Caller      : general, $self->get_adaptor
  Status      : At Risk
              : under development

=cut

sub set_AssemblySliceAdaptor {
  my $self = shift;
  my $assembly_slice_adaptor = shift;

  unless ($assembly_slice_adaptor and ref($assembly_slice_adaptor) and
    $assembly_slice_adaptor->isa('Bio::EnsEMBL::DBSQL::AssemblySliceAdaptor')) {
      throw("Need a Bio::EnsEMBL::AssemblySliceAdaptor.");
  }

  $self->{'adaptors'}->{'AssemblySlice'} = $assembly_slice_adaptor;
}


=head2 get_AssemblySliceAdaptor

  Example     : my $assembly_slice_adaptor = $msc->get_AssemblySliceAdaptor;
  Description : Gets a AssemblySliceAdaptor from this container. The adaptor can
                be used to attach MappedSlice for alternative assemblies.
  Return type : Bio::EnsEMBL::DBSQL::AssemblySliceAdaptor
  Exceptions  : thrown on wrong or missing argument
  Caller      : general, $self->get_adaptor
  Status      : At Risk
              : under development

=cut

sub get_AssemblySliceAdaptor {
  my $self = shift;

  unless ($self->{'adaptors'}->{'AssemblySlice'}) {
    warning("No AssemblySliceAdaptor attached to MappedSliceContainer.");
  }

  return $self->{'adaptors'}->{'AssemblySlice'};
}


# [todo]
sub set_AlignSliceAdaptor {
  throw("Not implemented yet!");
}


# [todo]
sub get_AlignSliceAdaptor {
  throw("Not implemented yet!");
}


# [todo]
sub set_StrainSliceAdaptor {
  my $self = shift;
  my $strain_slice_adaptor = shift;

  unless ($strain_slice_adaptor and ref($strain_slice_adaptor) and
    $strain_slice_adaptor->isa('Bio::EnsEMBL::DBSQL::StrainSliceAdaptor')) {
      throw("Need a Bio::EnsEMBL::StrainSliceAdaptor.");
  }

  $self->{'adaptors'}->{'StrainSlice'} = $strain_slice_adaptor;
}


# [todo]
sub get_StrainSliceAdaptor {
  my $self = shift;

  unless ($self->{'adaptors'}->{'StrainSlice'}) {
    warning("No StrainSliceAdaptor attached to MappedSliceContainer.");
  }

  return $self->{'adaptors'}->{'StrainSlice'};
}


=head2 attach_AssemblySlice

  Arg[1]      : String $version - assembly version to attach
  Example     : $msc->attach_AssemblySlice('NCBIM36');
  Description : Attaches a MappedSlice for an alternative assembly to this
                container.
  Return type : none
  Exceptions  : thrown on missing argument
  Caller      : general, Bio::EnsEMBL::DBSQL::AssemblySliceAdaptor
  Status      : At Risk
              : under development

=cut

sub attach_AssemblySlice {
  my $self = shift;
  my $version = shift;

  throw("Need a version.") unless ($version);

  my $asa = $self->get_AssemblySliceAdaptor;
  return unless ($asa);

  my @mapped_slices = @{ $asa->fetch_by_version($self, $version) };

  push @{ $self->{'mapped_slices'} }, @mapped_slices;
}


=head2 attach_StrainSlice

  Arg[1]      : String $strain - name of strain to attach
  Example     : $msc->attach_StrainSlice('Watson');
  Description : Attaches a MappedSlice for an alternative strain to this
                container.
  Return type : none
  Exceptions  : thrown on missing argument
  Caller      : general, Bio::EnsEMBL::DBSQL::StrainSliceAdaptor
  Status      : At Risk
              : under development

=cut

sub attach_StrainSlice {
  my $self = shift;
  my $strain = shift;

  throw("Need a strain.") unless ($strain);

  my $ssa = $self->get_StrainSliceAdaptor;
  return unless ($ssa);

  my @mapped_slices = @{ $ssa->fetch_by_name($self, $strain) };

  push @{ $self->{'mapped_slices'} }, @mapped_slices;
}



=head2 get_all_MappedSlices

  Example     : foreach my $mapped_slice (@{ $msc->get_all_MappedSlices }) {
                  print $mapped_slice->name, "\n";
                }
  Description : Returns all MappedSlices attached to this container.
  Return type : listref of Bio::EnsEMBL::MappedSlice
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_all_MappedSlices {
  my $self = shift;
  return $self->{'mapped_slices'};
}


# [todo]
sub sub_MappedSliceContainer {
  throw("Not implemented yet!");
}


=head2 ref_slice

  Arg[1]      : (optional) Bio::EnsEMBL::Slice - the reference slice to set
  Example     : my $ref_slice = $mapped_slice_container->ref_slice;
                print "This MappedSliceContainer is based on the reference
                  slice ", $ref_slice->name, "\n";
  Description : Getter/setter for the reference slice.
  Return type : Bio::EnsEMBL::Slice
  Exceptions  : thrown on wrong argument type
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub ref_slice {
  my $self = shift;
  
  if (@_) {
    my $slice = shift;
    
    unless (ref($slice) and ($slice->isa('Bio::EnsEMBL::Slice') or $slice->isa('Bio::EnsEMBL::LRGSlice'))) {
      throw("Need a Bio::EnsEMBL::Slice.");
    }
    
    $self->{'ref_slice'} = $slice;
  }

  return $self->{'ref_slice'};
}


=head2 container_slice

  Arg[1]      : (optional) Bio::EnsEMBL::Slice - the container slice to set
  Example     : my $container_slice = $mapped_slice_container->container_slice;
                print "The common slice used by this MappedSliceContainer is ",
                  $container_slice->name, "\n";
  Description : Getter/setter for the container slice. This is an artificial
                slice which defines the common coordinate system used by the
                MappedSlices attached to this container.
  Return type : Bio::EnsEMBL::Slice
  Exceptions  : thrown on wrong argument type
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub container_slice {
  my $self = shift;
  
  if (@_) {
    my $slice = shift;
    
    unless (ref($slice) and ($slice->isa('Bio::EnsEMBL::Slice') or $slice->isa('Bio::EnsEMBL::LRGSlice')) ) {
      throw("Need a Bio::EnsEMBL::Slice.");
    }
    
    $self->{'container_slice'} = $slice;
  }

  return $self->{'container_slice'};
}


=head2 mapper

  Arg[1]      : (optional) Bio::EnsEMBL::Mapper - the mapper to set
  Example     : my $mapper = Bio::EnsEMBL::Mapper->new('ref', 'mapped');
                $mapped_slice_container->mapper($mapper);
  Description : Getter/setter for the mapper to map between reference slice and
                the artificial container coord system.
  Return type : Bio::EnsEMBL::Mapper
  Exceptions  : thrown on wrong argument type
  Caller      : internal, Bio::EnsEMBL::MappedSlice->AUTOLOAD
  Status      : At Risk
              : under development

=cut

sub mapper {
  my $self = shift;
  
  if (@_) {
    my $mapper = shift;
    
    unless (ref($mapper) and $mapper->isa('Bio::EnsEMBL::Mapper')) {
      throw("Need a Bio::EnsEMBL::Mapper.");
    }
    
    $self->{'mapper'} = $mapper;
  }

  return $self->{'mapper'};
}


=head2 expanded

  Arg[1]      : (optional) Boolean - expanded mode to set
  Example     : if ($mapped_slice_container->expanded) {
                  # do more elaborate mapping than in collapsed mode
                  [...]
                }
  Description : Getter/setter for expanded mode.
                
                By default, MappedSliceContainer use collapsed mode, which
                means that no inserts in the reference sequence are allowed
                when constructing the MappedSlices. in this mode, the
                mapped_slice artificial coord system will be identical with the
                ref_slice coord system.
                
                By setting expanded mode, you allow inserts in the reference
                sequence.
  Return type : Boolean
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub expanded {
  my $self = shift;
  $self->{'expanded'} = shift if (@_);
  return $self->{'expanded'};
}

=head2 seq

  Example     : my $seq = $container->seq()
  Description : Retrieves the expanded sequence of the artificial container
                slice, including "-" characters where there are inserts in any
                of the attached mapped slices.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub seq {
  my $self = shift;
  
  my $container_seq = '';
  
  # check there's a mapper
  if(defined($self->mapper)) {
    my $start = 0;
    my $slice = $self->ref_slice();
    my $seq = $slice->seq();
    
    foreach my $coord($self->mapper->map_coordinates($slice->seq_region_name, $slice->start, $slice->end, $slice->strand, 'ref_slice')) {
      # if it is a normal coordinate insert sequence
      if(!$coord->isa('Bio::EnsEMBL::Mapper::IndelCoordinate')) {
        $container_seq .= substr($seq, $start, $coord->length());
        $start += $coord->length;
      }
      
      # if it is a gap or indel insert "-"
      else {
        $container_seq .= '-' x $coord->length();
      }
    }
  }
  
  return $container_seq;
}


1;

