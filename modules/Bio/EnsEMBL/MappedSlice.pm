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

Bio::EnsEMBL::MappedSlice - an object representing a mapped slice

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

This object represents a mapped slice, i.e. a slice that's attached
to a reference slice and a mapper to convert coordinates to/from the
reference. The attachment is done via a MappedSliceContainer which
has the reference slice and the "container slice" defining the common
coordinate system for all MappedSlices.

A MappedSlice is supposed to behave as close to a Bio::EnsEMBL::Slice
as possible. Most Slice methods are implemented in MappedSlice and will
return an equivalent value to what Slice does. There are some exceptions
of unimplemented methods, either because there is no useful equivalent
for a MappedSlice to do, or they are too complicated.

Not supported Bio::EnsEMBL::Slice methods:

  All deprecated methods
  All Bio::PrimarySeqI compliance methods
  expand
  get_generic_features
  get_seq_region_id
  seq_region_Slice

Not currently supported but maybe should/could:

  calculate_pi
  calculate_theta
  get_base_count
  get_by_Individual
  get_by_strain
  invert

Internally, a MappedSlice is a collection of Bio::EnsEMBL::Slices and
associated Bio::EnsEMBL::Mappers which map the slices to the common
container coordinate system.

MappedSlices are usually created and attached to a MappedSliceContainer
by an adaptor/factory.

=head1 METHODS

  new
  add_Slice_Mapper_pair
  get_all_Slice_Mapper_pairs
  adaptor
  container
  name
  seq_region_name
  start
  end
  strand
  length
  seq_region_length
  centrepoint
  coord_system
  coord_system_name
  is_toplevel
  seq (not implemented yet)
  subseq (not implemented yet)
  get_repeatmasked_seq (not implemented yet)
  sub_MappedSlice (not implemented yet)
  project (not implemented yet)

=head1 RELATED MODULES

  Bio::EnsEMBL::MappedSlice
  Bio::EnsEMBL::DBSQL::AssemblySliceAdaptor
  Bio::EnsEMBL::Compara::AlignSlice
  Bio::EnsEMBL::Compara::AlignSlice::Slice
  Bio::EnsEMBL::AlignStrainSlice
  Bio::EnsEMBL::StrainSlice

=cut

package Bio::EnsEMBL::MappedSlice;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Mapper;
use Scalar::Util qw(weaken);

use vars qw($AUTOLOAD);


=head2 new

  Arg [ADAPTOR]   : Adaptor $adaptor - an adaptor of the appropriate type
  Arg [CONTAINER] : Bio::EnsEMBL::MappedSliceContainer $container - the
                    container this MappedSlice is attached to
  Arg [NAME]      : String $name - name
  Example     : my $mapped_slice = Bio::EnsEMBL::MappedSlice->new(
                  -ADAPTOR   => $adaptor,
                  -CONTAINER => $container,
                  -NAME      => $name,
                );
  Description : Constructor. Usually you won't call this method manually, but
                the MappedSlice will be constructed by an adaptor/factory.
  Return type : Bio::EnsEMBL::MappedSlice
  Exceptions  : thrown on wrong or missing arguments
  Caller      : general, MappedSlice adaptors
  Status      : At Risk
              : under development

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($adaptor, $container, $name) =
    rearrange([qw(ADAPTOR CONTAINER NAME)], @_);

  # arguement check
  unless ($container and ref($container) and
          $container->isa('Bio::EnsEMBL::MappedSliceContainer')) {
    throw("Need a MappedSliceContainer.");
  }

  my $self = {};
  bless ($self, $class);

  #
  # initialise object
  #
  
  # need to weaken reference to prevent circular reference
  weaken($self->{'container'} = $container);

  $self->adaptor($adaptor) if (defined($adaptor));
  $self->{'name'} = $name if (defined($name));

  $self->{'slice_mapper_pairs'} = [];

  return $self;
}


=head2 add_Slice_Mapper_pair 

  Arg[1]      : Bio::EnsEMBL::Slice $slice - slice to add
  Arg[2]      : Bio::EnsEMBL::Mapper $mapper - the mapper for this slice
  Example     : $mapped_slice->add_Slice_Mapper_pair($slice, $mapper);
  Description : Adds a native slice and a corresponding mapper to map to/from
                the artificial container coord system.
  Return type : listref of Bio::EnsEMBL::MappedSlice
  Exceptions  : thrown on wrong or missing arguments
  Caller      : general, MappedSlice adaptors
  Status      : At Risk
              : under development

=cut

sub add_Slice_Mapper_pair {
  my $self = shift;
  my $slice = shift;
  my $mapper = shift;

  # argument check
  unless ($slice and ref($slice) and ($slice->isa('Bio::EnsEMBL::Slice') or $slice->isa('Bio::EnsEMBL::LRGSlice')) ) {
    throw("You must provide a slice.");
  }

  unless ($mapper and ref($mapper) and $mapper->isa('Bio::EnsEMBL::Mapper')) {
    throw("You must provide a mapper.");
  }

  push @{ $self->{'slice_mapper_pairs'} }, [ $slice, $mapper ];
  
  return $self->{'slice_mapper_pairs'};
}


=head2 get_all_Slice_Mapper_pairs 

  Example     : foreach my $pair (@{ $self->get_all_Slice_Mapper_pairs }) {
                  my ($slice, $mapper) = @$pair;

                  # get container coordinates
                  my @coords = $mapper->map_coordinates(
                    $slice->seq_region_name,
                    $slice->start,
                    $slice->end,
                    $slice->strand,
                    'mapped_slice'
                  );

                  # ....
                }
  Description : Gets all Slice/Mapper pairs this MappedSlice is composed of.
                Each slice (and features on it) can be mapped onto the
                artificial container coord system using the mapper.
  Return type : listref of listref of a Bio::EnsEMBL::Slice and
                Bio::EnsEMBL::Mapper pair
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_all_Slice_Mapper_pairs {
  my $self = shift;
  return $self->{'slice_mapper_pairs'};
}


=head2 adaptor

  Arg[1]      : (optional) Adaptor $adaptor - the adaptor/factory for this
                object
  Example     : $mapped_slice->adaptor($assembly_slice_adaptor);
  Description : Getter/setter for the adaptor/factory for this object.
  Return type : Adaptor of appropriate type
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub adaptor {
  my $self = shift;
  weaken($self->{'adaptor'} = shift) if (@_);
  return $self->{'adaptor'};
}


=head2 container

  Arg[1]      : (optional) Bio::EnsEMBL::MappedSliceContainer - the container
                this object is attached to
  Example     : my $container = $mapped_slice->container;
                print $container->ref_slice->name, "\n";
  Description : Getter/setter for the container this object is attached to. The
                container will give you access to the reference slice, a common
                artificial container slice, and a mapper to map to it from the
                container coord system.

                The implementation uses a weak reference to attach the container
                since the container holds a list of MappedSlices itself.
  Return type : Bio::EnsEMBL::MappedSliceContainer
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub container {
  my $self = shift;
  weaken($self->{'container'} = shift) if (@_);
  return $self->{'container'};
}


=head2 name

  Arg[1]      : String - the name of this object
  Example     : my $name = $mapped_slice->container->ref_slice->name .
                  ":mapped_" . $ident_string;
                $mapped_slice->name($name);
  Description : Getter/setter for this object's name
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub name {
  my $self = shift;
  $self->{'name'} = shift if (@_);
  return $self->{'name'};
}


=head2 seq_region_name

  Example     : my $sr_name = $mapped_slice->seq_region_name;
  Description : Returns the seq_region name of the reference slice.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub seq_region_name {
  my $self = shift;
  return $self->container->ref_slice->seq_region_name;
}


=head2 start

  Example     : my $start = $mapped_slice->start;
  Description : Returns the start of the container slice.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub start {
  my $self = shift;
  return $self->container->container_slice->start;
}


=head2 end

  Example     : my $end = $mapped_slice->end;
  Description : Returns the end of the container slice.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub end {
  my $self = shift;
  return $self->container->container_slice->end;
}


=head2 strand

  Example     : my $strand = $mapped_slice->strand;
  Description : Returns the strand of the container slice.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub strand {
  my $self = shift;
  return $self->container->container_slice->strand;
}


=head2 length

  Example     : my $length = $mapped_slice->length;
  Description : Returns the length of the container slice
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub length {
  my $self = shift;
  return $self->container->container_slice->length;
}


=head2 seq_region_length

  Example     : my $sr_length = $mapped_slice->seq_region_length;
  Description : Returns the seq_region length of the reference slice.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub seq_region_length {
  my $self = shift;
  return $self->container->ref_slice->seq_region_length;
}


=head2 centrepoint

  Example     : my $centrepoint = $mapped_slice->centrepoint;
  Description : Returns the centrepoint of the container slice.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub centrepoint {
  my $self = shift;
  return $self->container->container_slice->centrepoint;
}


=head2 coord_system

  Example     : my $cs = $mapped_slice->coord_system;
  Description : Returns the coord system of the container slice.
  Return type : Bio::EnsEMBL::CoordSystem
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub coord_system {
  my $self = shift;
  return $self->container->container_slice->coord_system;
}

=head2 coord_system_name

  Example     : my $cs_name = $mapped_slice->coord_system_name;
  Description : Returns the coord system name of the container slice.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub coord_system_name {
  my $self = shift;
  return $self->container->container_slice->coord_system_name;
}

=head2 is_toplevel

  Example     : my $toplevel_flag = $mapped_slice->is_toplevel;
  Description : Returns weather the container slice is toplevel.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub is_toplevel {
  my $self = shift;
  return $self->container->container_slice->is_toplevel;
}


=head2 seq

  Example     : my $seq = $mapped_slice->seq()
  Description : Retrieves the expanded sequence of this mapped slice,
                including "-" characters where there are inserts in any other
                mapped slices. This will align with the sequence returned by
                the container's seq() method.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub seq {
  my $self = shift;
  
  # create an empty string
  my $ms_seq = '';
  
  # this coord represents the current position in the MS sequence
  my $start = 0;

  # get slice/mapper pairs from mapped slice (usually only one anyway)
  foreach my $pair(@{$self->get_all_Slice_Mapper_pairs()}) {
    my ($s, $m) = @$pair;
    
    # make sure to send extra args
    # eg strain slices might need read coverage filtering
    my $seq = $s->seq(@_);
    
    # project from mapped slice to reference slice using the mapper
    foreach my $ref_coord($m->map_coordinates('mapped_slice', 1, CORE::length($seq), $s->strand, 'mapped_slice')) {
    
      # normal coord
      if(!$ref_coord->isa('Bio::EnsEMBL::Mapper::IndelCoordinate') && !$ref_coord->isa('Bio::EnsEMBL::Mapper::Gap')) {
      
        # project from reference slice to container slice using the container's mapper
        foreach my $ms_coord($self->container->mapper->map_coordinates($self->container->ref_slice->seq_region_name, $ref_coord->start, $ref_coord->end, $ref_coord->strand, 'ref_slice')) {
          
          # normal coord
          if(!$ms_coord->isa('Bio::EnsEMBL::Mapper::IndelCoordinate') && !$ms_coord->isa('Bio::EnsEMBL::Mapper::Gap')) {
            $ms_seq .= substr($seq, $start, $ms_coord->length);
            
            $start += $ms_coord->length();
          }
          
          # indel coord
          else {
            $ms_seq .= '-' x $ms_coord->length();
          }
        }
      }
      
      # indel / gap
      else {
        if($ref_coord->isa('Bio::EnsEMBL::Mapper::IndelCoordinate')) {
          if($ref_coord->gap_length > 0) {
            $ms_seq .= substr($seq, $start, $ref_coord->gap_length);
            $start += $ref_coord->gap_length;
          }
          $ms_seq .= '-' x ($ref_coord->length() - $ref_coord->gap_length());
        } else {
          $ms_seq .= '-' x $ref_coord->length();
        }
      }
    }
  }
  
  return $ms_seq;
}

sub subseq {
}

sub get_repeatmasked_seq {
}

sub sub_MappedSlice {
}

sub project {
}


=head2 AUTOLOAD

  Arg[1..N]   : Arguments passed on to the calls on the underlying slices.
  Example     : my @simple_features = @{ $mapped_slice->get_all_SimpleFeatures };
  Description : Aggregate data gathered from composing Slices.
                This will call Slice->get_all_* and combine the results.
                Coordinates will be transformed to be on the container slice
                coordinate system.

                Calls involving DAS features are skipped since the DAS adaptor
                handles coordinate conversions natively.
  Return type : listref of features (same type as corresponding Slice method)
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub AUTOLOAD {
  my $self = shift;
  
  my $method = $AUTOLOAD;
  $method =~ s/.*:://;
  
  # AUTOLOAD should only deal with get_all_* methods
  return unless ($method =~ /^get_all_/);

  # skip DAS methods
  return if ($method =~ /DAS/);
  
  my @mapped_features = ();

  foreach my $pair (@{ $self->get_all_Slice_Mapper_pairs }) {
    my ($slice, $mapper) = @$pair;
    #warn $slice->name;

    # call $method on each native slice composing the MappedSlice
    my @features = @{ $slice->$method(@_) };
    
    # map features onto the artificial container coordinate system
    foreach my $f (@features) {
      
      my @coords = $mapper->map_coordinates(
          $f->slice->seq_region_name,
          $f->start,
          $f->end,
          $f->strand,
          'mapped_slice'
      );

      # sanity check
      if (scalar(@coords) > 1) {
        warning("Got more than one Coordinate returned, expected only one!\n");
      }

      $f->start($coords[0]->start);
      $f->end($coords[0]->end);
      $f->strand($coords[0]->strand);
      $f->slice($self->container->container_slice);

      push @mapped_features, $f;
    }

  }

  return \@mapped_features;
}


1;

