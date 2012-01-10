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

Bio::EnsEMBL::DBSQL::AssemblySliceAdaptor - adaptor/factory for MappedSlices
representing alternative assemblies

=head1 SYNOPSIS

  my $slice =
    $slice_adaptor->fetch_by_region( 'chromosome', 14, 900000, 950000 );

  my $msc = Bio::EnsEMBL::MappedSliceContainer->new( -SLICE => $slice );

  my $asa = Bio::EnsEMBL::DBSQL::AssemblySliceAdaptor->new;

  my ($mapped_slice) = @{ $asa->fetch_by_version( $msc, 'NCBIM36' ) };

=head1 DESCRIPTION

NOTE: this code is under development and not fully functional nor tested
yet.  Use only for development.

This adaptor is a factory for creating MappedSlices representing
alternative assemblies and attaching them to a MappedSliceContainer. A
mapper will be created to map between the reference slice and the common
container slice coordinate system.

=head1 METHODS

  new
  fetch_by_version

=head1 REALTED MODULES

  Bio::EnsEMBL::MappedSlice
  Bio::EnsEMBL::MappedSliceContainer
  Bio::EnsEMBL::Compara::AlignSlice
  Bio::EnsEMBL::Compara::AlignSlice::Slice
  Bio::EnsEMBL::AlignStrainSlice
  Bio::EnsEMBL::StrainSlice

=cut

package Bio::EnsEMBL::DBSQL::AssemblySliceAdaptor;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::MappedSlice;
use Bio::EnsEMBL::Mapper;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

our @ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 new

  Example     : my $assembly_slice_adaptor =
                  Bio::EnsEMBL::DBSQL::AssemblySliceAdaptor->new;
  Description : Constructor.
  Return type : Bio::EnsEMBL::DBSQL::AssemblySliceAdaptor
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);
  
  return $self;
}


=head2 fetch_by_version

  Arg[1]      : Bio::EnsEMBL::MappedSliceContainer $container - the container
                  to attach MappedSlices to
  Arg[2]      : String $version - the assembly version to fetch
  Example     : my ($mapped_slice) = @{ $msc->fetch_by_version('NCBIM36') };
  Description : Creates a MappedSlice representing an alternative assembly
                version of the container's reference slice.
  Return type : listref of Bio::EnsEMBL::MappedSlice
  Exceptions  : thrown on wrong or missing arguments
  Caller      : general, Bio::EnsEMBL::MappedSliceContainer
  Status      : At Risk
              : under development

=cut

sub fetch_by_version {
  my $self = shift;
  my $container = shift;
  my $version = shift;

  # arguement check
  unless ($container and ref($container) and
          $container->isa('Bio::EnsEMBL::MappedSliceContainer')) {
    throw("Need a MappedSliceContainer.");
  }

  unless ($version) {
    throw("Need an assembly version.");
  }

  my $slice = $container->ref_slice;

  # project slice onto other assembly and construct MappedSlice for result
  my $mapped_slice = Bio::EnsEMBL::MappedSlice->new(
      -ADAPTOR   => $self,
      -CONTAINER => $container,
      -NAME      => $slice->name."\#mapped_$version",
  );

  my $cs_name = $slice->coord_system_name;

  foreach my $seg (@{ $slice->project($cs_name, $version) }) {
  
    my $proj_slice = $seg->to_Slice;
    
    # create a Mapper to map to/from the mapped_slice artificial coord system
    my $mapper = Bio::EnsEMBL::Mapper->new('mapped_slice', 'native_slice');  
    
    # tell the mapper how to map this segment
    $mapper->add_map_coordinates(
        'mapped_slice',
        $seg->from_start,
        $seg->from_end,
        ($slice->strand * $proj_slice->strand),
        $proj_slice->seq_region_name,
        $proj_slice->start,
        $proj_slice->end
    );
    
    # add the Slice/Mapper pair to the MappedSlice
    $mapped_slice->add_Slice_Mapper_pair($proj_slice, $mapper);
  }

  return [$mapped_slice];
}


=head2 fetch_by_name

  Arg[1]      : Bio::EnsEMBL::MappedSliceContainer $container - the container
                  to attach MappedSlices to
  Arg[2]      : String $name - the assembly name to fetch
  Arg[3]      : (optional) String $version -- the version for the new assembly
  Example     : my ($mapped_slice) = @{ $msc->fetch_by_name('LRG1','1') };
  Description : Creates a MappedSlice representing an alternative assembly
                version of the container's reference slice.
  Return type : listref of Bio::EnsEMBL::MappedSlice
  Exceptions  : thrown on wrong or missing arguments
  Caller      : general, Bio::EnsEMBL::MappedSliceContainer
  Status      : At Risk
              : under development

=cut

sub fetch_by_name {
  my $self = shift;
  my $container = shift;
  my $name = shift;
  my $version = shift; 

  # arguement check
  unless ($container and ref($container) and
          $container->isa('Bio::EnsEMBL::MappedSliceContainer')) {
    throw("Need a MappedSliceContainer.");
  }

  unless ($name) {
    throw("Need an assembly name.");
  }

  $version ||= '';
  my $slice = $container->ref_slice;

  # project slice onto other assembly and construct MappedSlice for result
  my $mapped_slice = Bio::EnsEMBL::MappedSlice->new(
      -ADAPTOR   => $self,
      -CONTAINER => $container,
      -NAME      => $slice->name."\#mapped_$name:$version",
  );


  foreach my $seg (@{ $slice->project($name, $version) }) {
  
    my $proj_slice = $seg->to_Slice;
    
    # create a Mapper to map to/from the mapped_slice artificial coord system
    my $mapper = Bio::EnsEMBL::Mapper->new('mapped_slice', 'native_slice');  
    
    # tell the mapper how to map this segment
    $mapper->add_map_coordinates(
        'mapped_slice',
        $seg->from_start,
        $seg->from_end,
        ($slice->strand * $proj_slice->strand),
        $proj_slice->seq_region_name,
        $proj_slice->start,
        $proj_slice->end
    );
    
    # add the Slice/Mapper pair to the MappedSlice
    $mapped_slice->add_Slice_Mapper_pair($proj_slice, $mapper);
  }

  return [$mapped_slice];
}


1;

