package Bio::EnsEMBL::DBSQL::AssemblySliceAdaptor;

=head1 NAME

Bio::EnsEMBL::DBSQL::AssemblySliceAdaptor

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS


=head1 REALTED MODULES

Bio::EnsEMBL::MappedSlice
Bio::EnsEMBL::Compara::AlignSlice
Bio::EnsEMBL::Compara::AlignSlice::Slice
Bio::EnsEMBL::AlignStrainSlice
Bio::EnsEMBL::StrainSlice

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut


use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::MappedSlice;
use Bio::EnsEMBL::Mapper;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

our @ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);
  
  return $self;
}


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

  my $slice = $container->ref_Slice;

  # project slice onto other assembly and construct MappedSlice for result
  my $mapped_slice = Bio::EnsEMBL::MappedSlice->new(
      -SLICE     => $slice,
      -CONTAINER => $container,
      -ADAPTOR   => $self,
      -NAME      => $slice->name.":mapped_$version",
  );

  my $cs_name = $slice->coord_system_name;

  foreach my $seg (@{ $slice->project($cs_name, $version) }) {
  
    my $proj_slice = $seg->to_Slice;
    
    # create a Mapper to map to/from ref_slice
    my $mapper = Bio::EnsEMBL::Mapper->new('ref_slice', 'mapped_slice');  
    
    # tell the mapper how to map this segment
    $mapper->add_map_coordinates(
        'ref_slice',
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

