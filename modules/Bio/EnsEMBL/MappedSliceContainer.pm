package Bio::EnsEMBL::MappedSliceContainer;

=head1 NAME

Bio::EnsEMBL::MappedSliceContainer - container for mapped slices

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


# define avalable adaptormajs to use with this container
my %adaptors = map { $_ => 1 } qw(assembly align strain);


sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($ref_slice) = rearrange([qw(SLICE)], @_);

  # argument check
  unless ($ref_slice and ref($ref_slice) and
          $ref_slice->isa('Bio::EnsEMBL::Slice')) {
    throw("You must provide a reference slice.");
  }
  
  my $self = {};
  bless ($self, $class);

  # initialise object
  $self->{'ref_slice'} = $ref_slice;
  $self->{'mapped_slices'} = [];

  return $self;
}


=head2 

  Arg[1]      : 
  Example     : 
  Description : Parameterisable wrapper for all methods that set adaptors
  Return type : 
  Exceptions  : 
  Caller      : 
  Status      : At Risk
              : under development

=cut

sub set_adaptor {
  my $self = shift;
  my $type = shift;

  # argument check
  unless ($type and $adaptors{$type}) {
    throw("Missing or unknown adaptor type.");
  }

  $type = ucfirst($type);
  my $method = "set_${type}SliceAdaptor";

  return $self->$method;
}


=head2 

  Arg[1]      : 
  Example     : 
  Description : Parameterisable wrapper for all methods that get adaptors
  Return type : 
  Exceptions  : 
  Caller      : 
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


sub set_AssemblySliceAdaptor {
  my $self = shift;
  my $assembly_slice_adaptor = shift;

  unless ($assembly_slice_adaptor and ref($assembly_slice_adatpor) and
          $assembly_slice_adaptor->isa('Bio::EnsEMBL::AssemblySliceAdaptor')) {
    throw("Need a Bio::EnsEMBL::AssemblySliceAdaptor.");
  }

  $self->{'adaptors'}->{'AssemblySlice'} = $assembly_slice_adaptor;
}


sub get_AssemblySliceAdaptor {
  my $self = shift;

  unless ($self->{'adaptors'}->{'AssemblySlice'}) {
    warning("No AssemblySliceAdaptor attached to MappedSliceContainer.");
  }

  return $self->{'adaptors'}->{'AssemblySlice'};
}


sub set_AlignSliceAdaptor {
  throw("Not implemented yet!");
}


sub get_AlignSliceAdaptor {
  throw("Not implemented yet!");
}


sub set_StrainSliceAdaptor {
  throw("Not implemented yet!");
}


sub get_StrainSliceAdaptor {
  throw("Not implemented yet!");
}


sub attach_AssemblySlice {
  my $self = shift;
  my $version = shift;

  throw("Need a version.") unless ($version);

  my $asa = $self->get_AssemblySliceAdaptor;
  return unless ($asa);

  my @mapped_slices = @{ $asa->fetch_by_version($self, $version) };

  push @{ $self->{'mapped_slices'} }, @mapped_slices;
}


sub ref_Slice {
  my $self = shift;
  $self->{'ref_slice'} = shift if (@_);
  return $self->{'ref_slice'};
}


sub get_all_MappedSlices {
  my $self = shift;
  return $self->{'mapped_slices'};
}


sub sub_MappedSliceContainer {
  throw("Not implemented yet!");
}


1;

