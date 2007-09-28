package Bio::EnsEMBL::MappedSlice;

=head1 NAME

Bio::EnsEMBL::MappedSlice

=head1 SYNOPSIS


=head1 DESCRIPTION

  Not supported Bio::EnsEMBL::Slice methods:

    All deprecated methods
    All Bio::PrimarySeqI compliance methods
    expand
    get_generic_features
    get_seq_region_id
    seq_region_Slice

  Not supported but maybe should/could:

    calculate_pi
    calculate_theta
    get_base_count
    get_by_Individual
    get_by_strain
    invert
    

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
use Bio::EnsEMBL::Mapper;
use Scalar::Util qw(weaken);

use vars qw($AUTOLOAD);


sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($adaptor, $container, $ref_slice, $name) =
    rearrange([qw(ADAPTOR CONTAINER SLICE NAME)], @_);

  # arguement check
  unless ($container and ref($container) and
          $container->isa('Bio::EnsEMBL::MappedSliceContainer')) {
    throw("Need a MappedSliceContainer.");
  }

  unless ($ref_slice and ref($ref_slice) and
          $ref_slice->isa('Bio::EnsEMBL::Slice')) {
    throw("You must provide a reference slice.");
  }
  
  my $self = {};
  bless ($self, $class);

  #
  # initialise object
  #
  
  # need to weaken reference to prevent circular reference
  weaken($self->{'container'} = $container);
  $self->{'ref_slice'} = $ref_slice;
  $self->{'adaptor'} = $adaptor if (defined($adaptor));
  $self->{'name'} = $name if (defined($name));

  $self->{'slice_mapper_pairs'} = [];

  return $self;
}


sub add_Slice_Mapper_pair {
  my $self = shift;
  my $slice = shift;
  my $mapper = shift;

  # argument check
  unless ($slice and ref($slice) and $slice->isa('Bio::EnsEMBL::Slice')) {
    throw("You must provide a slice.");
  }

  unless ($mapper and ref($mapper) and $mapper->isa('Bio::EnsEMBL::Mapper')) {
    throw("You must provide a mapper.");
  }

  push @{ $self->{'slice_mapper_pairs'} }, [ $slice, $mapper ];
  
  return $self->{'slice_mapper_pairs'};
}


sub get_all_Slice_Mapper_pairs {
  my $self = shift;
  return $self->{'slice_mapper_pairs'};
}


sub adaptor {
  my $self = shift;
  $self->{'adaptor'} = shift if (@_);
  return $self->{'adaptor'};
}


sub container {
  my $self = shift;
  weaken($self->{'container'} = shift) if (@_);
  return $self->{'container'};
}


sub name {
  my $self = shift;
  $self->{'name'} = shift if (@_);
  return $self->{'name'};
}

# return ref_slice->$_
sub seq_region_name {
  my $self = shift;
  return $self->container->ref_Slice->seq_region_name;
}

# return ref_slice->$_
sub start {
  my $self = shift;
  return $self->container->ref_Slice->start;
}

# return ref_slice->$_
sub end {
  my $self = shift;
  return $self->container->ref_Slice->end;
}

# return ref_slice->$_
sub strand {
  my $self = shift;
  return $self->container->ref_Slice->strand;
}

# container length?
sub length {
}

# return ref_slice->$_
sub seq_region_length {
  my $self = shift;
  return $self->container->ref_Slice->seq_region_length;
}

# return ref_slice->$_
sub coord_system {
  my $self = shift;
  return $self->container->ref_Slice->coord_system;
}

# return ref_slice->$_
sub coord_system_name {
  my $self = shift;
  return $self->container->ref_Slice->coord_system_name;
}

# return ref_slice->$_
sub is_toplevel {
  my $self = shift;
  return $self->container->ref_Slice->is_toplevel;
}

# container centrepoint
sub centrepoint {
}

sub seq {
}

sub subseq {
}

sub get_repeatmasked_seq {
}

sub sub_MappedSlice {
}

sub project {
}


sub get_all_DASFeatures_dsn {
  my $self = shift;

  foreach my $pair (@{ $self->get_all_Slice_Mapper_pairs }) {
    my ($slice, $mapper) = @$pair;

    # call $method on each native slice composing the MappedSlice
    my ($featref, $styleref, $segref) =
      @{ $slice->get_all_DASFeatures_dsn(@_) };
    
    use Data::Dumper;
    warn Data::Dumper::Dumper($featref);
    warn Data::Dumper::Dumper($styleref);
    warn Data::Dumper::Dumper($segref);

  }

}


sub get_all_DAS_Features {
  my $self = shift;

  foreach my $pair (@{ $self->get_all_Slice_Mapper_pairs }) {
    my ($slice, $mapper) = @$pair;
    
    print "\nFetching DAS features for ".$slice->name."\n";

    # call $method on each native slice composing the MappedSlice
    my ($featref, $styleref, $segref) = $slice->get_all_DAS_Features;
    return ($featref, $styleref, $segref);

  }

}


=head2 

  Arg[1]      : 
  Example     : 
  Description : Aggregate data gathered from composing Slices.
                This will call Slice->get_all_* and munge the results.
  Return type : 
  Exceptions  : 
  Caller      : 
  Status      : At Risk
              : under development

=cut

sub AUTOLOAD {
  my $self = shift;
  
  my $method = $AUTOLOAD;
  $method =~ s/.*:://;
  
  # AUTOLOAD should only deal with get_all_* methods
  return unless ($method =~ /^get_all_/);
  
  my @mapped_features = ();

  foreach my $pair (@{ $self->get_all_Slice_Mapper_pairs }) {
    my ($slice, $mapper) = @$pair;
    #warn $slice->name;

    # call $method on each native slice composing the MappedSlice
    my @features = @{ $slice->$method(@_) };
    
    # map features onto the MappedSlice coordinate system
    foreach my $f (@features) {
      # first map onto the ref_slice cs using the mapper from this pair
      my @ref_coords = $mapper->map_coordinates(
          $f->slice->seq_region_name,
          $f->start,
          $f->end,
          $f->strand,
          'mapped_slice'
      );

      # sanity check
      if (scalar(@ref_coords) > 1) {
        warning("Got more than one Coordinate returned, expected only one!\n");
      }

      $f->start($ref_coords[0]->start);
      $f->end($ref_coords[0]->end);
      $f->strand($f->strand*$ref_coords[0]->strand);
      $f->slice($self->container->ref_Slice);

      # then map from ref_slice to the common cs using the mapper in the 
      # MappedSliceContainer

      push @mapped_features, $f;      
    }

  }

  return \@mapped_features;
}


1;

