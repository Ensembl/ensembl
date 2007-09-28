package Bio::EnsEMBL::Utils::AssemblyProjector;

=head1 NAME

Bio::EnsEMBL::Utils::AssemblyProjector -
utility class to post-process projections from one assembly to another

=head1 SYNOPSIS

  my $assembly_projector = Bio::EnsEMBL::Utils::AssemblyProjector->new(
      -OLD_ASSEMBLY => NCBIM36,
      -NEW_ASSEMBLY => NCBIM37,
  );

  # fetch a slice on the old assembly
  my $slice = $slice_adaptor->fetch_by_region('chromosome', 1, undef, undef,
    undef, 'NCBIM36');

  my $new_slice = $assembly_projector->old_to_new($slice);

=head1 DESCRIPTION

This class implements some utility functions which apply sensible rules to the
results of projecting a feature or slice from one assembly to another.

=head1 METHODS


=head1 REALTED MODULES

The process of creating a whole genome alignment between two assemblies (which
is the basis for the use of the methods in this class) is done by a series of
scripts. Please see

  ensembl/misc-scripts/assembly/README

for a high-level description of this process, and POD in the individual scripts
for the details.

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
no warnings qw(uninitialized);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Slice;


=head2 new

  Arg [OLD_ASSEMBLY]    : name of the old assembly
  Arg [OLD_ASSEMBLY]    : name of the new assembly
  Arg [OBJECT_TYPE]     : (optional) object type ('slice' or 'feature')
  Arg [MERGE_FRAGMENTS] : (optional) Boolean - determines if segments are merged
                          to return a single object spannin all segments
                          (default: true)
  Arg [CHECK_LENGTH]    : (optional) Boolean - determines if projected objects
                          have to have same length as original (default: false)
  Example     : my $ap = Bio::EnsEMBL::Utils::AssemblyProjector->new(
                  -OLD_ASSEMBLY => NCBIM36,
                  -NEW_ASSEMBLY => NCBIM37,
                );
  Description : Constructor.
  Return type : a Bio::EnsEMBL::Utils::AssemblyProjector object
  Exceptions  : thrown on missing arguments
                thrown on invalid OBJECT_TYPE
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($old_assembly, $new_assembly, $object_type, $merge_fragments,
      $check_length) = rearrange([qw(OLD_ASSEMBLY NEW_ASSEMBLY OBJECT_TYPE
                                     MERGE_FRAGMENTS CHECK_LENGTH)], @_);

  unless ($old_assembly and $new_assembly) {
    throw("You must provide an old and new assembly name.");
  }

  if ($object_type and !
    (lc($object_type) eq 'feature' or lc($object_type) eq 'slice')) {
      throw("Type must be 'feature' or 'slice'");
  }

  my $self = {};
  bless ($self, $class);

  # initialise
  $self->{'old_assembly'} = $old_assembly;
  $self->{'new_assembly'} = $new_assembly;
  $self->{'object_type'} = $object_type;
  
  # by default, merge fragments
  $self->{'merge_fragments'} = $merge_fragments || 1;

  # by default, do not check length
  $self->{'check_length'} = $check_length || 0;

  return $self;
}


=head2 project

  Arg[1]      : Bio::EnsEMBL::Slice or Bio::EnsEMBL::Feature $object -
                the object to project
  Arg[2]      : String $to_assembly - assembly to project to
  Example     : my $new_slice = $assembly_projector->project($old_slice, 
                  'NCBIM37');
  Description : Projects a Slice or Feature to the specified assembly.

                Several tests are performed on the result to discard unwanted
                results. All projection segments have to be on the same
                seq_region and strand. If -MERGE_FRAGMENTS is set, gaps will be
                bridged by creating a single object from first_segment_start to
                last_segment_end. If -CHECK_LENGTH is set, the projected object
                will have to have the same length as the original.

                The return value of this method will always be a single object.
                
                Please see the comments in the code for details about these
                rules.
  Return type : same a Arg 1
  Exceptions  : thrown on invalid arguments
  Caller      : general, $self->old_to_new, $self->new_to_old
  Status      : At Risk
              : under development

=cut

sub project {
  my ($self, $object, $to_assembly) = @_;

  throw("Need an assembly version to project to.") unless ($to_assembly);
  throw("Need an object to project.") unless ($object and ref($object));

  my ($slice, $object_type);

  if ($object->isa('Bio::EnsEMBL::Feature')) {
    $slice = $object->feature_Slice;
    $object_type = 'feature';
  } elsif ($object->isa('Bio::EnsEMBL::Slice')) {
    $slice = $object;
    $object_type = 'slice';
  } else {
    throw("Need a Feature or Slice to project.");
  }

  # warn if trying to project to assembly version the object already is on
  if ($slice->coord_system->version eq $to_assembly) {
    warning("Assembly version to project to ($to_assembly) is the same as your object's assembly (".$slice->coord_system->version.").");
  }

  my $cs_name = $slice->coord_system_name;

  # [todo] project $slice instead? what is more efficient for Features?
  my @segments = @{ $object->project($cs_name, $to_assembly) };

  # apply rules to projection results
  # 
  # discard the projected feature/slice if
  #   1. it doesn't project at all (no segments returned)
  #   2. [unless MERGE_FRAGMENTS is set] the projection is fragmented (more
  #      than one segment)
  #   3. [if CHECK_LENGTH is set] the projection doesn't have the same length
  #      as the original feature/slice
  #   4. all segments are on same chromosome and strand

  # test (1)
  return undef unless (@segments);

  # test (2)
  return undef if (!($self->merge_fragments) and scalar(@segments) > 1);

  # test (3)
  my $first_slice = $segments[0]->to_Slice;
  my $last_slice = $segments[-1]->to_Slice;
  return undef if ($self->check_length and
    ($last_slice->end - $first_slice->start + 1) != $object->length);

  # test (4)
  my %sr_names = ();
  my %strands = ();
  foreach my $seg (@segments) {
    my $sl = $seg->to_Slice;
    $sr_names{$sl->seq_region_name}++;
    $strands{$sl->strand}++;
  }
  return undef if (scalar(keys %sr_names) > 1 or scalar(keys %strands) > 1);

  # everything looks fine, so adjust the coords of your feature/slice
  my $new_slice = $first_slice;
  $new_slice->{'end'} = $last_slice->end;

  if ($object_type eq 'slice') {
    return $new_slice;
  } else {
    $object->start($new_slice->start);
    $object->end($new_slice->end);
    $object->strand($new_slice->strand); # ????
    $object->slice($new_slice->seq_region_Slice);

    return $object;
  }
  
}

  
=head2 old_to_new

  Arg[1]      : Bio::EnsEMBL::Slice or Bio::EnsEMBL::Feature $object -
                the object to project
  Example     : my $new_slice = $assembly_projector->old_to_new($old_slice);
  Description : Projects a Slice or Feature from old to new assembly.
                This method is just a convenience wrapper for $self->project.
  Return type : same a Arg 1
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub old_to_new {
  my ($self, $object) = @_;
  return $self->project($object, $self->new_assembly);
}


=head2 new_to_old

  Arg[1]      : Bio::EnsEMBL::Slice or Bio::EnsEMBL::Feature $object -
                the object to project
  Example     : my $old_slice = $assembly_projector->new_to_old($new_slice, 1);
  Description : Projects a Slice or Feature from new to old assembly.
                This method is just a convenience wrapper for $self->project.
  Return type : same a Arg 1
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub new_to_old {
  my ($self, $object) = @_;
  return $self->project($object, $self->old_assembly);
}


#
# accessors
#
sub old_assembly {
  my $self = shift;
  $self->{'old_assembly'} = shift if (@_);
  return $self->{'old_assembly'};
}


sub new_assembly {
  my $self = shift;
  $self->{'new_assembly'} = shift if (@_);
  return $self->{'new_assembly'};
}


sub object_type {
  my $self = shift;
  $self->{'object_type'} = shift if (@_);
  return $self->{'object_type'};
}


sub merge_fragments {
  my $self = shift;
  $self->{'merge_fragments'} = shift if (@_);
  return $self->{'merge_fragments'};
}


sub check_length {
  my $self = shift;
  $self->{'check_length'} = shift if (@_);
  return $self->{'check_length'};
}


1;


