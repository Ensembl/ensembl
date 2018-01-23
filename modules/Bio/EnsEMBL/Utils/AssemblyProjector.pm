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

Bio::EnsEMBL::Utils::AssemblyProjector -
utility class to post-process projections from one assembly to another

=head1 SYNOPSIS

  # connect to an old database
  my $dba_old = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => 'ensembldb.ensembl.org',
    -port   => 3306,
    -user   => ensro,
    -dbname => 'mus_musculus_core_46_36g',
    -group  => 'core_old',
  );

  # connect to the new database containing the mapping between old and
  # new assembly
  my $dba_new = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => 'ensembldb.ensembl.org',
    -port   => 3306,
    -user   => ensro,
    -dbname => 'mus_musculus_core_47_37',
    -group  => 'core_new',
  );

  my $assembly_projector = Bio::EnsEMBL::Utils::AssemblyProjector->new(
    -OLD_ASSEMBLY    => 'NCBIM36',
    -NEW_ASSEMBLY    => 'NCBIM37',
    -ADAPTOR         => $dba_new,
    -EXTERNAL_SOURCE => 1,
    -MERGE_FRAGMENTS => 1,
    -CHECK_LENGTH    => 0,
  );

  # fetch a slice on the old assembly
  my $slice_adaptor = $dba_old->get_SliceAdaptor;
  my $slice =
    $slice_adaptor->fetch_by_region( 'chromosome', 1, undef, undef,
    undef, 'NCBIM36' );

  my $new_slice = $assembly_projector->old_to_new($slice);

  print $new_slice->name, " (", $assembly_projector->last_status, ")\n";

=head1 DESCRIPTION

This class implements some utility functions for converting coordinates
between assemblies. A mapping between the two assemblies has to present
the database for this to work, see the 'Related Modules' section below
on how to generate the mapping.

In addition to the "raw" projecting of features and slices, the methods
in this module also apply some sensible rules to the results of the
projection (like discarding unwanted results or merging fragmented
projections). These are the rules (depending on configuration):

Discard the projected feature/slice if:

  1. it doesn't project at all (no segments returned)
  2. [unless MERGE_FRAGMENTS is set] the projection is fragmented (more
     than one segment)
  3. [if CHECK_LENGTH is set] the projection doesn't have the same
     length as the original feature/slice
  4. all segments are on same chromosome and strand

If a projection fails any of these rules, undef is returned instead of
a projected feature/slice. You can use the last_status() method to find
out about the results of the rules tests.

Also note that when projecting features, only a shallow projection is
performed, i.e. other features attached to your features (e.g. the
transcripts of a gene) are not projected automatically, so it will be
the responsability of the user code project all levels of features
involved.

=head1 METHODS

  new
  project
  old_to_new
  new_to_old
  adaptor
  external_source
  old_assembly
  new_assembly
  merge_fragments
  check_length

=head1 RELATED MODULES

The process of creating a whole genome alignment between two assemblies
(which is the basis for the use of the methods in this class) is done by
a series of scripts. Please see

  ensembl/misc-scripts/assembly/README

for a high-level description of this process, and POD in the individual
scripts for the details.

=cut

package Bio::EnsEMBL::Utils::AssemblyProjector;

use strict;
use warnings;
no warnings qw(uninitialized);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Slice;
use Scalar::Util qw(weaken);

=head2 new

  Arg [ADAPTOR]         : Bio::EnsEMBL::DBSQL::DBAdaptor $adaptor - a db adaptor
                          for a database containing the assembly mapping
  Arg [EXTERNAL_SOURCE] : (optional) Boolean $external_source - indicates if
                          source is from a different database
  Arg [OLD_ASSEMBLY]    : name of the old assembly
  Arg [OLD_ASSEMBLY]    : name of the new assembly
  Arg [OBJECT_TYPE]     : (optional) object type ('slice' or 'feature')
  Arg [MERGE_FRAGMENTS] : (optional) Boolean - determines if segments are merged
                          to return a single object spanning all segments
                          (default: true)
  Arg [CHECK_LENGTH]    : (optional) Boolean - determines if projected objects
                          have to have same length as original (default: false)
  Example     : my $ap = Bio::EnsEMBL::Utils::AssemblyProjector->new(
                  -DBADAPTOR    => $dba,
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

  my ($adaptor, $external_source, $old_assembly, $new_assembly,
      $merge_fragments, $check_length) = rearrange([qw(ADAPTOR EXTERNAL_SOURCE 
        OLD_ASSEMBLY NEW_ASSEMBLY MERGE_FRAGMENTS CHECK_LENGTH)], @_);

  unless ($adaptor and ref($adaptor) and
          $adaptor->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    throw("You must provide a DBAdaptor to a database containing the assembly mapping.");
  }

  unless ($old_assembly and $new_assembly) {
    throw("You must provide an old and new assembly name.");
  }

  my $self = {};
  bless ($self, $class);

  # initialise
  $self->adaptor($adaptor);
  $self->{'old_assembly'} = $old_assembly;
  $self->{'new_assembly'} = $new_assembly;
  
  # by default, merge fragments
  $self->{'merge_fragments'} = $merge_fragments || 1;

  # by default, do not check length
  $self->{'check_length'} = $check_length || 0;

  # by default, features and slices are expected in same database as the 
  # assembly mapping
  $self->{'external_source'} = $external_source || 0;

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
                will have to have the same length as the original. You can use
                the last_status() method to find out what the result of some of
                these rule tests were. Please see the comments in the code for
                more details about these rules.

                The return value of this method will always be a single object,
                or undef if the projection fails any of the rules.
                
                Note that when projecting features, only a "shallow" projection
                is performed, i.e. attached features aren't projected
                automatically! (e.g. if you project a gene, its transcripts will
                have to be projected manually before storing the new gene)
  Return type : same a Arg 1, or undef if projection fails any of the rules
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
    $object_type = 'feature';
  } elsif ($object->isa('Bio::EnsEMBL::Slice')) {
    $object_type = 'slice';
  } else {
    throw("Need a Feature or Slice to project.");
  }

  # if the feature or slice is sourced from another db, we have to "transfer"
  # it to the db that contains the assembly mapping. the transfer is very 
  # shallow but that should do for our purposes
  if ($self->external_source) {
    my $slice_adaptor = $self->adaptor->get_SliceAdaptor;

    if ($object_type eq 'feature') {
    
      # createa a new slice from the target db
      my $f_slice = $object->slice;
      my $target_slice = $slice_adaptor->fetch_by_name($f_slice->name);
      
      # now change the feature so that it appears it's from the target db
      $object->slice($target_slice);
    
    } else {
    
      # createa a new slice from the target db
      $object = $slice_adaptor->fetch_by_name($object->name);
      
    }
  }

  if ($object_type eq 'feature') {
    $slice = $object->feature_Slice;
  } else {
    $slice = $object;
  }

  # warn if trying to project to assembly version the object already is on
  if ($slice->coord_system->version eq $to_assembly) {
    warning("Assembly version to project to ($to_assembly) is the same as your object's assembly (".$slice->coord_system->version.").");
  }

  # now project the slice
  my $cs_name = $slice->coord_system_name;
  my @segments = @{ $slice->project($cs_name, $to_assembly) };

  # we need to reverse the projection segment list if the orignial 
  if ($slice->strand == -1) {
    @segments = reverse(@segments);
  }

  # apply rules to projection results
  # 
  # discard the projected feature/slice if
  #   1. it doesn't project at all (no segments returned)
  #   2. [unless MERGE_FRAGMENTS is set] the projection is fragmented (more
  #      than one segment)
  #   3. [if CHECK_LENGTH is set] the projection doesn't have the same length
  #      as the original feature/slice
  #   4. all segments are on same chromosome and strand

  # keep track of the status of applied rules
  my @status = ();

  # test (1)
  return undef unless (@segments);
  #warn "DEBUG: passed test 1\n";

  # test (2)
  return undef if (!($self->merge_fragments) and scalar(@segments) > 1);
  push @status, 'fragmented' if (scalar(@segments) > 1);
  #warn "DEBUG: passed test 2\n";

  # test (3)
  my $first_slice = $segments[0]->to_Slice;
  my $last_slice = $segments[-1]->to_Slice;
  my $length_mismatch = (($last_slice->end - $first_slice->start + 1) !=
    $object->length);
  return undef if ($self->check_length and $length_mismatch);
  push @status, 'length_mismatch' if ($length_mismatch);
  #warn "DEBUG: passed test 3\n";

  # test (4)
  my %sr_names = ();
  my %strands = ();
  foreach my $seg (@segments) {
    my $sl = $seg->to_Slice;
    $sr_names{$sl->seq_region_name}++;
    $strands{$sl->strand}++;
  }
  return undef if (scalar(keys %sr_names) > 1 or scalar(keys %strands) > 1);
  #warn "DEBUG: passed test 4\n";

  # remember rule status
  $self->last_status(join('|', @status));

  # everything looks fine, so adjust the coords of your feature/slice
  my $new_slice = $first_slice;
  $new_slice->{'end'} = $last_slice->end;

  if ($object_type eq 'slice') {
    return $new_slice;
  } else {
  
    $object->start($new_slice->start);
    $object->end($new_slice->end);
    $object->strand($new_slice->strand);
    $object->slice($new_slice->seq_region_Slice);

    # undef dbID and adaptor so you can store the feature in the target db
    $object->dbID(undef);
    $object->adaptor(undef);

    return $object;
  }
  
}


=head2 old_to_new

  Arg[1]      : Bio::EnsEMBL::Slice or Bio::EnsEMBL::Feature $object -
                the object to project
  Example     : my $new_slice = $assembly_projector->old_to_new($old_slice);
  Description : Projects a Slice or Feature from old to new assembly.
                This method is just a convenience wrapper for $self->project.
  Return type : same a Arg 1, or undef
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
  Return type : same a Arg 1, or undef
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
sub adaptor {
  my $self = shift;
  weaken($self->{'adaptor'} = shift) if (@_);
  return $self->{'adaptor'};
}


sub external_source {
  my $self = shift;
  $self->{'external_source'} = shift if (@_);
  return $self->{'external_source'};
}


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


sub last_status {
  my $self = shift;
  $self->{'last_status'} = shift if (@_);
  return $self->{'last_status'};
}


1;


