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

Bio::EnsEMBL::CoordSystem

=head1 SYNOPSIS

  my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(...);

  my $csa = $db->get_CoordSystemAdaptor();

  #
  # Get all coord systems in the database:
  #
  foreach my $cs ( @{ $csa->fetch_all() } ) {
    my $str = join ':', $cs->name(), $cs->version(), $cs->dbID();
    print "$str\n";
  }

=head1 DESCRIPTION

This is a simple object which contains a few coordinate system attributes:
name, internal identifier, version.  A coordinate system is uniquely defined
by its name and version.  A version of a coordinate system applies to all
sequences within a coordinate system.  This should not be confused with
individual sequence versions.

Take for example the Human assembly.  The version 'NCBI33' applies to
to all chromosomes in the NCBI33 assembly (that is the entire 'chromosome'
coordinate system).  The 'clone' coordinate system in the same database would
have no version however.  Although the clone sequences have their own sequence
versions, there is no version which applies to the entire set of clones.

Coordinate system objects are immutable. Their name and version, and other
attributes may not be altered after they are created.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::CoordSystem;

use Bio::EnsEMBL::Storable;

use Bio::EnsEMBL::Utils::Argument  qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

  Arg [..]   : List of named arguments:
               -NAME      - The name of the coordinate system
               -VERSION   - (optional) The version of the coordinate system.
                            Note that if the version passed in is undefined,
                            it will be set to the empty string in the
                            resulting CoordSystem object.
               -RANK      - The rank of the coordinate system. The highest
                            level coordinate system should have rank 1, the
                            second highest rank 2 and so on.  An example of
                            a high level coordinate system is 'chromosome' an
                            example of a lower level coordinate system is
                            'clone'.
               -TOP_LEVEL - (optional) Sets whether this is a top-level coord
                            system. Default = 0. This should only be set to
                            true if you are creating an artificial toplevel
                            coordsystem by the name of 'toplevel'
               -SEQUENCE_LEVEL - (optional) Sets whether this is a sequence
                            level coordinate system. Default = 0
               -DEFAULT   - (optional)
                            Whether this is the default version of the 
                            coordinate systems of this name. Default = 0
               -DBID      - (optional) The internal identifier of this
                             coordinate system
               -ADAPTOR   - (optional) The adaptor which provides database
                            interaction for this object
  Example    : $cs = Bio::EnsEMBL::CoordSystem->new(-NAME    => 'chromosome',
                                                    -VERSION => 'NCBI33',
                                                    -RANK    => 1,
                                                    -DBID    => 1,
                                                    -ADAPTOR => adaptor,
                                                    -DEFAULT => 1,
                                                    -SEQUENCE_LEVEL => 0);
  Description: Creates a new CoordSystem object representing a coordinate
               system.
  Returntype : Bio::EnsEMBL::CoordSystem
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my ( $name, $version, $top_level, $sequence_level, $default, $rank ) =
    rearrange( [ 'NAME',      'VERSION',
                 'TOP_LEVEL', 'SEQUENCE_LEVEL',
                 'DEFAULT',   'RANK' ],
               @_ );

  $top_level      = ($top_level)      ? 1 : 0;
  $sequence_level = ($sequence_level) ? 1 : 0;
  $default        = ($default)        ? 1 : 0;
  $rank ||= 0;

  if ( $top_level == 1 ) {
    if ( $rank != 0 ) {
      throw('RANK argument must be 0 if TOP_LEVEL is 1');
    }

    if ( defined($name) ) {
      if ( $name ne 'toplevel' ) {
        throw('The NAME argument must be "toplevel" if TOP_LEVEL is 1');
      }
    } else {
      $name = 'toplevel';
    }

    if ( $sequence_level == 1 ) {
      throw("SEQUENCE_LEVEL argument must be 0 if TOP_LEVEL is 1");
    }

    $default = 0;

  } else {

    if ( $rank == 0 ) {
      throw("RANK argument must be non-zero unless TOP_LEVEL is 1");
    }

    if ( !defined($name) ) {
      throw('The NAME argument is required');
    } elsif ( $name eq 'toplevel' ) {
      throw(   "Cannot name coord system 'toplevel' "
             . "unless TOP_LEVEL is 1" );
    }

  }

  if ( $rank !~ /^\d+$/ ) {
    throw('The RANK argument must be a positive integer');
  }

  $self->{'version'}        = $version;
  $self->{'name'}           = $name;
  $self->{'top_level'}      = $top_level;
  $self->{'sequence_level'} = $sequence_level;
  $self->{'default'}        = $default;
  $self->{'rank'}           = $rank;

  return $self;
} ## end sub new


=head2 name

  Arg [1]    : (optional) string $name
  Example    : print $coord_system->name();
  Description: Getter for the name of this coordinate system
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub name {
  my $self = shift;
  return $self->{'name'};
}



=head2 version

  Arg [1]    : none
  Example    : print $coord->version();
  Description: Getter for the version of this coordinate system.  This
               will return an empty string if no version is defined for this
               coordinate system.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub version {
  my $self = shift;

  return '' unless defined $self->{'version'};
  return $self->{'version'};
}



=head2 species

  Arg [1]    : none
  Example    : print $coord->species();
  Description: Shortcut method to get the species this CoordSystem refers to.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub species {
  my $self = shift;
  return $self->adaptor->db->species;
}



=head2 equals

  Arg [1]    : Bio::EnsEMBL::CoordSystem $cs
               The coord system to compare to for equality.
  Example    : if($coord_sys->equals($other_coord_sys)) { ... }
  Description: Compares 2 coordinate systems and returns true if they are
               equivalent.  The definition of equivalent is sharing the same
               name and version.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub equals {
  my $self = shift;
  my $cs = shift;

  if(!$cs || !ref($cs) || !$cs->isa('Bio::EnsEMBL::CoordSystem')) {
    if ($cs->isa('Bio::EnsEMBL::ExternalData::DAS::CoordSystem')) {
      return $cs->equals($self);
    }
    throw('Argument must be a CoordSystem');
  }

  if($self->version() eq $cs->version() && $self->{'name'} eq $cs->name()) {
    return 1;
  }

  return 0;
}




=head2 is_top_level

  Arg [1]    : none
  Example    : if($coord_sys->is_top_level()) { ... }
  Description: Returns true if this is the toplevel pseudo coordinate system.
               The toplevel coordinate system is not a real coordinate system
               which is stored in the database, but it is a placeholder that
               can be used to request transformations or retrievals to/from
               the highest defined coordinate system in a given region.
  Returntype : 1 or 0
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub is_top_level {
  my $self = shift;
  return $self->{'top_level'};
}


=head2 is_sequence_level

  Arg [1]    : none
  Example    : if($coord_sys->is_sequence_level()) { ... }
  Description: Returns true if this is a sequence level coordinate system
  Returntype : 1 or 0
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub is_sequence_level {
  my $self = shift;
  return $self->{'sequence_level'};
}


=head2 is_default

  Arg [1]    : none
  Example    : if($coord_sys->is_default()) { ... }
  Description: Returns true if this coordinate system is the default
               version of the coordinate system of this name.
  Returntype : 1 or 0
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub is_default {
  my $self = shift;
  return $self->{'default'};
}




=head2 rank

  Arg [1]    : none
  Example    : if($cs1->rank() < $cs2->rank()) {
                 print $cs1->name(), " is a higher level coord system than",
                       $cs2->name(), "\n";
               }
  Description: Returns the rank of this coordinate system.  A lower number
               is a higher coordinate system.  The highest level coordinate
               system has a rank of 1 (e.g. 'chromosome').  The toplevel
               pseudo coordinate system has a rank of 0.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub rank {
  my $self = shift;
  return $self->{'rank'};
}

1;
