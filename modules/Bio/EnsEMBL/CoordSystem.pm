#
# EnsEMBL module for Bio::EnsEMBL::CoordSystem
#

=head1 NAME

Bio::EnsEMBL::CoordSystem

=head1 SYNOPSIS

  my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(...);

  my $csa = $db->get_CoordSystemAdaptor();

  #
  # Get all coord systems in the database:
  #
  foreach my $cs (@{$csa->fetch_all()}) {
    my $str = join ':', $cs->name(),$cs->version(),$cs->dbID();
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

Coordinate system objects are immutable. Their name and version attributes
may not be altered after they are created.

=head1 AUTHOR - Graham McVicker

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


use strict;
use warnings;

package Bio::EnsEMBL::CoordSystem;

use Bio::EnsEMBL::Storable;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

  Arg [..]   : List of named arguments:
               -NAME    - The name of the coordinate system
               -VERSION - (optional) The version of the coordinate system
               -TOP_LEVEL - (optional) Sets whether this is a top-level coord
                            system. Default = 0
               -SEQUENCE_LEVEL - (optional) Sets whether this is a sequence
                                 level coordinate system. Default = 0
               -DBID    - (optional) The internal identifier of this
                                     coordinate system
               -ADAPTOR - (optional) The adaptor which provides database
                                     interaction for this object
  Example    : $cs = Bio::EnsEMBL::CoordSystem->new(-NAME => 'chromosome',
                                                    -VERSION => 'NCBI33',
                                                    -DBID    => 1,
                                                    -ADAPTOR => adaptor);
  Description: Creates a new CoordSystem object representing a coordinate
               system.
  Returntype : Bio::EnsEMBL::CoordSystem
  Exceptions : none
  Caller     : general

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my ($name,$version, $top_level, $sequence_level) =
    rearrange(['NAME','VERSION','TOP_LEVEL', 'SEQUENCE_LEVEL'], @_);
  
  throw('The NAME argument is required') if(!$name);

  $version = '' if(!defined($version));

  $top_level ||= 0;
  $sequence_level ||= 0;

  $self->{'version'} = $version;
  $self->{'name'} = $name;
  $self->{'top_level'} = $top_level;
  $self->{'sequence_level'} = $sequence_level;

  return $self;
}


=head2 name

  Arg [1]    : (optional) string $name
  Example    : print $coord_system->name();
  Description: Getter for the name of this coordinate system
  Returntype : string
  Exceptions : none
  Caller     : general

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

=cut

sub version {
  my $self = shift;
  return $self->{'version'};
}




=head2 equals

  Arg [1]    : (optional) string version
  Example    : if($coord_sys->equals($other_coord_sys)) { ... }
  Description: Compares 2 coordinate systems
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub equals {
  my $self = shift;
  my $cs = shift;

  if(!$cs || !ref($cs) || !$cs->isa('Bio::EnsEMBL::CoordSystem')) {
    throw('Argument must be a Bio::EnsEMBL::CoordSystem');
  }

  if($self->{'version'} eq $cs->version() && $self->{'name'} eq $cs->name()) {
    return 1;
  }

  return 0;
}




=head2 is_top_level

  Arg [1]    : none
  Example    : if($coord_sys->is_top_level()) { ... }
  Description: Returns true if this is a top level coordinate system
  Returntype : 0 or 1
  Exceptions : none
  Caller     : general

=cut

sub is_top_level {
  my $self = shift;
  return $self->{'top_level'};
}


=head2 is_sequence_level

  Arg [1]    : none
  Example    : if($coord_sys->is_sequence_level()) { ... }
  Description: Returns true if this is a sequence level coordinate system
  Returntype : 0 or 1
  Exceptions : none
  Caller     : general

=cut

sub is_sequence_level {
  my $self = shift;
  return $self->{'top_level'};
}


1;
