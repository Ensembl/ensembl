#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::CoordSystemAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::DBSQL::CoordSystemAdaptor

=head1 SYNOPSIS

  my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(...);

  my $csa = $db->get_CoordSystemAdaptor();

  #
  # Get all coord systems in the database:
  #
  foreach my $cs (@{$csa->fetch_all()}) {
    my ($id, $name, $version) = @$cs;
    print "$name $version\n";
  }

  #
  # Fetching by name:
  #

  #use the default version of coord_system 'chromosome' (e.g. NCBI33):
  my ($id, $name, $version) = $csa->fetch_by_name('chromosome');

  #get an explicit version of coord_system 'chromosome':
  my ($id, $name, $version) = $csa->fetch_by_name('chromsome', 'NCBI34');

  #get all coord_systems of name 'chromosome':
  foreach my $cs (@{$csa->fetch_all_by_name('chromosome')}) {
    my ($id, $name, $version) = @$cs;
  }

  #
  # Fetching by top level:
  #

  #Get the default top_level coord system:
  my ($id, $name, $version) = $csa->fetch_top_level();

  #Get a particular version of a top_level coord system:
  my ($id, $name, $version) = $csa->fetch_top_level('NCBI34');

  #Get all top level coord systems:
  foreach my $cs (@{$csa->fetch_all_top_level()}) {
    my ($id, $name, $version) = @$cs;
  }

  #
  # Fetching by sequence level:
  #

  #Get the coord system which is used to store sequence:
  my ($id, $name, $version) = $csa->fetch_sequence_level();

  #
  # Fetching by id
  #
  my ($id, $name, $version) = $csa->fetch_by_dbID(1);


=head1 DESCRIPTION

This adaptor allows the querying of information from the coordinate system
adaptor.  There is no CoordSystem object since this is really a source
of meta information and the needed information can be easily
(and more speedily) represented as (id,name,version) triplets.

Note that many coordinate systems do not have a concept of a version
for the entire coordinate system (though they may have a per-sequence version).The 'chromosome' coordinate system usually has a version (i.e. the 
assembly version) but the clonal coordinate system does not (despite having
individual sequence versions).  In the case where a coordinate system does
not have a version a triplet is still used for representation, but the version
is an empty string ''.


=head1 AUTHOR - Graham McVicker

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

use strict;
use warnings;

package Bio::EnsEMBL::DBSQL::CoordSystemAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);



=head2 new

  Arg [1]    : See BaseAdaptor for arguments (none specific to this
               subclass)
  Example    : $cs = $db->get_CoordSystemAdaptor(); #better than new()
  Description: Creates a new CoordSystem adaptor and caches the contents
               of the coord_system table in memory.
  Returntype : Bio::EnsEMBL::DBSQL::CoordSystemAdaptor
  Exceptions : none
  Caller     :

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new();

  #
  # Cache the entire contents of the coord_system table cross-referenced
  # by dbID and name
  #

  #keyed on name, list of coord_system value
  $self->{'_name_cache'} = {};

  #keyed on id, coord_system value
  $self->{'_dbID_cache'} = {};

  #keyed on id, 1/undef values
  $self->{'_is_sequence_level'} = {};
  $self->{'_is_top_level'} = {};
  $self->{'_is_default_version'} = {};

  my $sth = $self->prepare('SELECT coord_system_id, name, version, attrib ' .
                           'FROM coord_system');

  $sth->execute();

  my ($dbID, $name, $version, $attrib);
  $sth->bind_columns(\$dbID, \$name, \$version, \$attrib);

  while($sth->fetch()) {
    my $cs = [$dbID, $name, $version || ''];
    my @attribs = split(',',$attrib);

    $self->{'_dbID_cache'}->{$dbID} = $cs;

    $self->{'_name_cache'}->{lc($name)} ||= [];
    push @{$self->{'_name_cache'}->{lc($name)}}, $cs;

    foreach my $attrib (@attribs) {
      $self->{"_is_$attrib"}->{$dbID} = 1;
    }
  }

  return $self;
}



=head2 fetch_all

  Arg [1]    : none
  Example    : foreach my $cs (@{$csa->fetch_all()}) {
                 my ($id, $name, $version) = @$cs;
                 print "$name $version\n";
               }
  Description: Retrieves every coordinate system defined in the DB
               formatted as a listref of [dbID,name,version] triplets.
  Returntype : listref of [dbID,name,version] triplets
  Exceptions : none
  Caller     : general

=cut

sub fetch_all {
  my $self = shift;

  my @coord_systems = values %{$self->{'_dbID_cache'}};

  return \@coord_systems;
}



=head2 fetch_by_name

  Arg [1]    : string $name
               The name of the coordinate system to retrieve
  Arg [2]    : string $version (optional)
               The version of the coordinate system to retrieve.  If not
               specified the default version will be used.
  Example    : ($id, $name, $ver) = $csa->fetch_by_name('clone');
               ($id, $name, $ver) = $csa->fetch_by_name('chromosome',
                                                        'NCBI33');
  Description: Retrieves a coordinate system by its name
  Returntype : (id, name, version) triplet
  Exceptions : throw if the requested coordinate system does not exist
               throw if no name argument provided
               warning if no version provided and default does not exist
  Caller     : general

=cut

sub fetch_by_name {
  my $self = shift;
  my $name = lc(shift); #case insensitve matching
  my $version = shift;

  throw('Name argument is required.') if(!$name);

  $version = lc($version) if($version);

  my @coord_systems = @{$self->{'_name_cache'}->{$name}};

  throw('Coord_system with name [$name] does not exist.') if(!@coord_systems);

  foreach my $cs (@coord_systems) {
    if($version) {
      return @$cs if(lc($cs->[2]) eq $version);
    } elsif($self->{'_is_default'}) {
      return @$cs;
    }
  }

  if($version) {
    throw("Coord_system [$name] does not exist with version [$version]");
  }

  #didn't find a default, just take first one
  my $cs =  shift @coord_systems;
  my $v = $cs->[2];
  warning("No default version for coord_system [$name] exists. " .
      "Using version [$v] arbitrarily");

  return @$cs;
}


=head2 fetch_all_by_name

  Arg [1]    : string $name
               The name of the coordinate system to retrieve
  Example    : foreach my $cs (@{$csa->fetch_all_by_name('chromosome')}){
                 my ($id, $name, $version) = @$cs;
               }
  Description: Retrieves all coordinate systems of a particular name
  Returntype : listref of [id, name, version triplets]
  Exceptions : throw if no name argument provided
  Caller     : general

=cut

sub fetch_all_by_name {
  my $self = shift;
  my $name = lc(shift); #case insensitive matching

  throw('Name argument is required') if(!$name);

  return $self->{'_name_cache'}->{$name} || [];
}



=head2 fetch_by_dbID

  Arg [1]    : int dbID
  Example    : ($dbID, $name, $version) = $csa->fetch_by_dbID(4);
  Description: Retrieve a coord_system via its internal
               identifier.  The coord_system is returned formatted as
               a (id, name, version) triplet.
  Returntype : (id, name, version) triplet
  Exceptions : thrown if no coord_system exists for specified dbID
  Caller     : general

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  throw('dbID argument is required') if(!$dbID);

  my $cs = $self->{'_dbID_cache'}->{'dbID'};

  throw ('Coord_system with dbID [$dbID] does not exist')  if(!$cs);

  return @$cs;
}



=head2 fetch_top_level

  Arg [1]    : string $version (optional)
               The version of the top-level to obtain. If not provided
               the default version top-level will be returned.
  Example    : ($id, $name, $version) = $csa->fetch_top_level();
               ($id, $name, $version) = $csa->fetch_top_level('NCBI34');
  Description: Retrieves the top level coordinate system.  It is possible
               to have multiple top-level coordinate systems if there is
               more than one version of the top-level system.  For
               example the top level of a homo_sapiens database could
               be 'chromosome' but this might have 'NCBI33', 'NCBI31',
               and 'NCBI34' versions. One of these will be defined as
               the default and will be returned if no specific version
               is requested. If a specific version is requested then it
               will be returned.
  Returntype : a (id, name, version) triplet
  Exceptions : throw if no top-level coord_system exists for specified
               version
               throw if no top-level coord_system exists at all
               warning if no version specified and no default version
               exists
  Caller     : general

=cut

sub fetch_top_level {
  my $self = shift;

  return $self->_fetch_by_attrib('top_level', @_);
}


=head2 fetch_all_top_level

  Arg [1]    : none
  Example    : foreach my $cs (@{$csa->fetch_all_top_level()}) {
                 my ($id, $name, $version) = @$cs;
               }
  Description: Retrieves all top level coordinate systems defined in
               the DB. It is possible to have multiple top-level 
               coordinate systems if there is more than one version
               of the top-level system. For example the top level of a
               homo_sapiens database could be 'chromosome' but this
               might have 'NCBI33', 'NCBI31', and 'NCBI34' versions.
  Returntype : listref of [id, name, version] triplets
  Exceptions : none
  Caller     : general

=cut

sub fetch_all_top_level {
  my $self = shift;

  return $self->_fetch_all_by_attrib('top_level', @_);
}



=head2 fetch_sequence_level

  Arg [1]    : none
  Example    : ($id, $name, $version) = $csa->fetch_sequence_level();
  Description: Retrieves the coordinate system at which sequence
               is stored at.
  Returntype : a (id, name, version) triplet
  Exceptions : throw if no sequence_level coord system exists at all
               throw if multiple sequence_level coord systems exists
  Caller     : general

=cut

sub fetch_sequence_level {
  my $self = shift;

  my @dbIDs= keys %{$self->{'_is_sequence_level'}};

  throw('No sequence_level coord_system is defined') if(!@dbIDs);

  if(@dbIDs > 1) {
    throw('Multiple sequence_level coord_systems are defined.' .
          'Only one is currently supported');
  }

  my $dbID = shift;

  return @{$self->{'_dbID_cache'}->{$dbID}};
}



sub _fetch_by_attrib {
  my $self = shift;
  my $attrib = shift;
  my $version = shift;

  $version = lc($version) if($version);

  my @dbIDs = keys %{$self->{"_is_$attrib"}};

  throw("No $attrib coordinate system defined") if(!@dbIDs);

  foreach my $dbID (@dbIDs) {
    my $cs = $self->{'_dbID_cache'}->{$dbID};
    if($version) {
      return @$cs if(lc($version) eq $cs->[2]);
    } elsif($self->{'_is_default'}->{$dbID}) {
      return @$cs;
    }
  }

  #specifically requested attrib system was not found
  if($version) {
    throw("$attrib coord_system with version [$version] does not exist");
  }

  #coordsystem with attrib exists but no default is defined:
  my $dbID = shift @dbIDs;
  my $cs = $self->{'_dbID_cache'}->{$dbID};
  my $v = $cs->[2];
  warning("No default version for $attrib coord_system exists. " .
          "Using version [$v] arbitrarily");

  return @$cs;
}

sub _fetch_all_by_attrib {
  my $self = shift;
  my $attrib = shift;

  my @coord_systems = ();
  foreach my $dbID (keys %{$self->{"_is_$attrib"}}) {
    push @coord_systems, $self->{"_dbID_cache"}->{$dbID};
  }

  return \@coord_systems;
}


1;




