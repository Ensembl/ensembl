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
    print $cs->name, ' ',  $cs->version, "\n";
  }

  #
  # Fetching by name:
  #

  #use the default version of coord_system 'chromosome' (e.g. NCBI33):
  $cs = $csa->fetch_by_name('chromosome');

  #get an explicit version of coord_system 'chromosome':
  $cs = $csa->fetch_by_name('chromsome', 'NCBI34');

  #get all coord_systems of name 'chromosome':
  foreach $cs (@{$csa->fetch_all_by_name('chromosome')}) {
     print $cs->name, ' ', $cs->version, "\n";
  }

  #
  # Fetching by top level:
  #

  #Get the default top_level coord system:
  $cs = $csa->fetch_top_level();

  #Get a particular version of a top_level coord system:
  $cs = $csa->fetch_top_level('NCBI34');

  #Get all top level coord systems:
  foreach $cs (@{$csa->fetch_all_top_level()}) {
    print $cs->name(), ' ', $cs->version, "\n";
  }

  #can also use an alias in fetch_by_name:
  $cs = $csa->fetch_by_name('toplevel');

  #
  # Fetching by sequence level:
  #

  #Get the coord system which is used to store sequence:
  $cs = $csa->fetch_sequence_level();

  #can also use an alias in fetch_by_name:
  $cs = $csa->fetch_by_name('seqlevel');

  #
  # Fetching by id
  #
  $cs = $csa->fetch_by_dbID(1);


=head1 DESCRIPTION

This adaptor allows the querying of information from the coordinate system
adaptor.

Note that many coordinate systems do not have a concept of a version
for the entire coordinate system (though they may have a per-sequence version).
The 'chromosome' coordinate system usually has a version (i.e. the
assembly version) but the clonal coordinate system does not (despite having
individual sequence versions).  In the case where a coordinate system does
not have a version an empty string ('') is used instead.

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
use Bio::EnsEMBL::CoordSystem;

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

  my $self = $class->SUPER::new(@_);

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
    my $seq_lvl = 0;
    my $top_lvl = 0;
    foreach my $attrib (split(',', $attrib)) {
      $self->{"_is_$attrib"}->{$dbID} = 1;
      if($attrib eq 'sequence_level') {
        $seq_lvl = 1;
      } elsif($attrib eq 'top_level') {
        $top_lvl = 1;
      }
    }

    my $cs = Bio::EnsEMBL::CoordSystem->new
      (-DBID           => $dbID,
       -ADAPTOR        => $self,
       -NAME           => $name,
       -VERSION        => $version,
       -TOP_LEVEL      => $top_lvl,
       -SEQUENCE_LEVEL => $seq_lvl);

    $self->{'_dbID_cache'}->{$dbID} = $cs;

    $self->{'_name_cache'}->{lc($name)} ||= [];
    push @{$self->{'_name_cache'}->{lc($name)}}, $cs;
  }
  $sth->finish();

  #
  # Retrieve the list of the coordinate systems that features are stored in
  # and cache them
  #
  $sth = $self->prepare('SELECT table_name, coord_system_id FROM meta_coord');
  $sth->execute();
  while(my ($table_name, $dbID) = $sth->fetchrow_array()) {
    $self->{'_feature_cache'}->{lc($table_name)} ||= [];
    my $cs = $self->{'_dbID_cache'}->{$dbID};
    if(!$cs) {
      throw("meta_coord table refers to non-existant coord_system id=[$dbID]");
    }
    push @{$self->{'_feature_cache'}->{lc($table_name)}}, $cs;
  }
  $sth->finish();

  #
  # Retrieve a list of available mappings from the meta table.
  # this may eventually be moved a table of its own if this proves too
  # cumbersome
  #

  my %mappings;
  my $mc = $self->db()->get_MetaContainer();
  foreach my $map_pair (@{$mc->list_value_by_key('assembly.mapping')}) {
    my ($asm,$cmp) = split(/\|/, $map_pair);
    if(!$cmp || !$cmp) {
      throw('incorrectly formatted assembly.mapping values in meta table');
    }
    my($asm_name, $asm_version) = split(/:/, $asm);
    my($cmp_name, $cmp_version) = split(/:/, $cmp);

    my $cmp_cs = $self->fetch_by_name($cmp_name,$cmp_version);
    my $asm_cs = $self->fetch_by_name($asm_name,$asm_version);

    $mappings{$asm_cs->dbID} ||= {};
    $mappings{$asm_cs->dbID}->{$cmp_cs->dbID} = 1;
  }

  $self->{'_mapping_paths'} = \%mappings;

  return $self;
}



=head2 fetch_all

  Arg [1]    : none
  Example    : foreach my $cs (@{$csa->fetch_all()}) {
                 print $cs->name(), ' ', $cs->version(), "\n";
               }
  Description: Retrieves every coordinate system defined in the DB
  Returntype : listref of Bio::EnsEMBL::CoordSystems
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
               The name of the coordinate system to retrieve.  Alternatively
               this may be an alias for a real coordinate system.  Valid
               aliases are 'toplevel' and 'seqlevel'.
  Arg [2]    : string $version (optional)
               The version of the coordinate system to retrieve.  If not
               specified the default version will be used.
  Example    : $coord_sys = $csa->fetch_by_name('clone');
               $coord_sys = $csa->fetch_by_name('chromosome', 'NCBI33');
               #toplevel is an alias for the highest coordinate system
               #such as the chromosome coordinate system
               $coord_sys = $csa->fetch_by_name('toplevel');
               $coord_sys = $csa->fetch_by_name('toplevel', 'NCBI31');
               #seqlevel is an alias for the sequence level coordinate system
               #such as the clone or contig coordinate system
               $coord_sys = $csa->fetch_by_name('seqlevel');
  Description: Retrieves a coordinate system by its name
  Returntype : Bio::EnsEMBL::CoordSystem
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


  if($name eq 'seqlevel') {
    return $self->fetch_sequence_level();
  } elsif($name eq 'toplevel') {
    return $self->fetch_top_level($version);
  }

  if(!exists($self->{'_name_cache'}->{$name})) {
    my $guess = '';
    if($name =~ /top/) {
      $guess = "\nDid you mean 'toplevel' instead of '$name'?";
    } elsif($name =~ /seq/) {
      $guess = "\nDid you mean 'seqlevel' instead of '$name'?";
    }
    throw("Coord_system with name [$name] does not exist.$guess");
  }

  my @coord_systems = @{$self->{'_name_cache'}->{$name}};

  throw("Coord_system with name [$name] does not exist.") if(!@coord_systems);

  foreach my $cs (@coord_systems) {
    if($version) {
      return $cs if(lc($cs->version()) eq $version);
    } elsif($self->{'_is_default_version'}->{$cs->dbID()}) {
      return $cs;
    }
  }

  if($version) {
    throw("Coord_system [$name] does not exist with version [$version]");
  }

  #didn't find a default, just take first one
  my $cs =  shift @coord_systems;
  my $v = $cs->version();
  warning("No default version for coord_system [$name] exists. " .
      "Using version [$v] arbitrarily");

  return $cs;
}


=head2 fetch_all_by_name

  Arg [1]    : string $name
               The name of the coordinate system to retrieve.  This can be
               the name of an actual coordinate system or an alias for a
               coordinate system.  Valid aliases are 'toplevel' and 'seqlevel'.
  Example    : foreach my $cs (@{$csa->fetch_all_by_name('chromosome')}){
                 print $cs->name(), ' ', $cs->version();
               }
  Description: Retrieves all coordinate systems of a particular name
  Returntype : listref of Bio::EnsEMBL::CoordSystem objects
  Exceptions : throw if no name argument provided
  Caller     : general

=cut

sub fetch_all_by_name {
  my $self = shift;
  my $name = lc(shift); #case insensitive matching

  throw('Name argument is required') if(!$name);

  if($name eq 'seqlevel') {
    return [$self->fetch_sequence_level()];
  } elsif($name eq 'toplevel') {
    return $self->fetch_all_top_level();
  }

  return $self->{'_name_cache'}->{$name} || [];
}




=head2 fetch_all_by_feature_table

  Arg [1]    : string $table - the name of the table to retrieve coord systems
               for
  Example    : my @coord_systems = $csa->fetch_by_feature_table('gene')
  Description: This retrieves the list of coordinate systems that features
               in a particular table are stored.  It is used internally by
               the API to perform queries to these tables and to ensure that
               features are only stored in appropriate coordinate systems.
  Returntype : listref of Bio::EnsEMBL::CoordSystem objects
  Exceptions : none
  Caller     : BaseFeatureAdaptor

=cut

sub fetch_all_by_feature_table {
  my $self = shift;
  my $table = lc(shift); #case insensitive matching

  throw('Name argument is required') unless $table;

  my $result = $self->{'_feature_cache'}->{$table};

  if(!$result) {
    throw("Feature table [$table] does not have a defined coordinate system" .
          " in the meta_coord table");
  }

  return $result;
}


=head2 add_feature_table

  Arg [1]    : Bio::EnsEMBL::CoordSystem $cs
               The coordinate system to associate with a feature table
  Arg [2]    : string $table - the name of the table in which features of
               a given coordinate system will be stored in
  Example    : $csa->add_feature_table($chr_coord_system, 'gene');
  Description: This function tells the coordinate system adaptor that
               features from a specified table will be stored in a certain
               coordinate system.  If this information is not already stored
               in the database it will be added.
  Returntype : none
  Exceptions : none
  Caller     : BaseFeatureAdaptor

=cut


sub add_feature_table {
  my $self = shift;
  my $cs   = shift;
  my $table = lc(shift);

  if(!ref($cs) || !$cs->isa('Bio::EnsEMBL::CoordSystem')) {
    throw('CoordSystem argument is required.');
  }

  if(!$table) {
    throw('Table argument is required.');
  }

  my $coord_systems = $self->{'_feature_cache'}->{$table};

  my ($exists) = grep {$_->equals($cs)} @$coord_systems;

  #do not do anything if this feature table is already associated with the
  #coord system
  return if($exists);

  #make sure to use a coord system stored in this database so that we
  #do not use the wrong coord_system_id
  if(!$cs->is_stored($self->db())) {
    $cs = $self->fetch_by_name($cs->name(), $cs->version);
    throw("CoordSystem not found in database.");
  }

  #store the new tablename -> coord system relationship in the db
  #ignore failures b/c during the pipeline multiple processes may try
  #to update this table and only the first will be successful
  my $sth = $self->prepare('INSERT IGNORE INTO meta_coord ' .
                              'SET coord_system_id = ?, ' .
                                  'table_name = ?');

  $sth->execute($cs->dbID, $table);

  #update the internal cache
  $self->{'_feature_cache'}->{$table} ||= [];
  push @{$self->{'_feature_cache'}->{$table}}, $cs;

  return;
}



=head2 fetch_by_dbID

  Arg [1]    : int dbID
  Example    : $cs = $csa->fetch_by_dbID(4);
  Description: Retrieves a coord_system via its internal
               identifier.
  Returntype : Bio::EnsEMBL::CoordSystem
  Exceptions : thrown if no coord_system exists for specified dbID
  Caller     : general

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  throw('dbID argument is required') if(!$dbID);

  my $cs = $self->{'_dbID_cache'}->{$dbID};

  throw("Coord_system with dbID [$dbID] does not exist")  if(!$cs);

  return $cs;
}



=head2 fetch_top_level

  Arg [1]    : string $version (optional)
               The version of the top-level to obtain. If not provided
               the default version top-level will be returned.
  Example    : $cs = $csa->fetch_top_level();
               $cs = $csa->fetch_top_level('NCBI34');
  Description: Retrieves the top level coordinate system.  It is possible
               to have multiple top-level coordinate systems if there is
               more than one version of the top-level system.  For
               example the top level of a homo_sapiens database could
               be 'chromosome' but this might have 'NCBI33', 'NCBI31',
               and 'NCBI34' versions. One of these will be defined as
               the default and will be returned if no specific version
               is requested. If a specific version is requested then it
               will be returned.
  Returntype : a Bio::EnsEMBL::CoordSystem object
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
  Returntype : listref of Bio::EnsEMBL::CoordSystem objects
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
  Returntype : Bio::EnsEMBL::CoordSystem
  Exceptions : throw if no sequence_level coord system exists at all
               throw if multiple sequence_level coord systems exists
  Caller     : general

=cut

sub fetch_sequence_level {
  my $self = shift;

  my @dbIDs = keys %{$self->{'_is_sequence_level'}};

  throw('No sequence_level coord_system is defined') if(!@dbIDs);

  if(@dbIDs > 1) {
    throw('Multiple sequence_level coord_systems are defined.' .
          'Only one is currently supported');
  }

  return $self->{'_dbID_cache'}->{$dbIDs[0]};
}




=head2 get_mapping_path

  Arg [1]    : Bio::EnsEMBL::CoordSystem $cs1
  Arg [2]    : Bio::EnsEMBL::CoordSystem $cs2
  Example    : foreach my $cs @{$csa->get_mapping_path($cs1,$cs2);
  Description: Given two coordinate systems this will return a mapping path
               between them.  The path is formatted as a list of coordinate
               systems starting with the assembled coord systems and
               descending through component systems.  For example, if the
               following mappings were defined in the meta table:
               chromosome -> clone
               clone -> contig

               And the contig and chromosome coordinate systems where
               provided as arguments like so:
               $csa->get_mapping_path($chr_cs,$ctg_cs);

               The return values would be:
               [$chr_cs, $clone_cs, $contig_cs]

               The return value would be the same even if the order of
               arguments was reversed.

	       This becomes a bit more problematic when the relationship is
               something like:
               chromosome -> contig
               clone      -> contig

               In this case the contig coordinate system is the component
               coordinate system for both mappings and for the following
               request:
               $csa->get_mappging_path($chr_cs, $cln_cs);

	       Either of the following mapping paths would be valid:
               [$chr_cs, $contig_cs, $clone_cs]
               or
               [$clone_cs, $contig_cs, $chr_cs]

	       Also note that the ordering of the above is not
               assembled to component but rather
               assembled -> component -> assembled.

               If no mapping path exists, an reference to an empty list is
               returned.

  Returntype : listref of coord_sytem ids ordered from assembled to smaller
               component coord_systems
  Exceptions : none
  Caller     : general

=cut

sub get_mapping_path {
  my $self = shift;
  my $cs1 = shift;
  my $cs2 = shift;
  my $seen = shift || {};

  $self->{'_shortest_path'} ||= {};

  my $cs1_id = $cs1->dbID();
  my $cs2_id = $cs2->dbID();

  # if this method has already been called with the same arguments
  # return the cached result
  if($self->{'_shortest_path'}->{"$cs1_id:$cs2_id"}) {
    return $self->{'_shortest_path'}->{"$cs1_id:$cs2_id"};
  }

  #if we have already seen this pair then there is some circular logic
  #encoded in the mappings.  This is not good.
  if($seen->{"$cs1_id:$cs2_id"}) {
    throw("Circular logic detected in defined assembly mappings");
  }

  #if there is a direct mapping between this coord system and other one
  #then path between cannot be shorter, just return the one step path
  if($self->{'_mapping_paths'}->{$cs1_id}->{$cs2_id}) {
    $self->{'_shortest_path'}->{"$cs1_id:$cs2_id"} = [$cs1,$cs2];
    $self->{'_shortest_path'}->{"$cs2_id:$cs1_id"} = [$cs1,$cs2];
    return [$cs1,$cs2];
  }
  if($self->{'_mapping_paths'}->{$cs2_id}->{$cs1_id}) {
    $self->{'_shortest_path'}->{"$cs1_id:$cs2_id"} = [$cs2,$cs1];
    $self->{'_shortest_path'}->{"$cs2_id:$cs1_id"} = [$cs2,$cs1];
    return [$cs2,$cs1];
  }

  $seen->{"$cs1_id:$cs2_id"} = 1;
  $seen->{"$cs2_id:$cs1_id"} = 1;

  # There is no direct mapping available, but there may be an indirect
  # path.  Call this method recursively on the components of both paths.
  my $shortest;

  #need to try both as assembled since we do not know which is the assembled
  #coord_system and which is the component
  foreach my $pair ([$cs1,$cs2], [$cs2,$cs1]) {
    my ($asm_cs, $target_cs) = @$pair;
    my $asm_cs_id = $asm_cs->dbID();

    foreach my $cmp_cs_id (keys %{$self->{'_mapping_paths'}->{$asm_cs_id}}) {
      my $cmp_cs = $self->fetch_by_dbID($cmp_cs_id);
      my $path = $self->get_mapping_path($cmp_cs, $target_cs, $seen);
      my $len = @$path;
      my $shortest;

      next if($len == 0);

      #Check whether the component was used as an assembled
      #or component in the next part of the path:
      if($cmp_cs_id == $path->[0]->dbID) {
        $path = [$asm_cs, @$path];
      } else {
        $path = [@$path, $asm_cs];
      }

      #shortest possible indirect, add to path so far and return
      if($len == 2) {
        $self->{'_shortest_path'}->{"$cs1_id:$cs2_id"} = $path;
        $self->{'_shortest_path'}->{"$cs2_id:$cs1_id"} = $path;
        return $path;
      } elsif(!defined($shortest) || $len+1 < @$shortest) {
        #keep track of the shortest path found so far,
        #there may yet be shorter..
        $shortest = $path;
      }
    }
    #use the shortest path found so far,
    #if no path was found continue, using the the other id as assembled
    if(defined($shortest)) {
      $self->{'_shortest_path'}->{"$cs1_id:$cs2_id"} = $shortest;
      $self->{'_shortest_path'}->{"$cs2_id:$cs1_id"} = $shortest;
      return $shortest;
    }
  }

  #did not find any possible path
  $self->{'_shortest_path'}->{"$cs1_id:$cs2_id"} = [];
  $self->{'_shortest_path'}->{"$cs2_id:$cs1_id"} = [];
  return [];
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
      return $cs if(lc($version) eq $cs->version());
    } elsif($self->{'_is_default_version'}->{$dbID}) {
      return $cs;
    }
  }

  #specifically requested attrib system was not found
  if($version) {
    throw("$attrib coord_system with version [$version] does not exist");
  }

  #coordsystem with attrib exists but no default is defined:
  my $dbID = shift @dbIDs;
  my $cs = $self->{'_dbID_cache'}->{$dbID};
  my $v = $cs->version();
  warning("No default version for $attrib coord_system exists. " .
          "Using version [$v] arbitrarily");

  return $cs;
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


#
# Called during db destruction to clean up internal cache structures etc.
#
sub deleteObj {
  my $self = shift;

  #break circular adaptor <-> db references
  $self->SUPER::deleteObj();

  #breack circular object <-> adaptor references
  delete $self->{'_feature_cache'};
  delete $self->{'_name_cache'};
  delete $self->{'_dbID_cache'};
  delete $self->{'_mapping_paths'};
  delete $self->{'_shortest_paths'};
}

1;




