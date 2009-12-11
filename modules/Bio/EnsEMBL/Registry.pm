=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Registry

=head1 SYNOPSIS

  use Bio::EnsEMBL::Registry;

  my $registry = 'Bio::EnsEMBL::Registry';

  $registry->load_all("configuration_file");

  $gene_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Gene' );


=head1 DESCRIPTION

All Adaptors are stored/registered using this module. This module should
then be used to get the adaptors needed.

The registry can be loaded from a configuration file using the load_all
method.

If a filename is passed to load_all then this is used.  Else if the
enviroment variable ENSEMBL_REGISTRY is set to the name on an existing
configuration file, then this is used.  Else if the file .ensembl_init
in your home directory exist, it is used.

For the Web server ENSEMBL_REGISTRY should be set in SiteDefs.pm.  This
will then be passed on to load_all.


The registry can also be loaded via the method load_registry_from_db
which given a database host will load the latest versions of the Ensembl
databases from it.

The four types of registries are for db adaptors, dba adaptors, dna
adaptors and the standard type.

=head2 db

These are registries for backwards compatibility and enable the
subroutines to add other adaptors to connections.

e.g. get_all_db_adaptors, get_db_adaptor, add_db_adaptor,
remove_db_adaptor are the old DBAdaptor subroutines which are now
redirected to the Registry.

So if before we had

  my $sfa = $self->adaptor()->db()->get_db_adaptor('blast');

We now want to change this to

  my $sfa =
    Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "blast" );


=head2 DBA

These are the stores for the DBAdaptors

The Registry will create all the DBConnections needed now if you set up
the configuration correctly. So instead of the old commands like

  my $db           = Bio::EnsEMBL::DBSQL::DBAdaptor->new(...);
  my $exon_adaptor = $db->get_ExonAdaptor;

we should now have just

  my $exon_adaptor =
    Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "exon" );


=head2 DNA

This is an internal Registry and allows the configuration of a dnadb.
An example here is to set the est database to get its dna data from the
core database.

  ## set the est db to use the core for getting dna data.
  # Bio::EnsEMBL::Utils::ConfigRegistry->dnadb_add( "Homo Sapiens",
  #   "core", "Homo Sapiens", "est" );


=head2 adaptors

This is the registry for all the general types of adaptors like
GeneAdaptor, ExonAdaptor, Slice Adaptor etc.

These are accessed by the get_adaptor subroutine i.e.

  my $exon_adaptor =
    Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "exon" );

=head1 METHODS

=cut

package Bio::EnsEMBL::Registry;

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning );
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::ConfigRegistry;
use DBI;

use vars qw(%registry_register);

my $API_VERSION = 56;

# This is a map from group names to Ensembl DB adaptors.  Used by
# load_all() and reset_DBAdaptor().
my %group2adaptor = (
  'blast'         => 'Bio::EnsEMBL::External::BlastAdaptor',
  'compara'       => 'Bio::EnsEMBL::Compara::DBSQL::DBAdaptor',
  'core'          => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
  'estgene'       => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
  'funcgen'       => 'Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor',
  'haplotype'     => 'Bio::EnsEMBL::ExternalData::Haplotype::DBAdaptor',
  'hive'          => 'Bio::EnsEMBL::Hive::DBSQL::DBAdaptor',
  'lite'          => 'Bio::EnsEMBL::Lite::DBAdaptor',
  'ontology'      => 'Bio::EnsEMBL::DBSQL::OntologyDBAdaptor',
  'otherfeatures' => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
  'pipeline'      => 'Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor',
  'snp'           => 'Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor',
  'variation'     => 'Bio::EnsEMBL::Variation::DBSQL::DBAdaptor',
  'vega'          => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
);


=head2 load_all

 Will load the registry with the configuration file which is
 obtained from the first in the following and in that order.

  1) If an argument is passed to this method, this is used as the
     name of the configuration file to read.

  2) If the enviroment variable ENSEMBL_REGISTRY is set, this is
     used as the name of the configuration file to read.

  3) If the file .ensembl_init exist in the home directory, it is
     used as the configuration file.

  Arg [1]    : (optional) string
               Name of file to load the registry from.

  Arg [2]    : (optional) integer
               If not 0, will print out all information.

  Arg [3]    : (optional) integer
               If not 0, the database connection will not be
               cleared, if 0 or if not set the database connections
               will be cleared (this is the default).

  Arg [4]:     (optional) boolean
               This option will turn off caching for slice features,
               so, every time a set of features is retrieved,
               they will come from the database instead of the
               cache.  This option is only recommended for advanced
               users, specially if you need to store and retrieve
               features.  It might reduce performance when querying
               the database if not used properly.  If in doubt, do
               not use it or ask in ensembl-dev.

  Example    : Bio::EnsEMBL::Registry->load_all();
  Returntype : none
  Exceptions : none
  Status     : Stable

=cut

sub load_all {
    my $class = shift;
    my ( $config_file, $verbose, $no_clear, $no_cache ) = @_;

    $config_file ||= $ENV{ENSEMBL_REGISTRY}
      || $ENV{HOME} . "/.ensembl_init";

    $verbose  ||= 0;
    $no_clear ||= 0;
    $no_cache ||= 0;

    if ( !defined($config_file) ) {
        if ($verbose) {
            print( STDERR
                   "No default registry configuration to load.\n" );
        }
    } elsif ( !-e $config_file ) {
        if ($verbose) {
            printf( STDERR "Configuration file '%s' does not exist. "
                      . "Registry configuration not loaded.\n",
                    $config_file );
        }
    } else {
        if ( defined( $registry_register{'seen'} ) ) {
            if ( !$no_clear ) {
                if ($verbose) {
                    print( STDERR "Clearing previously loaded "
                           . "registry configuration\n" );
                }
                $class->clear();
            }
        }
        $registry_register{'seen'} = 1;

        if ($verbose) {
            printf( STDERR
                      "Loading registry configuration from '%s'.\n",
                    $config_file );
        }

        my $cfg;

        eval { require Config::IniFiles };
        if ($@) {
          # The user does not have the 'Config::IniFiles' module.
          if ($verbose) {
            print( STDERR "No Config::IniFiles module found, "
                   . "assuming this is not an ini-file\n" );
          }
          # If the configuration file *is* an ini-file, we can expect a
          # load of compilation errors from the next eval...
        } else {
          # The user has the 'Config::IniFiles' module installed.  See
          # if this is an ini-file or not...
          $cfg = Config::IniFiles->new( -file => $config_file );
        }

        if ( defined $cfg ) {
		  my %default_adaptor_args = ();

            if ( $cfg->SectionExists('default') ) {
                # The 'default' section is special.  It contain default
                # values that should be implicit to all other section in
                # this configuration file.  Aliases are added if there
                # is also a 'species' setting.

                my $alias = $cfg->val( 'default', 'alias' );
                $cfg->delval( 'default', 'alias' );

                my $species = $cfg->val( 'default', 'species' );

                if ( defined($alias) && defined($species) ) {
                    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
                                     -species => $species,
                                     -alias => [ split( /\n/, $alias ) ]
                    );
                }

                %default_adaptor_args =
                  map { '-' . $_ => $cfg->val( 'default', $_ ) }
                  $cfg->Parameters('default');
            }

            foreach my $section ( $cfg->Sections() ) {
                if ( $section eq 'default' )
                {    # We have already done the 'default' section.
                    next;
                }

                my $group = $cfg->val( $section, 'group' )
                  || $cfg->val( 'default', 'group' );

                if ( !defined($group) ) {
                    printf( STDERR "Key 'group' is undefined "
                              . "for configuration section '%s', "
                              . "skipping this section.\n",
                            $section );
                    next;
                }

                my $adaptor = $group2adaptor{ lc($group) };
                if ( !defined($adaptor) ) {
                    printf( STDERR "Unknown group '%s' "
                              . "for configuration section '%s', "
                              . "skipping this section.\n",
                            $group, $section );
                    next;
                }

                # Handle aliases.  A section must have both an 'alias'
                # setting and a 'species' setting for aliases to be
                # added.  The 'species' setting might be inherited from
                # the 'default' section.

                my $alias = $cfg->val( $section, 'alias' );
                $cfg->delval( $section, 'alias' );

                my $species = $cfg->val( $section, 'species' )
                  || $cfg->val( 'default', 'species' );

                if ( defined($alias) && defined($species) ) {
                    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
                                     -species => $species,
                                     -alias => [ split( /\n/, $alias ) ]
                    );
                }

                # Fill in the adaptor initialization arguments.
                # We trust the user to provide sensible key-value pairs.
                my %adaptor_args = %default_adaptor_args;
                foreach my $parameter ( $cfg->Parameters($section) ) {
                  $adaptor_args{ '-' . $parameter } =
                    $cfg->val( $section, $parameter );

                  # when set, do not use the feature cache in the
                  # different adaptors
                  if ($no_cache) {
                    $adaptor_args{'-no_cache'} = 1;
                  }
                }

                if ($verbose) {
                    printf( "Configuring adaptor '%s' "
                              . "for configuration section '%s'...\n",
                            $adaptor, $section );
                }

                eval "require $adaptor";
                if ($@) { die($@) }

                $adaptor->new(%adaptor_args);

            } ## end foreach my $section ( $cfg->Sections...
        } else {
            # This is probably no ini-file but an old style piece
            # of configuration written in Perl.  We need to try to
            # require() it.

            eval { require($config_file) };
            if ($@) { die($@) }

            # To make the web code avoid doing this again:
            delete $INC{$config_file};
        }
    } ## end else [ if ( !defined($config_file...
} ## end sub load_all

=head2 clear

 Will clear the registry and disconnect from all databases.

  Example    : Bio::EnsEMBL::Registry->clear();
  Returntype : none
  Exceptions : none
  Status     : Stable

=cut

sub clear{
  my ($self);
  
  foreach my $dba (@{$registry_register{'_DBA'}}){
    if($dba->dbc->connected){
      $dba->dbc->db_handle->disconnect();
    }
  }
  %registry_register = ();
}

#
# db adaptors. (for backwards compatibility)
#

=head2 add_db

  Arg [1]    : db (DBAdaptor) to add adaptor to.
  Arg [2]    : name of the name to add the adaptor to in the registry.
  Arg [3]    : The adaptor to be added to the registry.
  Example    : Bio::EnsEMBL::Registry->add_db($db, "lite", $dba);
  Returntype : none
  Exceptions : none
  Status     : At Risk.
             : This is here for backwards compatibility only and may
             : be removed eventually.  Solution is to make sure the
             : db and the adaptor have the same species and the call
             : is then no longer needed.

=cut

sub add_db {
  my ( $class, $db, $name, $adap ) = @_;

  if ( lc( $db->species() ) ne lc( $adap->species ) ) {
    $registry_register{_SPECIES}{ lc( $db->species() ) }
      { lc( $db->group() ) }{'_special'}{ lc($name) } = $adap;
  }
}

=head2 remove_db

  Arg [1]    : db (DBAdaptor) to remove adaptor from.
  Arg [2]    : name to remove the adaptor from in the registry.
  Example    : my $db = Bio::EnsEMBL::Registry->remove_db($db, "lite");
  Returntype : adaptor
  Exceptions : none
  Status     : At Risk.
             : This is here for backwards compatibility only and may
             : be removed eventually.  Solution is to make sure the
             : db and the adaptor have the same species and the call
             : is then no longer needed.

=cut

sub remove_db {
  my ( $class, $db, $name ) = @_;

  my $ret =
    $registry_register{_SPECIES}{ lc( $db->species() ) }
    { lc( $db->group() ) }{'_special'}{ lc($name) };

  $registry_register{_SPECIES}{ lc( $db->species() ) }
    { lc( $db->group() ) }{'_special'}{ lc($name) } = undef;

  return $ret;
}

=head2 get_db

  Arg [1]    : db (DBAdaptor) to get adaptor from.
  Arg [2]    : name to get the adaptor for in the registry.
  Example    : my $db = Bio::EnsEMBL::Registry->get_db("Human", "core", "lite");
  Returntype : adaptor
  Exceptions : none
  Status     : At Risk.
             : This is here for backwards compatibility only and may
             : be removed eventually.  Solution is to make sure the
             : db and the adaptor have the same species then call
             : get_DBAdaptor instead.

=cut

sub get_db {
  my ( $class, $db, $name ) = @_;

  my $ret = Bio::EnsEMBL::Registry->get_DBAdaptor( lc( $db->species ),
    lc($name) );

  if ( defined($ret) ) { return $ret }

  return $registry_register{_SPECIES}{ lc( $db->species() ) }
    { lc( $db->group() ) }{'_special'}{ lc($name) };
}

=head2 get_all_db_adaptors

  Arg [1]    : db (DBAdaptor) to get all the adaptors from.
  Example    : my $db = Bio::EnsEMBL::Registry->get_all_db_adaptors($db);
  Returntype : adaptor
  Exceptions : none
  Status     : At Risk.
             : This is here for backwards compatibility only and
             : may be removed eventually.  Solution is to make
             : sure the dbs all have the same species then call
             : get_all_DBAdaptors(-species => "human");


=cut

sub get_all_db_adaptors {
  my ( $class, $db ) = @_;
  my %ret = ();

  # we now also want to add all the DBAdaptors for the same species.
  # as add_db_adaptor does not add if it is from the same species.

  foreach my $dba ( @{ $registry_register{'_DBA'} } ) {
    if ( lc( $dba->species() ) eq lc( $db->species() ) ) {
      $ret{ $dba->group() } = $dba;
    }
  }

  foreach my $key (
    keys %{
      $registry_register{_SPECIES}
        { $class->get_alias( $db->species() ) }{ lc( $db->group() ) }
        {'_special'} } )
  {
    $ret{$key} =
      $registry_register{_SPECIES}
      { $class->get_alias( $db->species() ) }{ lc( $db->group() ) }
      {'_special'}{$key};
  }

  return \%ret;
} ## end sub get_all_db_adaptors


#
# DBAdaptors
#

=head2 add_DBAdaptor

  Arg [1]    : name of the species to add the adaptor to in the registry.
  Arg [2]    : name of the group to add the adaptor to in the registry.
  Arg [3]    : The DBAaptor to be added to the registry.
  Example    : Bio::EnsEMBL::Registry->add_DBAdaptor("Human", "core", $dba);
  Returntype : none
  Exceptions : none
  caller     : internal
  Status     : Stable

=cut

sub add_DBAdaptor {
  my ( $class, $species, $group, $adap ) = @_;

  if ( !( $class->alias_exists($species) ) ) {
    $class->add_alias( $species, $species );
  }

  $species = $class->get_alias($species);

  $registry_register{_SPECIES}{$species}{ lc($group) }{'_DB'} = $adap;

  if ( !defined( $registry_register{'_DBA'} ) ) {
    my @list = ();
    push( @list, $adap );
    $registry_register{'_DBA'} = \@list;
  } else {
    push( @{ $registry_register{'_DBA'} }, $adap );
  }

}



=head2 get_DBAdaptor

  Arg [1]    : name of the species to get the adaptor for in the registry.
  Arg [2]    : name of the group to get the adaptor for in the registry.
  Example    : $dba = Bio::EnsEMBL::Registry->get_DBAdaptor("Human", "core");
  Returntype : DBAdaptor
  Exceptions : none
  Status     : Stable

=cut

sub get_DBAdaptor {
  my ( $class, $species, $group ) = @_;

  $species = $class->get_alias($species);

  return $registry_register{_SPECIES}{$species}{ lc($group) }{'_DB'};
}

=head2 get_all_DBAdaptors

  Arg [SPECIES]: (optional) string 
                  species name to get adaptors for
  Arg [GROUP]  : (optional) string 
                  group name to get adaptors for
  Example      : 
                @dba =
                  @{ Bio::EnsEMBL::Registry->get_all_DBAdaptors() };

                @human_dbas =
                  @{ Bio::EnsEMBL::Registry->get_all_DBAdaptors(
                    -species => 'human'
                  ) };

  Returntype   : list of DBAdaptors
  Exceptions   : none
  Status       : Stable

=cut

sub get_all_DBAdaptors {
  my ( $class, @args ) = @_;

  my ( $species, $group ) = rearrange( [qw(SPECIES GROUP)], @args );

  if ( defined($species) ) { $species = $class->get_alias($species) }

  my @ret;
  foreach my $dba ( @{ $registry_register{'_DBA'} } ) {
    if ( ( !defined($species) || lc($species) eq lc( $dba->species() ) )
      && ( !defined($group) || lc($group) eq lc( $dba->group() ) ) )
    {
      push( @ret, $dba );
    }
  }

  return \@ret;
}

=head2 get_all_DBAdaptors_by_connection

  Arg [1]    : DBConnection used to find DBAdaptors
  Returntype : reference to list of DBAdaptors
  Exceptions : none
  Example    : @dba = @{ Bio::EnsEMBL::Registry
                  ->get_all_DBAdaptors_by_connection($dbc) };
  Status     : Stable

=cut

sub get_all_DBAdaptors_by_connection {
  my ( $self, $dbc_orig ) = @_;

  my @return;

  foreach my $dba ( @{ $registry_register{'_DBA'} } ) {
    my $dbc = $dba->dbc();

    if (    defined($dbc)
         && $dbc->can('equals')
         && $dbc->equals($dbc_orig) )
    {
      push( @return, $dba );
    }
  }

  return \@return;
}

=head2 get_all_DBAdaptors_by_dbname

  Arg [1]    : string, name of database
  Returntype : reference to list of DBAdaptors
  Exceptions : none
  Example    : @dba = @{ Bio::EnsEMBL::Registry
                  ->get_all_DBAdaptors_by_dbname($dbname) };
  Status     : Stable

=cut

sub get_all_DBAdaptors_by_dbname {
  my ( $self, $dbname ) = @_;

  my @return;

  foreach my $dba ( @{ $registry_register{'_DBA'} } ) {
    my $dbc = $dba->dbc();

    if ( defined($dbc) && $dbc->dbname() eq $dbname ) {
      push( @return, $dba );
    }
  }

  return \@return;
}

=head2 remove_DBAdaptor

  Arg [1]    : name of the species to get the adaptor for in the registry.
  Arg [2]    : name of the group to get the adaptor for in the registry.
  Example    : $dba = Bio::EnsEMBL::Registry->remove_DBAdaptor("Human", "core");
  Returntype : none
  Exceptions : none
  Status     : At risk

=cut

sub remove_DBAdaptor {
  my ( $class, $species, $group ) = @_;

  $species = $class->get_alias($species);

  delete $registry_register{_SPECIES}{$species}{$group};
  # This will remove the DBAdaptor and all the other adaptors

  # Now remove if from the _DBA array
  my $index;

  foreach my $i ( 0 .. $#{ $registry_register{'_DBA'} } ) {
    my $dba = $registry_register{'_DBA'}->[$i];

    if ( ( $dba->species eq $species )
      && $dba->group eq $group )
    {
      $index = $i;
      last;
    }
  }

  # Now remove from _DBA cache
  if ( defined($index) ) {
    splice( @{ $registry_register{'_DBA'} }, $index, 1 );
  }

} ## end sub remove_DBAdaptor



=head2 reset_DBAdaptor

  Arg [1]:     string - species e.g. homo_sapiens
  Arg [2]:     string - DB group e.g. core
  Arg [3]:     string - new dbname
  Args [4-7]:  string - optional DB parameters, defaults to current db params if omitted
  Arg [8]:     hashref - Hash ref of additional parameters e.g. eFG dnadb params for auto selecting dnadb
  Usage :      $reg->reset_registry_db( 'homo_sapiens', 'core',
                  'homo_sapiens_core_37_35j' );
  Description: Resets a DB within the registry.
  Exceptions:  Throws if mandatory params not supplied
               Throws if species name is not already seen by the registry
               Throws if no current DB for species/group available
  Status :     At risk

=cut

sub reset_DBAdaptor {
  my (
    $self, $species, $group, $dbname, $host,
    $port, $user,    $pass,  $params
  ) = @_;

  # Check mandatory params
  if ( !( defined $species && defined $group && defined $dbname ) ) {
    throw(
      'Must provide at least a species, group, and dbname parameter '
        . 'to redefine a DB in the registry' );
  }

  # Validate species here
  my $alias = $self->get_alias($species);
  throw("Could not find registry alias for species:\t$species")
    if ( !defined $alias );

  # Get all current defaults if not defined

  my $db = $self->get_DBAdaptor( $alias, $group );
  my $class;

  if ($db) {
    $class = ref($db);
    $host ||= $db->dbc->host;
    $port ||= $db->dbc->port;
    $user ||= $db->dbc->username;
    $pass ||= $db->dbc->password;
  } else {
    #Now we need to test mandatory params
    $class = $group2adaptor{ lc($group) };

    if ( !( $host && $user ) ) {
      throw("No comparable $alias $group DB present in Registry. "
          . "You must pass at least a dbhost and dbuser" );
    }
  }

  $self->remove_DBAdaptor( $alias, $group );

  #my @adaptors = @{$self->get_all_adaptors};
  #This is causing a loop as it was constantly trying to reset the db
  #and never getting there.
  #I think this was left over from testing

  # ConfigRegistry should automatically add this to the Registry

  $db = $class->new(
    -user    => $user,
    -host    => $host,
    -port    => $port,
    -pass    => $pass,
    -dbname  => $dbname,
    -species => $alias,
    -group   => $group,
    %{$params} );

  return $db;
} ## end sub reset_DBAdaptor


#
# DNA Adaptors
#

=head2 add_DNAAdaptor

  Arg [1]    : name of the species to add the adaptor to in the registry.
  Arg [2]    : name of the group to add the adaptor to in the registry.
  Arg [3]    : name of the species to get the dna from
  Arg [4]    : name of the group to get the dna from
  Example    : Bio::EnsEMBL::Registry->add_DNAAdaptor("Human", "estgene", "Human", "core");
  Returntype : none
  Exceptions : none
  Status     : Stable

=cut

sub add_DNAAdaptor {
  my ( $class, $species, $group, $dnadb_species, $dnadb_group ) = @_;

  $species       = $class->get_alias($species);
  $dnadb_species = $class->get_alias($dnadb_species);
  if ( $dnadb_group->isa('Bio::EnsEMBL::DBSQL::DBAdaptor') ) {
    deprecated("");
  } else {
    $registry_register{_SPECIES}{$species}{ lc($group) }{'_DNA'} =
      $dnadb_group;
    $registry_register{_SPECIES}{$species}{ lc($group) }{'_DNA2'} =
      $dnadb_species;
  }
}

=head2 get_DNAAdaptor

  Arg [1]    : name of the species to get the adaptor for in the registry.
  Arg [2]    : name of the group to get the adaptor for in the registry.
  Example    : $dnaAdap = Bio::EnsEMBL::Registry->get_DNAAdaptor("Human", "core");
  Returntype : adaptor
  Exceptions : none
  Status     : Stable

=cut

sub get_DNAAdaptor {
  my ( $class, $species, $group ) = @_;

  $species = $class->get_alias($species);
  my $new_group =
    $registry_register{_SPECIES}{$species}{ lc($group) }{'_DNA'};
  my $new_species =
    $registry_register{_SPECIES}{$species}{ lc($group) }{'_DNA2'};

  if ( defined $new_group ) {
    return $class->get_DBAdaptor( $new_species, $new_group );
  }

  return undef;
}

#
# General Adaptors
#

=head2 add_adaptor

  Arg [1]    : name of the species to add the adaptor to in the registry.
  Arg [2]    : name of the group to add the adaptor to in the registry.
  Arg [3]    : name of the type to add the adaptor to in the registry.
  Arg [4]    : The DBAaptor to be added to the registry.
  Arg [5]    : (optional) if set okay to overwrite.
  Example    : Bio::EnsEMBL::Registry->add_adaptor("Human", "core", "Gene", $adap);
  Returntype : none
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub add_adaptor {
  my ( $class, $species, $group, $type, $adap, $reset ) = @_;

  $species = $class->get_alias($species);

  # Since the adaptors are not stored initially, only their class paths
  # when the adaptors are obtained, we need to store these instead.  It
  # is not necessarily an error if the registry is overwritten without
  # the reset set but it is an indication that we are overwriting a
  # database which should be a warning for now

  if ( defined($reset) )
  {    # JUST REST THE HASH VALUE NO MORE PROCESSING NEEDED
    $registry_register{_SPECIES}{$species}{ lc($group) }{ lc($type) } =
      $adap;
    return;
  }

  if (
    defined(
      $registry_register{_SPECIES}{$species}{ lc($group) }{ lc($type) }
    ) )
  {
  # print STDERR (
  #      "Overwriting Adaptor in Registry for $species $group $type\n");
    $registry_register{_SPECIES}{$species}{ lc($group) }{ lc($type) } =
      $adap;
    return;
  }
  $registry_register{_SPECIES}{$species}{ lc($group) }{ lc($type) } =
    $adap;

  if ( !defined( $registry_register{_SPECIES}{$species}{'list'} ) ) {
    $registry_register{_SPECIES}{$species}{'list'} = [$type];
  } else {
    push( @{ $registry_register{_SPECIES}{$species}{'list'} }, $type );
  }

  if ( !defined( $registry_register{_TYPE}{ lc($type) }{$species} ) ) {
    $registry_register{_TYPE}{ lc($type) }{$species} = [$type];
  } else {
    push( @{ $registry_register{_TYPE}{ lc($type) }{$species} },
      $adap );
  }

} ## end sub add_adaptor


=head2 get_adaptor

  Arg [1]    : name of the species to add the adaptor to in the registry.
  Arg [2]    : name of the group to add the adaptor to in the registry.
  Arg [3]    : name of the type to add the adaptor to in the registry.
  Example    : $adap = Bio::EnsEMBL::Registry->get_adaptor("Human", "core", "Gene");
  Returntype : adaptor
  Exceptions : none
  Status     : Stable

=cut

sub get_adaptor {
  my ( $class, $species, $group, $type ) = @_;

  $species = $class->get_alias($species);

  my %dnadb_adaptors = (
    'sequence'                 => 1,
    'assemblymapper'           => 1,
    'karyotypeband'            => 1,
    'repeatfeature'            => 1,
    'coordsystem'              => 1,
    'assemblyexceptionfeature' => 1
  );

  ## warn "$species, $group, $type";

  $type = lc($type);

  my $dnadb_group =
    $registry_register{_SPECIES}{$species}{ lc($group) }{'_DNA'};

  if ( defined($dnadb_group)
    && defined( $dnadb_adaptors{ lc($type) } ) )
  {
    $species =
      $registry_register{_SPECIES}{$species}{ lc($group) }{'_DNA2'};
    $group = $dnadb_group;
  }

  my $ret =
    $registry_register{_SPECIES}{$species}{ lc($group) }{ lc($type) };

  if ( !defined($ret) ) { return undef }
  if ( ref($ret) )      { return $ret }

  # Not instantiated yet

  my $dba = $registry_register{_SPECIES}{$species}{ lc($group) }{'_DB'};
  my $module = $ret;

  eval "require $module";
  if ($@) {
    warning("'$module' cannot be found.\nException $@\n");
    return undef;
  }

  if (
    !defined(
      $registry_register{_SPECIES}{$species}{ lc($group) }{'CHECKED'} )
    )
  {
    $registry_register{_SPECIES}{$species}{ lc($group) }{'CHECKED'} = 1;
    $class->version_check($dba);
  }

  my $adap = "$module"->new($dba);
  Bio::EnsEMBL::Registry->add_adaptor( $species, $group, $type, $adap,
    'reset' );
  $ret = $adap;

  return $ret;
} ## end sub get_adaptor

=head2 get_all_adaptors

  Arg [SPECIES] : (optional) string 
                  species name to get adaptors for
  Arg [GROUP] : (optional) string 
                  group name to get adaptors for
  Arg [TYPE] : (optional) string 
                  type to get adaptors for
  Example    : @adaps = @{Bio::EnsEMBL::Registry->get_all_adaptors()};
  Returntype : ref to list of adaptors
  Exceptions : none
  Status     : Stable

=cut

sub get_all_adaptors{
  my ($class,@args)= @_;
  my ($species, $group, $type);
  my @ret=();
  my (%species_hash, %group_hash, %type_hash);


  if(@args == 1){ # Old species only one parameter
    warn("-SPECIES argument should now be used to get species adaptors");
    $species = $args[0];
  }
  else{
    # new style -SPECIES, -GROUP, -TYPE
    ($species, $group, $type) =
      rearrange([qw(SPECIES GROUP TYPE)], @args);
  }

  if(defined($species)){
    $species_hash{$species} = 1;
  }
  else{
    # get list of species
    foreach my $dba (@{$registry_register{'_DBA'}}){
      $species_hash{lc($dba->species())} = 1;
    }
  }
  if(defined($group)){
    $group_hash{$group} = 1;
  }
  else{
    foreach my $dba (@{$registry_register{'_DBA'}}){
      $group_hash{lc($dba->group())} = 1;
    }
  }

  if ( defined($type) ) {
    $type_hash{$type} = 1;
  } else {
    foreach my $dba ( @{ $registry_register{'_DBA'} } ) {
      foreach my $ty (
        @{ $registry_register{_SPECIES}{ lc( $dba->species ) }{'list'} }
        )
      {
        $type_hash{ lc($ty) } = 1;
      }
    }
  }

  ### NOW NEED TO INSTANTIATE BY CALLING get_adaptor
  foreach my $sp ( keys %species_hash ) {
    foreach my $gr ( keys %group_hash ) {
      foreach my $ty ( keys %type_hash ) {
        my $temp = $class->get_adaptor( $sp, $gr, $ty );
        if ( defined($temp) ) {
          push @ret, $temp;
        }
      }
    }
  }

  return (\@ret);
}


=head2 add_alias

  Arg [1]    : name of the species to add alias for
  Arg [2]    : name of the alias
  Example    : Bio::EnsEMBL::Registry->add_alias("Homo Sapiens","Human");
  Description: add alternative name for the species.
  Returntype : none
  Exceptions : none
  Status     : Stable

=cut

sub add_alias{
  my ($class, $species,$key) = @_;

  $registry_register{'_ALIAS'}{lc($key)} = lc($species);
}

=head2 get_alias

  Arg [1]    : name of the possible alias to get species for
  Example    : Bio::EnsEMBL::Registry->get_alias("Human");
  Description: get proper species name.
  Returntype : species name
  Exceptions : none
  Status     : Stable

=cut

sub get_alias{
  my ($class, $key) = @_;

  if(!defined($registry_register{'_ALIAS'}{lc($key)})){
    return $key;
  }
  return $registry_register{'_ALIAS'}{lc($key)};
}

=head2 alias_exists

  Arg [1]    : name of the possible alias to get species for
  Example    : Bio::EnsEMBL::Registry->alias_exists("Human");
  Description: does the species name exist.
  Returntype : 1 if exists else 0
  Exceptions : none
  Status     : Stable

=cut

sub alias_exists {
  my ( $class, $key ) = @_;

  return defined( $registry_register{'_ALIAS'}{ lc($key) } );
}

=head2 set_disconnect_when_inactive

  Example    : Bio::EnsEMBL::Registry->set_disconnect_when_inactive();
  Description: Set the flag to make sure that the database connection is dropped if
               not being used on each database.
  Returntype : none
  Exceptions : none
  Status     : Stable

=cut

sub set_disconnect_when_inactive{
  foreach my $dba ( @{get_all_DBAdaptors()}){
    my $dbc = $dba->dbc;
    # Disconnect if connected
    $dbc->disconnect_if_idle() if $dbc->connected();
    $dbc->disconnect_when_inactive(1);
  }
}


=head2 disconnect_all

  Example    : Bio::EnsEMBL::Registry->disconnect_all();
  Description: disconnect from all the databases.
  Returntype : none
  Exceptions : none
  Status     : Stable

=cut

sub disconnect_all {
  foreach my $dba ( @{get_all_DBAdaptors()||[]} ){
    my $dbc = $dba->dbc;
    next unless $dbc;
    # Disconnect if connected
    $dbc->disconnect_if_idle() if $dbc->connected();
  }
}

=head2 change_access

  Will change the username and password for a set of databases.
  if host,user or database names are missing then these are not checked.
  So for example if you do not specify a database then ALL databases on
  the specified  host and port will be changed.

  Arg [1]    : name of the host to change access on
  Arg [2]    : port number to change access on
  Arg [3]    : name of the user to change access on
  Arg [4]    : name of the database to change access on
  Arg [5]    : name of the new user
  Arg [6]    : new password

  Example    : Bio::EnsEMBL::Registry->get_alias("Human");
  Description: change username and password on one or more databases
  Returntype : none
  Exceptions : none
  Status     : Stable

=cut

sub change_access{
my $self = shift;
    my ($host,$port,$user,$dbname,$new_user,$new_pass) = @_;
    foreach my $dba ( @{$registry_register{'_DBA'}}){
	my $dbc = $dba->dbc;
	if((!defined($host) or $host eq $dbc->host) and
	   (!defined($port) or $port eq $dbc->port) and
	   (!defined($user) or $user eq $dbc->username) and
	   (!defined($dbname) or $dbname eq $dbc->dbname)){
	    if($dbc->connected()){
		$dbc->db_handle->disconnect();
		$dbc->connected(undef);
	    }
	    # over write the username and password
	    $dbc->username($new_user);
	    $dbc->password($new_pass);
	}
    }
}



=head2 load_registry_from_url

  Arg [1] : string $url
  Arg [2] : (optional) integer
            If not 0, will print out all information.
  Arg [3] : (optional) integer
          This option will turn off caching for slice features,
          so, every time a set of features is retrieved, they
          will come from the database instead of the cache. This
          option is only recommended for advanced users, specially
          if you need to store and retrieve features. It might
          reduce performance when querying the database if not used
          properly. If in doubt, do not use it or ask in ensembl-dev

  Example : load_registry_from_url(
            'mysql://anonymous@ensembldb.ensembl.org:3306');

  Description: Will load the correct versions of the ensembl
               databases for the software release it can find on
               a database instance into the registry. Also adds
               a set of standard aliases. The url format is:
               mysql://[[username][:password]@]hostname[:port].  You
               can also request a specific version for the databases
               by adding a slash and the version number but your
               script may crash as the API version won't match the
               DB version.

  Exceptions : None.
  Status     : Stable
 
=cut

sub load_registry_from_url {
  my ( $self, $url, $verbose, $no_cache ) = @_;

  if ( $url =~ /mysql\:\/\/([^\@]+\@)?([^\:\/]+)(\:\d+)?(\/\d+)?/ ) {
    my $user_pass = $1;
    my $host      = $2;
    my $port      = $3;
    my $version   = $4;

    $user_pass =~ s/\@$//;
    my ( $user, $pass ) = $user_pass =~ m/([^\:]+)(\:.+)?/;
    $pass    =~ s/^\:// if ($pass);
    $port    =~ s/^\:// if ($port);
    $version =~ s/^\/// if ($version);

    $self->load_registry_from_db(
      -host       => $host,
      -user       => $user,
      -pass       => $pass,
      -port       => $port,
      -db_version => $version,
      -verbose    => $verbose,
      -no_cache   => $no_cache
    );
  } else {
    throw("Only MySQL URLs are accepted at the moment");
  }
} ## end sub load_registry_from_url


=head2 load_registry_from_db

  Arg [HOST] : string
                The domain name of the database host to connect to.

  Arg [USER] : string
                The name of the database user to connect with.

  Arg [PASS] : (optional) string
                The password to be used to connect to the database.

  Arg [PORT] : (optional) integer
                The port to use when connecting to the database.

  Arg [VERBOSE]: (optional) boolean
                Whether to print database messages.

  Arg [DB_VERSION]: (optional) integer
                By default, only databases corresponding to this API
                version are loaded. This allows the script to use
                databases from another version although it might not
                work properly.  This option should only be used for
                production or testing purposes and if you really
                know what you are doing.

  Arg [WAIT_TIMEOUT]: (optional) integer
                Time in seconds for the wait timeout to happen.
                Time after which the connection is deleted if not
                used.  By default this is 28800 (8 hours), so set
                this to greater than this if your connection are
                getting deleted.  Only set this if you are having
                problems and know what you are doing.

   Arg [-NO_CACHE]: (optional) int 1
                This option will turn off caching for slice
                features, so, every time a set of features is
                retrieved, they will come from the database instead
                of the cache.  This option is only recommended for
                advanced users, specially if you need to store and
                retrieve features.  It might reduce performance when
                querying the database if not used properly.  If in
                doubt, do not use it or ask in ensembl-dev.

  Example :

    $registry->load_registry_from_db(
      -host    => 'ensembldb.ensembl.org',
      -user    => 'anonymous',
      -verbose => '1'
    );

  Description: Will load the correct versions of the ensembl
               databases for the software release it can find on a
               database instance into the registry.  Also adds a set
               of standard aliases.

  Exceptions : None.
  Status     : Stable

=cut

sub load_registry_from_db {
  my ( $self, @args ) = @_;

  my ( $host, $port, $user, $pass, $verbose, $db_version, $wait_timeout,
    $no_cache )
    = rearrange( [
      'HOST',    'PORT',       'USER',         'PASS',
      'VERBOSE', 'DB_VERSION', 'WAIT_TIMEOUT', 'NO_CACHE'
    ],
    @args
    );

  my $go_version       = 0;
  my $ontology_version = 0;

  $user ||= "ensro";
  if ( !defined($port) ) {
    $port = 3306;
    if ( $host eq "ensembldb.ensembl.org" ) {
      if ( !defined($db_version) or $db_version >= 48 ) {
        $port = 5306;
      }
    }
  }

  $wait_timeout ||= 0;

  my $dbh =
    DBI->connect( "DBI:mysql:host=$host;port=$port", $user, $pass );

  my $res = $dbh->selectall_arrayref('SHOW DATABASES');
  my @dbnames = map { $_->[0] } @$res;

  my %temp;
  my $software_version = $self->software_version();

  if ( defined($db_version) ) {
    $software_version = $db_version;
  }

  if ($verbose) {
    printf( "Will only load v%d databases\n", $software_version );
  }

  for my $db (@dbnames) {
    if ( $db =~ /^(\w+)_(collection_core_(?:\d+_)?(\d+)_(\w+))/ )
    {    # NEEDS TO BE FIRST
      if ( $3 eq $software_version ) {
        $temp{$1} = $2;
      }
    } elsif ( $db =~ /^(.+)_(userdata)$/ ) {
      $temp{$1} = $2;
    } elsif ( $db =~ /^(ensembl_compara(?:_\w+)*?)_(\d+)$/ ) {
      if ( $2 eq $software_version ) {
        $temp{$1} = $2;
      }
    } elsif ( $db =~ /^(ensembl_ancestral(?:_\w+?)*?)_(\d+)$/ ) {
      if ( $2 eq $software_version ) {
        $temp{$1} = $2;
      }
    } elsif ( $db =~ /^ensembl_go_(\d+)/ ) {
      if ( $1 eq $software_version ) {
        $go_version = $1;
      }
    } elsif ( $db =~ /^(ensembl_ontology)_(\d+)/ ) {
      if ( $2 eq $software_version ) {
        $ontology_version = $2;
      }
    } elsif (
      $db =~ /^([a-z]+_[a-z]+_[a-z]+(?:_\d+)?)_(\d+)_(\w+)/ )
    {
      if ( $2 eq $software_version ) {
        $temp{$1} = $2 . "_" . $3;
      }
    } else {
      # warn( sprintf( "Skipping database '%s'\n", $db ) );
    }
  } ## end for my $db (@dbnames)

  @dbnames = ();

  foreach my $key ( keys %temp ) {
    push @dbnames, $key . "_" . $temp{$key};
  }

  # Register Core databases

  my @core_dbs = grep { /^[a-z]+_[a-z]+_core_(?:\d+_)?\d+_/ } @dbnames;

  foreach my $coredb (@core_dbs) {
    next if ($coredb =~ /collection/);  # Skip multi-species databases

    my ( $species, $num ) =
      ( $coredb =~ /(^[a-z]+_[a-z]+)_core_(?:\d+_)?(\d+)/ );

    my $dba =
      Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                         -group        => "core",
                                         -species      => $species,
                                         -host         => $host,
                                         -user         => $user,
                                         -pass         => $pass,
                                         -port         => $port,
                                         -dbname       => $coredb,
                                         -wait_timeout => $wait_timeout,
                                         -no_cache     => $no_cache );

    if ($verbose) {
      printf( "Species '%s' loaded from database '%s'\n",
              $species, $coredb );
    }
  }

  # Register multi-species databases
  my @multi_dbs = grep { /^\w+_collection_core_\w+$/ } @dbnames;

  foreach my $multidb (@multi_dbs) {
    my $sth =
      $dbh->prepare(
                 sprintf( 'SELECT species_id, meta_value FROM %s.meta ',
                          $dbh->quote_identifier($multidb) )
                   . "WHERE meta_key = 'species.db_name'" );

    $sth->execute();

    my ( $species_id, $species );
    $sth->bind_columns( \( $species_id, $species ) );

    while ( $sth->fetch() ) {
      my $dba =
        Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                         -group      => "core",
                                         -species    => $species,
                                         -species_id => $species_id,
                                         -multispecies_db => 1,
                                         -host            => $host,
                                         -user            => $user,
                                         -pass            => $pass,
                                         -port            => $port,
                                         -dbname          => $multidb,
                                         -wait_timeout => $wait_timeout,
                                         -no_cache     => $no_cache );

      if ($verbose) {
        printf( "Species '%s' (id:%d) loaded from database '%s'\n",
                $species, $species_id, $multidb );
      }
    }
  } ## end foreach my $multidb (@multi_dbs)

  # register cdna databases

  my @cdna_dbs = grep { /^[a-z]+_[a-z]+_cdna_(?:\d+_)?\d+_/ } @dbnames;

  for my $cdnadb (@cdna_dbs) {
    my ( $species, $num ) =
      ( $cdnadb =~ /(^[a-z]+_[a-z]+)_cdna_(?:\d+_)?(\d+)_/ );
    my $dba =
      Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                         -group        => "cdna",
                                         -species      => $species,
                                         -host         => $host,
                                         -user         => $user,
                                         -pass         => $pass,
                                         -port         => $port,
                                         -dbname       => $cdnadb,
                                         -wait_timeout => $wait_timeout,
                                         -no_cache     => $no_cache );

    if ($verbose) {
      printf( "%s loaded\n", $cdnadb );
    }
  }

  my @vega_dbs = grep { /^[a-z]+_[a-z]+_vega_\d+_/ } @dbnames;

  for my $vegadb (@vega_dbs) {
    my ( $species, $num ) =
      ( $vegadb =~ /(^[a-z]+_[a-z]+)_vega_(\d+)/ );
    my $dba =
      Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                         -group        => "vega",
                                         -species      => $species,
                                         -host         => $host,
                                         -user         => $user,
                                         -pass         => $pass,
                                         -port         => $port,
                                         -wait_timeout => $wait_timeout,
                                         -dbname       => $vegadb,
                                         -no_cache     => $no_cache );

    if ($verbose) {
      printf( "%s loaded\n", $vegadb );
    }
  }

  # Otherfeatures

  my @other_dbs = grep { /^[a-z]+_[a-z]+_otherfeatures_(?:\d+_)?\d+_/ } @dbnames;

  for my $other_db (@other_dbs) {
    my ( $species, $num ) =
      ( $other_db =~ /(^[a-z]+_[a-z]+)_otherfeatures_(?:\d+_)?(\d+)_/ );
    my $dba =
      Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                         -group   => "otherfeatures",
                                         -species => $species,
                                         -host    => $host,
                                         -user    => $user,
                                         -pass    => $pass,
                                         -port    => $port,
                                         -wait_timeout => $wait_timeout,
                                         -dbname       => $other_db,
                                         -no_cache     => $no_cache );

    if ($verbose) {
      printf( "%s loaded\n", $other_db );
    }
  }

  # User upload DBs

  my @userupload_dbs = grep { /_userdata$/ } @dbnames;
  for my $userupload_db (@userupload_dbs) {
    my ($species) = ( $userupload_db =~ /(^.+)_userdata$/ );
    my $dba =
      Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                         -group        => "userupload",
                                         -species      => $species,
                                         -host         => $host,
                                         -user         => $user,
                                         -pass         => $pass,
                                         -port         => $port,
                                         -wait_timeout => $wait_timeout,
                                         -dbname   => $userupload_db,
                                         -no_cache => $no_cache );

    if ($verbose) {
      printf( "%s loaded\n", $userupload_db );
    }
  }

  # Variation

  eval "require Bio::EnsEMBL::Variation::DBSQL::DBAdaptor";
  if ($@) {
    # Ignore variations as code required not there for this
    if ($verbose) {
      print(
           "Bio::EnsEMBL::Variation::DBSQL::DBAdaptor module not found "
             . "so variation databases will be ignored if found\n" );
    }
  } else {
    my @variation_dbs =
      grep { /^[a-z]+_[a-z]+_variation_(?:\d+_)?\d+_/ } @dbnames;

    for my $variation_db (@variation_dbs) {
      my ( $species, $num ) =
        ( $variation_db =~ /(^[a-z]+_[a-z]+)_variation_(?:\d+_)?(\d+)_/ );
      my $dba =
        Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(
                                         -group        => "variation",
                                         -species      => $species,
                                         -host         => $host,
                                         -user         => $user,
                                         -pass         => $pass,
                                         -port         => $port,
                                         -wait_timeout => $wait_timeout,
                                         -dbname       => $variation_db,
                                         -no_cache     => $no_cache );

      if ($verbose) {
        printf( "%s loaded\n", $variation_db );
      }
    }
  }

  eval "require Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor";
  if ($@) {
    if ($verbose) {
      # Ignore funcgen DBs as code required not there for this
      print( "Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor module not found "
         . "so functional genomics databases will be ignored if found\n"
      );
    }
  } else {
    my @funcgen_dbs = grep { /^[a-z]+_[a-z]+_funcgen_(?:\d+_)?\d+_/ } @dbnames;

    for my $funcgen_db (@funcgen_dbs) {
      my ( $species, $num ) =
        ( $funcgen_db =~ /(^[a-z]+_[a-z]+)_funcgen_(?:\d+_)?(\d+)_/ );
      my $dba =
        Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
                                         -group        => "funcgen",
                                         -species      => $species,
                                         -host         => $host,
                                         -user         => $user,
                                         -pass         => $pass,
                                         -port         => $port,
                                         -wait_timeout => $wait_timeout,
                                         -dbname       => $funcgen_db,
                                         -no_cache     => $no_cache );

      if ($verbose) {
        printf( "%s loaded\n", $funcgen_db );
      }
    }
  }

  # Compara

  my @compara_dbs = grep { /^ensembl_compara/ } @dbnames;

  if (@compara_dbs) {
    eval "require Bio::EnsEMBL::Compara::DBSQL::DBAdaptor";
    if ($@) {
      # Ignore Compara as code required not there for this
      if ($verbose) {
        printf(
          "Bio::EnsEMBL::Compara::DBSQL::DBAdaptor "
            . "not found so the following compara "
            . "databases will be ignored: %s\n",
          join( ', ', @compara_dbs ) );
      }
    } else {
      foreach my $compara_db (@compara_dbs) {
        # Looking for EnsEMBL Genomes Comparas.
        # ensembl_compara_bacteria_2_53 is registered as
        # 'bacteria', ensembl_compara_pan_homology_2_53 is
        # registered as 'pan_homology', ensembl_compara_53 is
        # registered as 'multi', and the alias 'compara' still
        # operates.

        my ($species) =
          $compara_db =~ /^ensembl_compara_(\w+)(?:_\d+){2}$/xm;

        $species ||= 'multi';

        my $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
          -group        => 'compara',
          -species      => $species,
          -host         => $host,
          -user         => $user,
          -pass         => $pass,
          -port         => $port,
          -wait_timeout => $wait_timeout,
          -dbname       => $compara_db,
          -no_cache     => $no_cache
        );

        if ($verbose) {
          printf( "%s loaded\n", $compara_db );
        }
      } ## end foreach my $compara_db (@compara_dbs)
    } ## end else [ if ($@)
  } elsif ($verbose) {
    print("No Compara databases found\n");
  }

  # Ancestral sequences

  my @ancestral_dbs =
    sort grep { /^ensembl_ancestral/ } @dbnames;

  if (@ancestral_dbs) {
    my $ancestral_db = shift @ancestral_dbs;

    my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
      -group        => 'core',
      -species      => 'Ancestral sequences',
      -host         => $host,
      -user         => $user,
      -pass         => $pass,
      -port         => $port,
      -wait_timeout => $wait_timeout,
      -dbname       => $ancestral_db,
      -no_cache     => $no_cache
    );

    if ($verbose) {
      printf( "%s loaded\n", $ancestral_db );

      if (@ancestral_dbs) {
        # If we still had some more then report the problem.
        printf(
          "Multiple ancestral databases found.\n"
            . "Ignoring the following: %s\n",
          join( ', ', @ancestral_dbs ) );
      }
    }
  } elsif ($verbose) {
    print("No ancestral database found\n");
  }

  # GO

  if ($go_version) {
    eval "require Bio::EnsEMBL::ExternalData::GO::GOAdaptor";
    if ($@) {
      #ignore go as code required not there for this
      #      print $@;
      if ($verbose) {
        print "GO software not installed "
          . "so GO database ensembl_go_$go_version will be ignored\n";
      }
    } else {
      my $go_db = "ensembl_go_" . $go_version;
      my $dba =
        Bio::EnsEMBL::ExternalData::GO::GOAdaptor->new(
                                                  -group    => "go",
                                                  -species  => "multi",
                                                  -host     => $host,
                                                  -user     => $user,
                                                  -pass     => $pass,
                                                  -port     => $port,
                                                  -dbname   => $go_db,
                                                  -no_cache => $no_cache
        );

      if ($verbose) {
        printf( "%s loaded\n", $go_db );
      }
    }
  } elsif ($verbose) {
    print("No GO database found\n");
  }

  # Ontology

  if ( $ontology_version != 0 ) {
    require Bio::EnsEMBL::DBSQL::OntologyDBAdaptor;

    my $ontology_db =
      sprintf( "ensembl_ontology_%d", $ontology_version );

    my $dba = Bio::EnsEMBL::DBSQL::OntologyDBAdaptor->new(
      '-species' => 'multi',
      '-group'   => 'ontology',
      '-host'    => $host,
      '-port'    => $port,
      '-user'    => $user,
      '-pass'    => $pass,
      '-dbname'  => $ontology_db,
    );

    if ($verbose) {
      printf( "%s loaded\n", $ontology_db );
    }
  } elsif ($verbose) {
    print("No ontology database found\n");
  }

  Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
    -species => 'multi',
    -alias   => ['compara'] );

  Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
    -species => 'multi',
    -alias   => ['go'] );

  Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
    -species => 'multi',
    -alias   => ['ontology'] );

  # Register aliases as found in adaptor meta tables.
  $self->find_and_add_aliases( '-handle' => $dbh );
  $dbh->disconnect();

} ## end sub load_registry_from_db

=head2 find_and_add_aliases

  Arg [ADAPTOR] : (optional) Bio::EnsEMBL::DBSQL::DBAdaptor
                  The adaptor to use to retrieve aliases from.

  Arg [GROUP]   : (optional) string
                  The group you want to find aliases for. If not
                  given assumes all types.

  Arg [HANDLE]  : (optional) DBI database handle
                  A connected database handle to use instead of the
                  database handles stored in the DBAdaptors.  Bypasses
                  the use of MetaContainer.

  Example       : Bio::EnsEMBL::Registry->find_and_add_aliases(
                    -ADAPTOR => $dba,
                    -GROUP   => 'core'
                  );

  Description   : Looks in the meta container for each database for
                  an entry called "species.alias".  If any are found
                  then the species adaptor is registered to that
                  set of aliases.  This can work across any adaptor
                  which has a MetaContainer.  If no MetaContainer
                  can be returned from a given adaptor then no alias
                  searching is performed.

  Return type   : none
  Exceptions    : Throws if an alias is found in more than one species.
  Status        : Stable

=cut

sub find_and_add_aliases {
  my $class = shift @_;

  my ( $adaptor, $group, $dbh ) =
    rearrange( [ 'ADAPTOR', 'GROUP', 'HANDLE' ], @_ );

  my @dbas;
  if ( defined($adaptor) ) {
    @dbas = ($adaptor);
  } else {
    @dbas = @{ $class->get_all_DBAdaptors( '-GROUP' => $group ) };
  }

  foreach my $dba (@dbas) {
    my @aliases;
    my $species = $dba->species();

    if ( defined($dbh) ) {
      my $dbname = $dba->dbc()->dbname();
      my $sth    = $dbh->prepare(
        sprintf(
          "SELECT meta_value FROM %s.meta "
            . "WHERE meta_key = 'species.alias' "
            . "AND species_id = ?",
          $dbh->quote_identifier($dbname) ) );

      # Execute, and don't care about errors (there will be errors for
      # databases without a 'meta' table.
      $sth->{'PrintError'} = 0;
      $sth->{'RaiseError'} = 0;
      if ( !$sth->execute( $dba->species_id() ) ) { next }
      $sth->{'PrintError'} = $dbh->{'PrintError'};
      $sth->{'RaiseError'} = $dbh->{'RaiseError'};

      my $alias;
      $sth->bind_columns( \$alias );
      while ( $sth->fetch() ) {
        push( @aliases, $alias );
      }
    } else {
      my $meta_container = eval { $dba->get_MetaContainer() };

      if ( defined($meta_container) ) {
        foreach my $key ( qw(
          species.alias
          species.taxonomy_id
          assembly.name
          species.common_name
          ) )
        {
          push( @aliases,
            @{ $meta_container->list_value_by_key($key) } );
        }
      }

      # Need to disconnect so we do not spam the MySQL servers trying to
      # get aliases.  Can only call disonnect if dbc was defined.
      if ( defined( $dba->dbc() ) ) {
        $dba->dbc()->disconnect_if_idle();
      }
    }

    foreach my $alias (@aliases) {
      if ( !$class->alias_exists($alias) ) {
        $class->add_alias( $species, $alias );
      } elsif ( $species ne $class->get_alias($alias) ) {
        throw(
          sprintf(
            "Trying to add alias '%s' to species '%s', "
              . " but it is already registrered for species '%s'\n",
            $alias, $species, $class->get_alias($alias) ) );
      }
    }

  } ## end foreach my $dba (@dbas)

} ## end sub find_and_add_aliases


=head2 load_registry_from_multiple_dbs

  Arg [1]   : Array of hashes, each hash being a set of arguments to
              load_registry_from_db() (see above).

  Example   :

    $registry->load_registry_from_multiple_dbs( {
        '-host'    => 'ensembldb.ensembl.org',
        '-user'    => 'anonymous',
        '-verbose' => '1'
      },
      {
        '-host'     => 'server.example.com',
        '-user'     => 'anonymouse',
        '-password' => 'cheese',
        '-verbose'  => '1'
      } );

  Description:  Will call load_registry_from_db() (see above)
                multiple times and merge the resulting registries
                into one, effectively allowing a user to connect to
                databases on multiple database servers from within
                one program.

                If a database is found on more than one server, the
                first found instance of that database will be used.

=cut

sub load_registry_from_multiple_dbs {
  my ( $self, @args ) = @_;

  my %merged_register;

  foreach my $arg (@args) {
    local %registry_register;

    my $verbose;

    ($verbose) = rearrange( ['VERBOSE'], %{$arg} );

    $self->load_registry_from_db( %{$arg} );

    #
    # Merge the localized %registry_register into %merged_register.
    #

    # Merge the _SPECIES and _ALIAS sections of %registry_register.
    foreach my $section ( 'Species', 'Alias' ) {
      my $section_key = '_' . uc($section);

      while ( my ( $key, $value ) =
        each( %{ $registry_register{$section_key} } ) )
      {
        if ( !exists( $merged_register{$section_key}{$key} ) ) {
          $merged_register{$section_key}{$key} = $value;
        } elsif ($verbose) {
          printf( "%s '%s' found on multiple servers, "
              . "using first found\n",
            $section, $key );
        }
      }
    }
  } ## end foreach my $arg (@args)

  # Add the DBAs from the _SPECIES section into the _DBA section.
  foreach my $species_hash ( values( %{ $merged_register{_SPECIES} } ) )
  {
    foreach my $group_hash ( values( %{$species_hash} ) ) {
      if ( ref($group_hash) eq 'HASH' && exists( $group_hash->{_DB} ) )
      {
        push( @{ $merged_register{_DBA} }, $group_hash->{_DB} );
      }
    }
  }

  %registry_register = %merged_register;

} ## end sub load_registry_from_multiple_dbs

#
# Web specific routines
#

=head2 DEPRECATED load_registry_with_web_adaptors

  DEPRECATED: Use load_registry_from_db instead.

=cut

sub load_registry_with_web_adaptors{
  my $class = shift;

  deprecate('Use the load_registry_from_db instead'); 
  eval{ require SiteDefs };
  if ($@){ die "Can't use SiteDefs.pm - $@\n"; }
    SiteDefs->import(qw(:ALL));

  eval{ require SpeciesDefs };
  if ($@){ die "Can't use SpeciesDefs.pm - $@\n"; }
  my $conf = new SpeciesDefs();
  
  my %species_alias = %{$SiteDefs::ENSEMBL_SPECIES_ALIASES};

  foreach my $spec (keys %species_alias){
    Bio::EnsEMBL::Registry->add_alias($species_alias{$spec},$spec);
  }

}

=head2 set_default_track

  Sets a flag to say that that this species/group are a default track and do not
  need to be added as another web track.

  Arg [1]    : name of the species to get the adaptors for in the registry.
  Arg [2]    : name of the type to get the adaptors for in the registry.
  Example    : $merged = Bio::EnsEMBL::Registry->set_default_track("Human","core");
  Returntype : none
  Exceptions : none
  Status     : At Risk.

=cut

sub set_default_track {
  my ( $class, $species, $group ) = @_;

  $species = get_alias($species);
  $registry_register{'def_track'}{$species}{ lc($group) } = 1;
}

=head2 default_track

  Check flag to see if this is a default track

  Arg [1]    : name of the species to get the adaptors for in the registry.
  Arg [2]    : name of the type to get the adaptors for in the registry.
  Example    : $merged = Bio::EnsEMBL::Registry->set_default_track("Human","core");
  Returntype : int 
  Exceptions : none
  Status     : At Risk.

=cut

sub default_track {
  my ( $class, $species, $group ) = @_;

  $species = get_alias($species);
  if (
    defined( $registry_register{'def_track'}{$species}{ lc($group) } ) )
  {
    return 1;
  }

  return 0;
}


=head2 add_new_tracks

  Will add new gene tracks to the configuration of the WEB server if they are
  not of the type default and the configuration already has genes in the display.

  Arg [1]    : hash of the default configuration of the web page
  Returntype : none
  Exceptions : none
  Called by  : UserConfig.pm
  Status     : At Risk.
  
=cut

sub add_new_tracks{
  my($class, $conf, $pos) = @_;

  my $start = 0;
  my $reg = $class;
  my $species_reg = $reg->get_alias($conf->{'species'},"nothrow");
  my %pars;
#  print STDERR "Species $species_reg check for default tracks\n";
  if(defined($species_reg)){
    foreach my $dba (@{$reg->get_all_DBAdaptors()}){
      if(!$reg->default_track($dba->species,$dba->group)){
	$pars{'available'} = "species ".$reg->get_alias($dba->species());
	$pars{'db_alias'} = $dba->group();
#	print STDERR "Adding new track for ".$dba->species."\t".$dba->group."\n";
	$conf->add_new_track_generictranscript('',$dba->group(), "black",$pos,%pars);
	$pos++;
      }
    }
  }
  return $pos;

}

=head2 software_version
  
  get the software version.
  
  Args       : none
  ReturnType : int
  Status     : At Risk
  
=cut
  
sub software_version{
  my ($self) = @_;
  return $API_VERSION;
}
  
=head2 no_version_check
  
  getter/setter for whether to run the version checking
  
  Arg[0]     : (optional) int
  Returntype : int or undef if not set
  Exceptions : none
  Status     : At Risk.

=cut
  
sub no_version_check {
  my ( $self, $arg ) = @_;
  ( defined $arg )
    && ( $registry_register{'_no_version_check'} = $arg );

  return $registry_register{'_no_version_check'};
}

  
=head2 version_check
  
  run the database/API code version check for a DBAdaptor
  
  Arg[0]     : DBAdaptor to check
  Returntype : int 1 if okay, 0 if not the same 
  Exceptions : none
  Status     : At Risk.

=cut
  
  
sub version_check {
  my ( $self, $dba ) = @_;

  # Check the datbase and versions match
  # give warning if they do not.
  my $check = no_version_check();

  if ( (
      defined( $ENV{HOME} )
      and ( -e $ENV{HOME} . "/.ensemblapi_no_version_check" ) )
    or ( defined($check) and ( $check != 0 ) ) )
  {
    return 1;
  }

  my $mca =
    $self->get_adaptor( $dba->species(), $dba->group(),
    "MetaContainer" );

  my $database_version = 0;
  if ( defined($mca) ) {
    $database_version = $mca->get_schema_version();
  }

  if ( $database_version == 0 ) {
    # Try to work out the version
    if ( $dba->dbc()->dbname() =~ /^_test_db_/ ) {
      return 1;
    }
    if ( $dba->dbc()->dbname() =~ /(\d+)_\S+$/ ) {
      $database_version = $1;
    } elsif ( $dba->dbc()->dbname() =~ /ensembl_compara_(\d+)/ ) {
      $database_version = $1;
    } elsif ( $dba->dbc()->dbname() =~ /ensembl_go_(\d+)/ ) {
      $database_version = $1;
    } elsif ( $dba->dbc()->dbname() =~ /ensembl_help_(\d+)/ ) {
      $database_version = $1;
    } elsif ( $dba->dbc()->dbname() =~ /ensembl_ontology_(\d+)/ ) {
      $database_version = $1;
    } else {
      warn(
        sprintf(
          "No database version for database %s "
            . ". You must be using a post version 34 database "
            . "with version 34 or later code.\n"
            . "You need to update your database "
            . "or use the appropriate Ensembl software release "
            . "to ensure your script does not crash\n",
          $dba->dbc()->dbname() ) );
    }
  } ## end if ( $database_version...

  if ( $database_version != $API_VERSION ) {
    warn(
      sprintf(
        "For %s there is a difference in the software release (%s) "
          . "and the database release (%s). "
          . "You should update one of these to ensure that your script "
          . "does not crash.\n",
        $dba->dbc()->dbname(),
        $API_VERSION, $database_version
      ) );
    return 0;
  }

  return 1;    # Ok
} ## end sub version_check


=head2 get_species_and_object_type
  
  get the species, ensembl type (Gene, Transcript , Translation or Exon) 
  and database type (core, vega, otherfeatures) for a given stable_id
  
  Arg[1]     : stable_id to find species and ensembl type for.
  Arg[2]     : (optional) integer. force searching of other databases that do not have a set format
               for the stable_id by connecting to the databases and trying to retrieve each of the
               ensembl object fetching using the stabke_id.
  Example    : my ($species, $type, $db_type) = Bio::EnsEMBL::Registry->get_species_and_object_type(ENST00000326632);
  Returntype : array. consisting of the species name, ensembl type and database type (core,vega,otherfeatures). undef, undef, undef returned if not found
  Exceptions : none
  Status     : At Risk.

=cut

#hashes containing codes of species and different database types to return
  
our %ensembl_type = qw(T Transcript G Gene P Translation E Exon);
our %ensembl_species = qw(
			  ENS     Homo_sapiens
			  ENSRNO  Rattus_norvegicus
			  ENSMUS  Mus_musculus
			  ENSGAL  Gallus_gallus
			  ENSBTA  Bos_taurus
			  ENSDAR  Danio_rerio
			  ENSCAF  Canis_familiaris
			  ENSPTR  Pan_troglodytes
			  ENSCPO  Cavia_porcellus
			  ENSCIN  Ciona_intestinalis
			  ENSCSAV Ciona_savignyi
			  ENSDNO  Dasypus_novemcinctus
			  ENSETE  Echinops_telfairi
			  ENSEEU  Erinaceus_europaeus
			  ENSFCA  Felis_catus
			  ENSGAC  Gasterosteus_aculeatus
			  ENSLAF  Loxodonta_africana
			  ENSMMU  Macaca_mulatta
			  ENSMOD  Monodelphis_domestica
			  ENSMLU  Myotis_lucifugus
			  ENSOAN  Ornithorhynchus_anatinus
			  ENSOCU  Oryctolagus_cuniculus
			  ENSORL  Oryzias_latipes
			  ENSSAR  Otolemur_garnettii
			  ENSSTO  Spermophilus_tridecemlineatus
			  ENSTBE  Tupaia_belangeri
			  SINFRU  Takifugu_rubripes
			  ENSXET  Xenopus_tropicalis
			  ENSMEU  Macropus_eugenii
			  ENSSSC  Sus_scrofa
			  ENSCJA  Callithrix_jaccus
			  );
our %vega_species = qw(
		       OTTHUM  Homo_sapiens
		       OTTMUS  Mus_musculus
		       );
our %ensembl_db_type = qw(EST otherfeatures);
our %vectorbase_species = qw(
			     AAEL Aedes_aegypti
			     AGAP Anopheles_gambiae
			     );
our %vectorbase_type = qw(R Transcript P Translation);

sub get_species_and_object_type{
  my ($self, $stable_id, $force) = @_;
 
## Check for Ensembl/Vega style identifiers...
  if( $stable_id =~ /(\w+?)(EST)?([GTPE])\d/ ) {
      return $ensembl_species{$1},$ensembl_type{$3}, $ensembl_db_type{$2} ||'core' if $ensembl_species{$1};
      return $vega_species{$1}, $ensembl_type{$3}, 'vega' if $vega_species{$1};
  }
## Check for Vector base style identifiers
  if( $stable_id =~ /^([A-Z]+)\d+-?(RP)?\w?/ && $vectorbase_species{$1} ) {
      return $vectorbase_species{$1}, $vectorbase_type{$2}||'Gene','core';
  }
  if( $stable_id =~ /^([A-Z]+)?(\.e|E)\d+$/ && ($1 eq '' || $vectorbase_species{$1}) ) {
      return $vectorbase_species{$1||'AGAP'}, 'Exon', 'core';
  }
  return unless defined($force) && $force;
## Finally if "force" is passed see if we can find a non-standard string...
  foreach my $species (qw(saccharomyces_cerevisiae tetraodon_nigroviridis
                          drosophila_melanogaster caenorhabditis_elegans)){
      foreach my $type (qw(Transcript Gene Translation Exon)){
	  my $adaptor = $self->get_adaptor($species, 'core', $type);
	  next unless defined($adaptor);
	  my $entity = $adaptor->fetch_by_stable_id($stable_id);
	  return ucfirst( $species ), $type,'core' if defined $entity;
      }
  }
  return;

}

1;
