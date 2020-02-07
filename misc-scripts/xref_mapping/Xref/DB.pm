=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

=head1 DESCRIPTION

A database class that auto-instantiates a DBIC schema

Convenience methods prevent having to delve into DBIC guts for common activities

=head1 SYNOPSIS

my $db = Xref::DB->new(
  config => {
    host => 'db.com',
    port => 3306,
    user => 'me',
    pass => 'Go on then',
    driver => 'mysql',
    db => 'name_for_db',
    create => 1, # Deploys the schema to the DB for first use
  }
);

$db = Xref::DB->new(
  config_file => 'db.conf' # config options in Config::General format
);

my $dbh = $db->dbh; # $dbh is a DBI database handle borrowed for direct SQL
$dbh->prepare('DROP TABLE dependent_xref');

$db->create_db_row('Xref',{
  xref_id => 1,
  accession => 'YAY',
  description => 'Sample new Xref',
  source_id => 1,
  ...
});

=cut


package Xref::DB;

use strict;
use warnings;

use Moose;
use namespace::autoclean;
use Config::General;
use Config::IniFiles; # FIXME normalise around one config format
use Carp;
use Xref::Schema;
use DBI;

has schema => (
  isa => 'Xref::Schema',
  is => 'ro',
  builder => '_init_db'
);

has config_file => (
  isa => 'Str',
  is => 'rw',
  builder => '_guess_config'
);

has config => (
  isa => 'HashRef',
  is => 'rw'
);

has now_function => (
  isa => 'Str',
  default => 'now()',
  is => 'rw',
);


=head2 _init_db
  Arg [1]    : HashRef of configuation parameters (driver, db, host, port, user, pass)
  Description: Initialise the core database.
  Return type: schema
  Caller     : internal

=cut

sub _init_db {
  my $self = shift;

  $self->_init_config if ! defined $self->config;
  $self->_validate_config($self->config);
  my %conf = %{ $self->config };
  my $enable_unicode = $conf{enable_unicode} // 0;
  my %opts;
  $opts{mysql_enable_utf8} = $enable_unicode if ($conf{driver} eq 'mysql');
  $opts{mysql_auto_reconnect} = 1 if ($conf{driver} eq 'mysql');
  $opts{sqlite_unicode} = $enable_unicode if($conf{driver} eq 'SQLite');
  my $dsn;
  if ($conf{driver} eq 'SQLite') {
    $dsn = sprintf 'dbi:%s:database=%s',$conf{driver},$conf{file};
    $self->now_function("date('now')");
  } else {
    $dsn = sprintf 'dbi:%s:database=%s;host=%s;port=%s', $conf{driver}, $conf{db}, $conf{host}, $conf{port};
  }

  my %deploy_opts = ();
  # Example deploy option $deploy_opts{add_drop_table} = 1;
  my $schema = Xref::Schema->connect($dsn, $conf{user}, $conf{pass}, \%opts);

  if ($conf{create} == 1 && $conf{driver} eq 'mysql') {
    my $dbh = DBI->connect(
      sprintf('DBI:%s:database=;host=%s;port=%s', $conf{driver}, $conf{host}, $conf{port}), $conf{user}, $conf{pass}, \%opts
    );

    # Remove database if already exists
    my %dbs = map {$_->[0] => 1} @{$dbh->selectall_arrayref('SHOW DATABASES')};
    my $dbname = $conf{db};
    if ($dbs{$dbname}) {
      $dbh->do( "DROP DATABASE $dbname;" );
    }

    my $db_collation = ( $enable_unicode ) ? 'utf8_general_ci' : 'latin1_swedish_ci';
    $dbh->do("CREATE DATABASE $dbname DEFAULT COLLATE ${db_collation};");

    $dbh->disconnect;
  }

  $schema->deploy(\%deploy_opts) if $conf{create} == 1;

  return $schema;
} ## end sub _init_db


=head2 _guess_config
  Description: Don't want production use to guess at least at the moment.
               This mainly exists so TestDB can override and replace with a
               useful default
  Return type: undef
  Caller     : internal

=cut

sub _guess_config {
  return;
} ## end sub _guess_config


=head2 _init_config
  Arg [1]    : HashRef of configuation parameters (driver, db, host, port, user, pass)
  Description: Initialisae the loading of the configuration file.
  Return type: HashRef - $self->config
  Caller     : internal

=cut

sub _init_config {
  my $self = shift;

  if (defined $self->config_file) {
    my $conf = Config::General->new($self->config_file);
    my %opts = $conf->getall();
    $self->config(\%opts);
  } else {
    confess 'No config or config_file provided to new(). Cannot execute';
  }

  return $self->config;
} ## end sub _init_config


=head2 _validate_config
  Arg [1]    : HashRef of configuation parameters (driver, db, host, port, user, pass)
  Description: Configuration file parameter validation
  Return type: DBI database handle
  Caller     : internal

=cut

sub _validate_config {
  my ($self,$config) = @_;
  my @required_keys = qw/driver/;
  if ($config->{driver} eq 'mysql') {
    push @required_keys, qw/db host port user pass/;
  } elsif ($config->{driver} eq 'SQLite') {
    push @required_keys, qw/file/;
  } else {
    confess q(TestDB config requires parameter 'driver' with value mysql or SQLite);
  }
  my @errors;
  foreach my $constraint (@required_keys) {
    if (! exists $config->{$constraint}) {
      push @errors, "Missing argument '$constraint'";
    }
  }
  if (scalar @errors > 0) {
    confess sprintf "%s \n%s",
      ($self->config_file) ? 'Missing options in '.$self->config_file. ': ' : 'Missing options in supplied config: ',
      join ';',@errors;
  }
} ## end sub _validate_config


=head2 dbh
  Description: Shortcut for accessing a database handle directly. I get the
               impression we might be doing this a lot.
  Return type: DBI database handle
  Caller     : internal

=cut

sub dbh {
  my $self = shift;
  return $self->schema->storage->dbh;
} ## end sub dbh


=head2 create_db_row
  Arg [1]    : model
  Arg [2]    : arguments : These should be key-value pairs matching the rows in
                           the table
  Description: Shortcut for creating things on the fly
  Return type:
  Caller     : internal

=cut

sub create_db_row {
  my ($self,$model, $params) = @_;
  my $source = $self->schema->resultset($model)->create(
    $params
  );
  return $source;
} ## end sub create_db_row


=head2 populate_metadata

  Arg [1]    : Config file path, normally xref_config.ini
  Description: Loads species and source information into the schema from the
               supplied file in Config::IniFiles format
  Caller     : User

=cut

# TODO: Provide species AND division in order to limit quantity of madness

sub populate_metadata {
  my ($self, $config_path) = @_;

  my $config = $self->_load_xref_config($config_path);

  # Populate species table with species taxa
  print "Iterating over species groups\n";

  my %sources;
  # First build up records for each potential source
  foreach my $section ( $config->GroupMembers('source')) {

    my ( $source_name ) = $section =~ /
      \A
      source\s+(.*)
      \Z
    /x;
    $sources{$source_name} = $self->_mangle_source_block($section, $config);

  }

  my %compiled_config;
  # Next, build a config hash from the species entries in the config file
  # and populate them with relevant source information from the %sources
  foreach my $section ( $config->GroupMembers('species') ) {
    my ( $species_name ) = $section =~ /
      \A
      species\s+(\S+)
      \s*
      \Z
    /x;

    my $taxon_id = $config->val( $section, 'taxonomy_id' );

    my @source_names = $config->val( $section, 'source', ());
    my %sources_for_species;

    foreach my $source_entry (@source_names) {

      if (! exists $sources{$source_entry}) {
        confess 'Species config references a source that is not defined in the config: '.$source_entry;
      }

      $sources_for_species{$source_entry} = $sources{$source_entry};
    }

    $compiled_config{$species_name} = {
      species_id => $taxon_id, # species_id == taxon_id in xref system,
      taxonomy_id => $taxon_id,
      sources => \%sources_for_species
    }
  }

  # Now populate the database with the result each source entry in the config with source-specific info

  foreach my $species ( keys %compiled_config ) {

      my $species_record = $self->schema->resultset('Species')->create({
        species_id => $compiled_config{$species}{species_id},
        name => $species,
        taxonomy_id => $compiled_config{$species}{taxonomy_id}
      });

      my $compiled_sources = $compiled_config{$species}{sources};

      foreach my $source ( keys %$compiled_sources ) {
        my $parser = delete $compiled_sources->{$source}->{parser};
        my $source_record = $self->schema->resultset('Source')->find_or_create(
          $compiled_sources->{$source}  # once trimmed, we can pump the source hash straight into DBIC
        );
        delete $sources{$source};

        $source_record->create_related(
           'source_url',
           {
             species_id => $compiled_config{$species}{species_id},
             parser => $parser
           }
         );
      } # End foreach source


  } # End foreach species

  # Add any sources which are not explicitly linked to a species

  foreach my $source ( keys %sources ) {

    my $parser = delete $sources{$source}->{parser};
    my $source_record = $self->schema->resultset('Source')->find_or_create(
      $sources{$source}
    );

  }

  return;
}

=head2 _load_xref_config

Arg [1]    : path to xref_config.ini
Description: Load the config file and sanity-check to ensure content is
             properly formatted for loading.
             FIXME: make parsing into a parser? Better yet switch to a format
             which can express this grammar correctly without micro-formatting
Returntype : Hashref config, as returned by Config::IniFiles

=cut

sub _load_xref_config {
  my ($self, $config_path) = @_;

  if (! -e $config_path) {
    confess "Unable to open config file $config_path";
  }
  my $config = Config::IniFiles->new(-file => $config_path);

  if (! defined $config) {
    confess "Errors in $config_path, unable to parse: ". join ',', @Config::IniFiles::errors;
  }
  my $source_id = 0;
  my %source_ids; # A tracker for sections we've seen

  foreach my $section ( $config->GroupMembers('source')) {
    my ( $spaces, $source_name ) = $section =~ /
      \A
      source(\s+)(\S+)
      \s*
      \Z
    /x;

    # This validation is how we used to validate this file
    if ( length($spaces) > 1 ) {
      confess(
        sprintf(
          "Too many spaces between the words 'source' and '%s'\nwhile reading source section '[%s]'\n",
          $source_name,
          $section
        )
      );
    }

    if ( index( $config->val( $section, 'name' ), "\n" ) != -1 ) {
      confess(
        sprintf( "The source section '[%s]' occurs more than once in \n", $section )
      );
    }

    $source_ids{$section} = ++$source_id; # Record existence of species->source
  }

  foreach my $section ( $config->GroupMembers('species')) {
    my ( $spaces, $species_name ) = $section =~ /
      \A
      species(\s+)(\S+)
      \s*
      \Z
    /x;

    if ( length($spaces) > 1 ) {
      confess(
        sprintf(
          "Too many spaces between the words 'species' and '%s'\nwhile reading species section '[%s]'\n",
          $species_name,
          $section
        )
      );
    }

    foreach my $source_name (split qr{ \n }msx, $config->val($section, 'source') ) {

      $source_name =~ s{ \s\z }{}msx; # Config file can easily contain trailling whitespace
      my $source_section = "source $source_name";
      # Check integrity of species to source mentions
      if ( !exists $source_ids{$source_section} ) {
        confess(
          sprintf( "Can not find source section '[%s]'\nwhile reading species section '[%s]'\n",
                   $source_section,
                   $section
                 )
        );
      }

    }

  }
  return $config;
}

=head2 _mangle_source_block

Arg 1      : String - section header, the name for the section that identifies
                      the section we need to extract values from
Arg 2      : Hashref - config file content, the output of _load_xref_config()
Description: Takes a single source record from xref_config.ini and creates a
             hashref. Some arguments are massaged to get them into the schema
             easily
Returntype : Hashref of a single source's properties
Caller     : Internal

=cut

sub _mangle_source_block {
  my ($self, $section, $config) = @_;

  my $source_config;
  # Get the easy ones done
  foreach my $key (qw/name priority parser/) {
    my $value = $config->val($section, $key);
    if (defined $value) {
      $source_config->{$key} = $value;
    }
  }

  my $priority_description = $config->val($section, 'prio_descr');
  if (defined $priority_description) {
    $source_config->{priority_description} = $priority_description;
  }

  return $source_config;
}

__PACKAGE__->meta->make_immutable;

1;
