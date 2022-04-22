
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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

package XrefParser::Database;

use strict;
use warnings;
use Carp;
use DBI;
use Bio::EnsEMBL::DBSQL::DBConnection;
use File::Spec::Functions;
use IO::File;
use English;
use vars qw(@ISA);
use strict;

@ISA = qw( Bio::EnsEMBL::DBSQL::DBConnection);

sub new {
  my ( $proto, $arg_ref ) = @_;

  my $class = ref $proto || $proto;
  my $self = bless {}, $class;

  $self->dbc(
        ref $arg_ref eq 'Bio::EnsEMBL::DBSQL::DBConnection' ? $arg_ref :
          Bio::EnsEMBL::DBSQL::DBConnection->new(
                                     -HOST   => $arg_ref->{host},
                                     -DBNAME => $arg_ref->{dbname},
                                     -USER   => $arg_ref->{user},
                                     -PASS   => $arg_ref->{pass} || '',
                                     -PORT => $arg_ref->{port} || '3306'
          ) );
  $self->verbose( $arg_ref->{verbose} );

  return $self;
}

sub verbose {
  my ( $self, $arg ) = @_;

  ( defined $arg ) && ( $self->{_verbose} = $arg );
  return $self->{_verbose};
}

sub dbc {
  my $self = shift;

  if (@_) {
    my $arg = shift;

    if ( defined($arg) ) {
      croak "$arg is not a DBConnection"
        unless $arg->isa('Bio::EnsEMBL::DBSQL::DBConnection');
      $self->{_dbc} = $arg;
    }
  }

  return $self->{_dbc};
}

sub host {
  my ( $self, $arg ) = @_;
  $self->dbc->host($arg) if defined $arg;
  return $self->dbc->host;
}

sub dbname {
  my ( $self, $arg ) = @_;
  $self->dbc->dbname($arg) if defined $arg;
  return $self->dbc->dbname;
}

sub user {
  my ( $self, $arg ) = @_;
  $self->dbc->user($arg) if defined $arg;
  return $self->dbc->user;
}

sub pass {
  my ( $self, $arg ) = @_;
  $self->dbc->pass($arg) if defined $arg;
  return $self->dbc->pass;
}

sub port {
  my ( $self, $arg ) = @_;
  $self->dbc->port($arg) if defined $arg;
  return $self->dbc->port;
}

sub dbi {
  my $self = shift;
  my $dbi;

  if ( !defined $dbi || !$dbi->ping() ) {
    my $connect_string =
      sprintf( "dbi:mysql:host=%s;port=%s;database=%s",
               $self->host, $self->port, $self->dbname );

    $dbi =
      DBI->connect( $connect_string, $self->user, $self->pass,
                    { 'RaiseError' => 1 } ) or
      croak( "Can't connect to database: " . $DBI::errstr );
    $dbi->{'mysql_auto_reconnect'} = 1;    # Reconnect on timeout
  }
  return $dbi;
}

# Create database if required.
# Assumes sql/table.sql and sql/populate_metadata.sql are present.
sub create {
  my ( $self, $sql_dir, $force, $drop_db, $preparse, $xref_source_dbi) = @_;
  $self->recreate_database( $force, $drop_db );
  $self->populate( $sql_dir, $force, $preparse, $xref_source_dbi);
}

sub populate {
  my ( $self, $sql_dir, $force, $preparse, $xref_source_dbi ) = @_;
  my $table_file = catfile( $sql_dir, 'sql', 'table.sql' );
  my $metadata_file;
  if ($xref_source_dbi) {
    $metadata_file = $self->copy_source_metadata_file($xref_source_dbi);
  } else {
    $metadata_file = $self->prepare_metadata_file( $sql_dir, $force, $preparse );
  }
  $self->populate_with_file($table_file);
  $self->populate_with_file($metadata_file);
}

sub copy_source_metadata_file {
  my ( $self, $xref_source_dbi) = @_;

  my $metadata_file = "metadata_config.sql";
  open(CONFIG, ">$metadata_file") or die "Can't open $metadata_file: $!";

  my $select_species_sth = $xref_source_dbi->prepare("SELECT species_id, taxonomy_id, name, aliases FROM species;");
  my ($species_id, $taxonomy_id, $name, $aliases);
  $select_species_sth->execute();
  $select_species_sth->bind_columns(\$species_id, \$taxonomy_id, \$name, \$aliases);
  while ($select_species_sth->fetch()) {
    print CONFIG "INSERT INTO species (species_id, taxonomy_id, name, aliases) ".
    " VALUES ('$species_id', '$taxonomy_id', '$name', '$aliases');";
    print CONFIG "\n\n";
  }
  $select_species_sth->finish();

  my $select_source_sth = $xref_source_dbi->prepare("SELECT source_id, name, source_release, ordered, priority, priority_description, status FROM source");
  my ($source_id, $source_release, $ordered, $priority, $priority_description, $status);
  $select_source_sth->execute();
  $select_source_sth->bind_columns(\$source_id, \$name, \$source_release, \$ordered, \$priority, \$priority_description, \$status);
  while ($select_source_sth->fetch()) {
    print CONFIG 'INSERT INTO source (source_id, name, source_release, ordered, priority, priority_description, status) ';
    printf CONFIG (' VALUES (%d, "%s", "%s", %d, %d, "%s", "%s") ;', 
    	    $source_id, $name, $source_release, $ordered, $priority, $priority_description, $status);
    print CONFIG "\n\n";
  }
  $select_source_sth->finish();

  my $select_dependent_sth = $xref_source_dbi->prepare("SELECT master_source_id, dependent_name FROM dependent_source");
  my ($master_source_id, $dependent_name);
  $select_dependent_sth->execute();
  $select_dependent_sth->bind_columns(\$master_source_id, \$dependent_name);
  while ($select_dependent_sth->fetch()) {
    print CONFIG "INSERT IGNORE INTO dependent_source (master_source_id, dependent_name) ".
    " VALUES ('$master_source_id', '$dependent_name');";
    print CONFIG "\n\n";
  }
  $select_dependent_sth->finish();

  my $select_source_url_sth = $xref_source_dbi->prepare("SELECT source_url_id, source_id, species_id, parser FROM source_url");
  my ($source_url_id, $parser);
  $select_source_url_sth->execute();
  $select_source_url_sth->bind_columns(\$source_url_id, \$source_id, \$species_id, \$parser);
  while ($select_source_url_sth->fetch()) {
    print CONFIG "INSERT INTO source_url (source_url_id, source_id, species_id, parser) ".
    " VALUES ('$source_url_id', '$source_id', '$species_id', '$parser');";
    print CONFIG "\n\n";
  }
  $select_source_url_sth->finish();

  close (CONFIG);
  return $metadata_file;
}

sub prepare_metadata_file {
  my ( $self, $sql_dir, $force, $preparse ) = @_;
  my $metadata_file =
    catfile( $sql_dir, 'sql', 'populate_metadata.sql' );
  my $ini_file = catfile( $sql_dir, 'xref_config.ini' );

  local $| = 1;    # flush stdout

  # Figure out whether to run 'xref_config2sql.pl' or not by comparing
  # the timestamps on 'xref_config.ini' and 'sql/populate_metadata.sql'.
  my $ini_tm  = ( stat $ini_file )[9];
  my $meta_tm = ( stat $metadata_file )[9];

  $preparse = 0 if (!defined($preparse));

  if ( !defined($meta_tm) || $ini_tm > $meta_tm ) {
    my $reply;
    if ($force) {
      $reply = 'y';
    }
    else {
      printf( "==> Your copy of 'xref_config.ini' is newer than '%s'\n",
              catfile( 'sql', 'populate_metadata.sql' ) );
      print(
           "==> Should I re-run 'xref_config2sql.pl' for you? [y/N]: ");

      $reply = <ARGV>;
      chomp $reply;
    }
    if ( lc( substr( $reply, 0, 1 ) ) eq 'y' ) {
      my $cmd = sprintf( "perl %s %s %u >%s",
                         catfile( $sql_dir, 'xref_config2sql.pl' ),
                         $ini_file, $preparse, $metadata_file );

      if ( system($cmd) == 0 ) {
        print("==> Done.\n") if ( $self->verbose );
      }
      else {
        if ( $? == -1 ) {
          croak("Failed to execute: $!\n");
        }
        elsif ( $? & 127 ) {
          croak( sprintf( "Command died with signal %d, %s coredump\n",
                          ( $? & 127 ),
                          ( $? & 128 ) ? 'with' : 'without' ) );
        }
        else {
          croak( sprintf( "Command exited with value %d\n", $? >> 8 ) );
        }
      }

    }
  } ## end if ( !defined($meta_tm...))
  return $metadata_file;
} ## end sub prepare_metadata_file

sub recreate_database {
  my ( $self, $force, $drop_db ) = @_;
  my $user   = $self->user;
  my $dbname = $self->dbname;
  my $host   = $self->host;
  my $pass   = $self->pass;
  my $port   = $self->port;
  my $dbh    = DBI->connect( "DBI:mysql:host=$host:port=$port",
                          $user, $pass, { 'RaiseError' => 1 } );
  # check to see if the database already exists
  my %dbs = map { $_->[0] => 1 }
    @{ $dbh->selectall_arrayref('SHOW DATABASES') };

  if ( $dbs{$dbname} ) {

    if ($drop_db) {
      $dbh->do("DROP DATABASE $dbname");
      print "Database $dbname dropped\n" if ( $self->verbose );
    }
    else {

      my $p;
      if ($force) {
        $p = 'yes';
      }
      else {
        print
"WARNING: about to drop database $dbname on $host:$port; yes to confirm, otherwise exit: ";
        $p = <ARGV>;
      }
      chomp $p;
      if ( $p eq 'yes' ) {
        $dbh->do("DROP DATABASE $dbname");
        print "Removed existing database $dbname\n"
          if ( $self->verbose );
      }
      else {
        print "$dbname NOT removed\n";
        exit(1);
      }

    }
  } ## end if ( $dbs{$dbname} )

  $dbh->do( 'CREATE DATABASE ' . $dbname );
} ## end sub recreate_database

sub populate_with_file {
  my ( $self, $sql_file ) = @_;
  my $previous_input_record_separator = $INPUT_RECORD_SEPARATOR;
  $INPUT_RECORD_SEPARATOR = ";";
  open( my $sql_fh, "<", $sql_file ) or die $sql_file;
  while (<$sql_fh>) {
    s/#(.*?)\n//g;
    next if /^\s+$/;
    $self->dbc->do($_);
  }
  $INPUT_RECORD_SEPARATOR = $previous_input_record_separator;
}

1;
