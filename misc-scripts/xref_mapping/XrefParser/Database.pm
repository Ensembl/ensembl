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

package XrefParser::Database;

use strict;
use warnings;
use Carp;
use DBI;
use Bio::EnsEMBL::DBSQL::DBConnection;
use File::Spec::Functions;
use IO::File;
use English;

sub new {
  my ($proto, $arg_ref) = @_;

  my $class = ref $proto || $proto;
  my $self =  bless {}, $class;
  
  $self->dbc(ref $arg_ref eq 'Bio::EnsEMBL::DBSQL::DBConnection' 
    ? $arg_ref : Bio::EnsEMBL::DBSQL::DBConnection->new(
    -HOST => $arg_ref->{host},
    -DBNAME => $arg_ref->{dbname},
    -USER => $arg_ref->{user},
    -PASS => $arg_ref->{pass} || '',
    -PORT => $arg_ref->{port} || '3306'
  ));
  $self->verbose($arg_ref->{verbose});

  return $self;
}

sub verbose {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_verbose} = $arg );
  return $self->{_verbose};
}

sub dbc {
  my $self = shift;

  if(@_) {
    my $arg = shift;

    if(defined($arg)) {
      croak "$arg is not a DBConnection" unless $arg->isa('Bio::EnsEMBL::DBSQL::DBConnection');
      $self->{_dbc} = $arg;
    }
  }

  return $self->{_dbc};
}

sub host {
  my ($self, $arg) = @_;
  $self->dbc->host($arg) if defined $arg;
  return $self->dbc->host;
}

sub dbname {
  my ($self, $arg) = @_;
  $self->dbc->dbname($arg) if defined $arg;
  return $self->dbc->dbname;
}

sub user {
  my ($self, $arg) = @_;
  $self->dbc->user($arg) if defined $arg;
  return $self->dbc->user;
}

sub pass {
  my ($self, $arg) = @_;
  $self->dbc->pass($arg) if defined $arg;
  return $self->dbc->pass;
}

sub port {
  my ($self, $arg) = @_;
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
            { 'RaiseError' => 1 } )
          or croak( "Can't connect to database: " . $DBI::errstr );
        $dbi->{'mysql_auto_reconnect'} = 1; # Reconnect on timeout
    }
    return $dbi;
}

# Create database if required. 
# Assumes sql/table.sql and sql/populate_metadata.sql are present.
sub create {
  my ($self, $sql_dir, $force, $drop_db) = @_;
  $self->recreate_database($force,$drop_db);
  $self->populate($sql_dir, $force);
}

sub populate {
  my ($self, $sql_dir, $force) = @_;
  my $table_file = catfile( $sql_dir, 'sql', 'table.sql' );
  my $metadata_file = $self->prepare_metadata_file($sql_dir, $force);
  $self->populate_with_file($table_file);
  $self->populate_with_file($metadata_file);
}

sub prepare_metadata_file {
  my ($self, $sql_dir, $force) = @_;
  my $metadata_file =
    catfile( $sql_dir, 'sql', 'populate_metadata.sql' );
  my $ini_file = catfile( $sql_dir, 'xref_config.ini' );

  local $| = 1;    # flush stdout

  # Figure out whether to run 'xref_config2sql.pl' or not by comparing
  # the timestamps on 'xref_config.ini' and 'sql/populate_metadata.sql'.
  my $ini_tm  = ( stat $ini_file )[9];
  my $meta_tm = ( stat $metadata_file )[9];

  if ( !defined($meta_tm) || $ini_tm > $meta_tm ) {
    my $reply;
    if($force){
      $reply = 'y';
    }
    else{
      printf( "==> Your copy of 'xref_config.ini' is newer than '%s'\n",
	      catfile( 'sql', 'populate_metadata.sql' ) );
      print("==> Should I re-run 'xref_config2sql.pl' for you? [y/N]: ");

      $reply = <ARGV>;
      chomp $reply;
    }
    if ( lc( substr( $reply, 0, 1 ) ) eq 'y' ) {
      my $cmd = sprintf( "perl %s %s >%s",
                         catfile( $sql_dir, 'xref_config2sql.pl' ),
                         $ini_file, $metadata_file );

      if ( system($cmd) == 0 ) {
        print("==> Done.\n") if($self->verbose);
      } else {
        if ( $? == -1 ) {
          croak("Failed to execute: $!\n");
        } elsif ( $? & 127 ) {
          croak(
                 sprintf( "Command died with signal %d, %s coredump\n",
                          ( $? & 127 ),
                          ( $? & 128 ) ? 'with' : 'without'
                 ) );
        } else {
          croak( sprintf( "Command exited with value %d\n", $? >> 8 ) );
        }
      }

    }
  } ## end if ( !defined($meta_tm...
  return $metadata_file;
}

sub recreate_database {
  my ($self,$force, $drop_db) = @_;
  my $user   = $self->user;
  my $dbname = $self->dbname;
  my $host   = $self->host;
  my $pass   = $self->pass;
  my $port = $self->port;
  my $dbh = DBI->connect( "DBI:mysql:host=$host:port=$port", $user, $pass,
                          {'RaiseError' => 1});
  # check to see if the database already exists
  my %dbs = map {$_->[0] => 1} @{$dbh->selectall_arrayref('SHOW DATABASES')};

  if ($dbs{$dbname}) {

    if ( $drop_db ) {
	$dbh->do( "DROP DATABASE $dbname" );
	print "Database $dbname dropped\n" if($self->verbose) ;
    } else {

      my $p;
      if($force){
	$p = 'yes';
      }
      else{
	print "WARNING: about to drop database $dbname on $host:$port; yes to confirm, otherwise exit: ";
	$p = <ARGV>;
      }
      chomp $p;
      if ($p eq 'yes') {
	$dbh->do( "DROP DATABASE $dbname" );
	print "Removed existing database $dbname\n" if($self->verbose);
      } else {
	print "$dbname NOT removed\n";
	exit(1);
      }

    }
  }

  $dbh->do( 'CREATE DATABASE ' . $dbname );
}

sub populate_with_file {
  my ($self, $sql_file) = @_;
  my $previous_input_record_separator = $INPUT_RECORD_SEPARATOR;
  $INPUT_RECORD_SEPARATOR = ";";
  open(my $sql_fh, "<", $sql_file) or die $sql_file;
  while(<$sql_fh>){
    s/#(.*?)\n//g;
    next if /^\s+$/;
    $self->dbc->do($_);
  }
  $INPUT_RECORD_SEPARATOR = $previous_input_record_separator;
}

1;
