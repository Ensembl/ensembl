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

use File::Spec::Functions;
use IO::File;

sub new
{
  my ($proto, $arg_ref) = @_;

  my $class = ref $proto || $proto;
  my $self =  bless {}, $class;

  $self->host($arg_ref->{host});
  $self->dbname($arg_ref->{dbname});
  $self->user($arg_ref->{user});
  $self->pass($arg_ref->{pass} || '');
  $self->port($arg_ref->{port} || '3306');
  $self->verbose($arg_ref->{verbose});

  return $self;
}

sub verbose {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_verbose} = $arg );
  return $self->{_verbose};
}

sub host {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_host} = $arg );
  return $self->{_host};
}

sub dbname {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_dbname} = $arg );
  return $self->{_dbname};
}

sub user {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_user} = $arg );
  return $self->{_user};
}

sub pass {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_pass} = $arg );
  return $self->{_pass};
}

sub port {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_port} = $arg );
  return $self->{_port};
}

sub dbi
{
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


# Create database if required. Assumes sql/table.sql and sql/populate_metadata.sql
# are present.

sub create {
  my ($self, $sql_dir, $force, $drop_db) = @_;

  my $user   = $self->user;
  my $dbname = $self->dbname;
  my $host   = $self->host;
  my $pass   = $self->pass;
  my $port = $self->port;

  my $dbh = DBI->connect( "DBI:mysql:host=$host:port=$port", $user, $pass,
                          {'RaiseError' => 1});

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

  my $table_file = catfile( $sql_dir, 'sql', 'table.sql' );

  printf( "Creating %s from %s\n", $dbname, $table_file ) if($self->verbose);
  if ( !-e $table_file ) {
    croak( "Cannot open  " . $table_file );
  }

  my $cmd =
    "mysql -u $user -p'$pass' -P $port -h $host $dbname < $table_file";
  system($cmd) == 0
    or croak("Cannot execute the following command (exit $?):\n$cmd\n");

  printf( "Populating metadata in %s from %s\n",
          $dbname, $metadata_file ) if($self->verbose);
  if ( !-e $metadata_file ) {
    croak( "Cannot open " . $metadata_file );
  }

  $cmd = "mysql -u $user -p'$pass' -P $port -h $host "
    . "$dbname < $metadata_file";
  system($cmd) == 0
    or croak("Cannot execute the following command (exit $?):\n$cmd\n");
  return;
}

1;
