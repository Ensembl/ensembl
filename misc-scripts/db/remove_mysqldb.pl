#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


package Script;

use strict;
use warnings;

use Carp;
use DBI;
use Getopt::Long qw/:config no_ignore_case auto_version bundling_override/;
use Pod::Usage;

sub args {
  my ($self) = @_;
  my $opts = { port => 3306, verbose => 0 };
  GetOptions(
    $opts, qw/
      host|hostname|h=s
      port|P=i
      username|user|u=s
      password|pass|p=s
      databases|database|db=s@
      verbose|v
      log=s
      help
      man
      /
  ) or pod2usage(-verbose => 1, -exitval => 1);
  pod2usage(-verbose => 1, -exitval => 0) if $opts->{help};
  pod2usage(-verbose => 2, -exitval => 0) if $opts->{man};
  $self->{opts} = $opts;
  return;
}

my $rcsid = '$Revision$';
our ($VERSION) = $rcsid =~ /(\d+\.\d+)/;

sub run {
  my ($class) = @_;
  my $self = bless({}, $class);
  $self->args();
  $self->logging();
  $self->process();

  if ($self->{oldfh}) {
    select($self->{oldfh});
  }
  
  return;
}

sub defaults {
  my ($self) = @_;
  my $o = $self->opts();
  
  foreach my $required (qw/username databases host/) {
    pod2usage(-msg => sprintf('No -%s specified', $required), -verbose => 2, -exitval => 1);
  }
  #Processing -opt 1 -opt 2,3 into opt => [1,2,3]
  $self->_cmd_line_to_array('databases') if $o->{databases};

  $self->v(q{Working with %d database(s)}, scalar(@{$o->{databases}}));
  
  return;
}

sub logging {
  my ($self) = @_;
  my $o = $self->opts();
  if ($o->{log}) {
    $o->{verbose} = 1;
    my $file = $o->{log};
    open my $fh, '>', $file or die "Cannot open log file '${file}' for writing: $!";
    my $oldfh = select($fh);
    $self->{oldfh} = $oldfh;
  }
  return;
}

sub process {
  my ($self) = @_;
  my $databases = $self->opts()->{databases};
  foreach my $db (@{$databases}) {
    $self->v('Working with database %s', $db);
    $self->dbh($db);
    $self->delete_tables();
    my $ddl = sprintf('drop database %s', $db);
    $self->v(qq{\tRunning: '%s'}, $ddl);
    $self->dbh()->do($ddl);
    $self->clear_tables();
    $self->clear_dbh();
  }
  return;
}

sub delete_tables {
  my ($self) = @_;
  my @tables = keys %{ $self->tables() };
  foreach my $table (@tables) {
    $self->v(q{Processing '%s'}, $table);
    my @sql;
    if($self->is_view($table)) {
      $self->v(q{'%s' is a view; just dropping}, $table);
      push(@sql, 'drop view '.$table);
    }
    else {
      $self->v(q{'%s' is being truncated and then dropped}, $table);
      push(@sql,
        'truncate table '.$table,
        'drop table '.$table,
      );
    }
    
    foreach my $statement (@sql) {
      $self->v(qq{\tRunning: '%s'}, $statement);
      $self->dbh()->do($statement);
    }
  }
  return;
}

sub v {
  my ($self, $msg, @args) = @_;
  return unless $self->opts()->{verbose};
  my $s_msg = sprintf($msg, @args);
  my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) =
    localtime(time());
  print sprintf("[%02d-%02d-%04d %02d:%02d:%02d] %s\n",
                $mday, $mon, $year + 1900,
                $hour, $min, $sec, $s_msg);
  return;
}

sub dbh {
  my ($self, $database) = @_;
  if (!exists $self->{'dbh'}) {
    my $o = $self->opts();
    my %args = (host => $o->{host}, port => $o->{port});
    $args{database} = $database if defined $database;

    my $dsn =
      'DBI:mysql:' . join(q{;}, map { $_ . '=' . $args{$_} } keys %args);
    $self->v('DBI connection URI %s', $dsn);
    my $dbh =
      DBI->connect($dsn, $o->{username}, $o->{password}, { RaiseError => 1 });

    $self->{dbh} = $dbh;
  }
  return $self->{'dbh'};
}

sub clear_dbh {
  my ($self) = @_;
  if (exists $self->{dbh}) {
    $self->{dbh}->disconnect();
    delete $self->{dbh};
  }
  return;
}

sub is_view {
  my ($self, $table) = @_;
  return ($self->tables()->{$table} eq 'VIEW') ? 1 : 0;
}

sub tables {
  my ($self) = @_;
  if (!exists $self->{tables}) {
    my $array =
      $self->dbh()->selectcol_arrayref(
'select TABLE_NAME, TABLE_TYPE from information_schema.TABLES where TABLE_SCHEMA = DATABASE()',
      { Columns => [ 1, 2 ] }
      );
    my %hits = @{$array};
    $self->{tables} = \%hits;
  }
  return $self->{tables};
}

sub clear_tables {
  my ($self) = @_;
  delete $self->{tables};
  return;
}

sub _cmd_line_to_array {
  my ($self, $key) = @_;
  my $array = $self->opts()->{$key};
  $array = (ref($array) && ref($array) eq 'ARRAY') ? $array : [$array];
  my $string = join(q{,}, @{$array});
  my @new_array = split(/,/, $string);
  $self->opts()->{$key} = \@new_array;
  return;
}

Script->run();

1;
__END__

=pod

=head1 NAME

remove_mysqldb.pl

=head1 SYNOPSIS

  #Basic
  ./remove_mysqldb.pl -username USER -password PASS -host HOST [-port PORT] -database DB [-verbose] [-help | -man]
  
  #Advanced
  ./remove_mysqldb.pl -username USER -password PASS -host HOST -port PORT -database DB -database DBTWO -verbose

=head1 DESCRIPTION

A script which is used to drop a database in the most lock friendly way.
If you issue a C<drop database DB> command to MySQL it first aquires
a C<LOCK_open> which is a global mutex always applied for when you open
and close a file. Once the drop command has the lock no other process can
query the MySQL DB until the lock is relinquished once the DB has been dropped.
If your DB is very large or has a large number of tables this can take a long
time and has you performing a denile of service attack on your own DB server.

=head1 OPTIONS

=over 8

=item B<--username | --user | -u>

REQUIRED. Username of the connecting account. Must be able to perform 
drop and truncate calls.

=item B<--password | -pass | -p>

Password of the connecting user.

=item B<--host | --hostname | -h>

REQUIRED. Host name of the database to connect to

=item B<--port | -P>

Optional integer of the database port. Defaults to 3306

=item B<--databases | --database | --db>

Allows database name specification and can be used more than once. Comma 
separated values are supported.

=item B<--verbose>

Makes the program give more information about what is going on. Otherwise
the program is silent.

=item B<--log>

If given the script will write all logs to output. Switches on C<--verbose>

=item B<--help>

Help message

=item B<--man>

Man page

=back

=head1 REQUIREMENTS

=over 8

=item Perl 5.8+

=item DBI

=item DBD::mysql

=back

=cut
