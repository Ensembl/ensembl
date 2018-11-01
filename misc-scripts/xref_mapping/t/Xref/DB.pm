package Xref::DB;

use Moose;
use Config::General;
use Carp;
use File::Temp qw/tempdir/;
use Xref::Schema;
use DBI;

has schema => (
  isa => 'Xref::Schema',
  is => 'ro',
  lazy=> 1,
  builder => '_init_db'
);

has config_file => (
  isa => 'Str',
  is => 'rw'
);

has config => (
  isa => 'HashRef',
  is => 'ro',
  lazy => 1,
  builder => '_init_config'
);

sub _init_db {
  my $self = shift;
  $self->_validate_config($self->config);
  my %conf = %{ $self->config };
  my %opts;
  $opts{mysql_enable_utf8} = 1 if ($conf{driver} eq 'mysql');
  $opts{mysql_auto_reconnect} = 1 if ($conf{driver} eq 'mysql');
  $opts{sqlite_unicode} = 1 if($conf{driver} eq 'SQLite');
  my $dsn; 
  if ($conf{driver} eq 'SQLite') { 
    $dsn = sprintf("dbi:%s:database=%s",$conf{driver},$conf{file}); 
  } else {
    $dsn = sprintf("dbi:%s:database=%s;host=%s;port=%s",$conf{driver},$conf{db},$conf{host},$conf{port});
  }
  
  my %deploy_opts = ();
  # Example deploy option $deploy_opts{add_drop_table} = 1;
  print STDERR 'Connecting: '.$dsn."\n";
  my $schema = Xref::Schema->connect($dsn, $conf{user},$conf{pass}, \%opts);
  
  if ($conf{create} == 1 && $conf{driver} eq 'mysql') {
    my $dbh = DBI->connect(sprintf("DBI:%s:database=;host=%s;port=%s",$conf{driver},$conf{host},$conf{port}),$conf{user},$conf{pass}, \%opts);
    $dbh->do('CREATE DATABASE '.$conf{db}.';');
    $dbh->disconnect;
  }
  $schema->deploy(\%deploy_opts) if $conf{create} == 1;

  return $schema;
}

sub _init_config {
  my $self = shift;
  if (! $self->config_file) { confess 'No config or config_file provided to new(). Cannot execute'; }
  my $conf = Config::General->new($self->config_file);
  my %opts = $conf->getall();
  return \%opts;
}

sub _validate_config {
  my ($self,$config) = @_;
  my @required_keys = qw/driver/;
  if ($config->{driver} eq 'mysql') {
    push @required_keys, qw/db host port user pass/;
  } elsif ($config->{driver} eq 'SQLite') {
    push @required_keys, qw/file/;
  } else {
    confess "TestDB config requires parameter 'driver' with value mysql or SQLite";
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
}

# Shortcut for accessing a database handle directly. I get the impression we might be doing this a lot.
sub dbh {
  my $self = shift;
  return $self->schema->storage->dbh;
}

# Shortcut for creating things on the fly
sub create_db_row {
  my ($self,$model, $params) = @_;
  my $source = $self->schema->resultset($model)->create(
    $params
  );
  return $source;
}

__PACKAGE__->meta->make_immutable;

1;