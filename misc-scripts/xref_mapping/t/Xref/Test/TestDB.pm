package Xref::Test::TestDB;

use Moose;
use Config::General;
use Carp;
use File::Temp qw/tempdir/;
use Xref::Schema;

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

has reuse => (
  isa => 'Bool',
  is => 'ro',
  default => 0,
);

sub _init_db {
  my $self = shift;
  my %conf = %{ $self->config };
  $self->_validate_config(\%conf);
  my %opts;
  $opts{mysql_enable_utf8} = 1 if ($conf{driver} eq 'mysql');
  $opts{sqlite_unicode} = 1 if($conf{driver} eq 'SQLite');
  my $dsn; 
  if ($conf{driver} eq 'SQLite') { 
    $dsn = sprintf("dbi:%s:database=%s",$conf{driver},$conf{file}); 
  } else {
    $dsn = sprintf("dbi:%s:database=%s;host=%s;port=%s",$conf{driver},$conf{db},$conf{host},$conf{port}); 
  }
  print STDERR 'Connecting: '.$dsn."\n";

  my $schema = Xref::Schema->connect($dsn, $conf{user},$conf{pass}, \%opts);
  my %deploy_opts = ();
  $deploy_opts{add_drop_table} = 1 if $conf{driver} eq 'mysql';
  $schema->deploy(%deploy_opts);

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

sub DEMOLISH {
  my $self = shift;
  if ($self->reuse == 0 && $self->config->{driver} eq 'SQLite') {
    unlink $self->config->{'file'};
  }
}

__PACKAGE__->meta->make_immutable;

1;