package Xref::Test::TestDB;

use Moose;
use File::Temp qw/tempdir/;

extends 'Xref::DB';


has reuse => (
  isa => 'Bool',
  is => 'ro',
  default => 0,
);


sub BUILD {
  my ($self) = @_;
  # auto-configure a DB connection in the absence of user supplied options
  $self->config_file('testdb.conf') if ! defined $self->config_file;
}

override '_validate_config' => sub {
  my ($self,$config) = @_;
  if (! exists $config->{db}) {
    $config->{db} = sprintf '%s_xref_test_%s',$ENV{USER},int(rand(100000));
    $config->{create} = 1;
  }
  super();
};
# On top of default config validator, inject a randomised test DB name


sub DEMOLISH {
  my $self = shift;
  if ($self->reuse == 0) {
    if ( $self->config->{driver} eq 'SQLite') {
      unlink $self->config->{file};
    } elsif ($self->config->{driver} eq 'mysql') {
      $self->schema->storage->dbh->do('drop database '.$self->config->{db});
    }
  }
}

__PACKAGE__->meta->make_immutable;

1;