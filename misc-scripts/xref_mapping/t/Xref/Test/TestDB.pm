package Xref::Test::TestDB;

use Moose;
use File::Temp qw/tempdir/;

extends 'Xref::DB';


has reuse => (
  isa => 'Bool',
  is => 'ro',
  default => 0,
);

sub DEMOLISH {
  my $self = shift;
  if ($self->reuse == 0) {
    if ( $self->config->{driver} eq 'SQLite') {
      unlink $self->config->{'file'};
    } elsif ($self->config->{driver} eq 'mysql') {
      $self->schema->undeploy;
    }
  }
}

__PACKAGE__->meta->make_immutable;

1;