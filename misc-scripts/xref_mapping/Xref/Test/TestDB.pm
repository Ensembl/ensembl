=head1 LICENSE

See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.
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

A database class that creates the DBIC schema in a test database, and then
cleans up after itself.

Convenience methods prevent having to delve into DBIC guts for common activities

=head1 SYNOPSIS

my $testdb = Xref::Test::TestDB->new(
  config_file => 'filename_if_not_default.conf',
  reuse => 1 # This prevents the test DB cleanup, so you can debug database content
);

=cut

package Xref::Test::TestDB;

use strict;
use warnings;
use Moose;
use namespace::autoclean;

extends 'Xref::DB';

has reuse => (
  isa => 'Bool',
  is => 'ro',
  default => 0,
);

sub _guess_config {
  return 'testdb.conf';
}

# On top of default config validator, inject a randomised test DB name
around '_init_config' => sub {
  my ($sub, $self, @args) = @_;

  my $proto_config = $self->$sub();

  if (! exists $proto_config->{db}) {
    $proto_config->{db} = sprintf '%s_xref_test_%s',$ENV{USER},int(rand(100000));
    $proto_config->{create} = 1;
  }
  # return $proto_config;
  $self->config($proto_config);
  return;
};


=head2 DEMOLISH

Description: It's a destructor. It cleans up databases left behind by the test
             Behaviour is overridden with $self->reuse(1)

=cut
sub DEMOLISH {
  my $self = shift;
  if ($self->reuse == 0 && defined $self->config) {
    if ( $self->config->{driver} eq 'SQLite') {
      unlink $self->config->{file};
    } elsif ($self->config->{driver} eq 'mysql') {
      $self->schema->storage->dbh->do('drop database '.$self->config->{db});
    }
  }
  return;
}

__PACKAGE__->meta->make_immutable;

1;
