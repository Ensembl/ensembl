=head1 LICENSE

See the NOTICE file distributed with this work for additional
information regarding copyright ownership.

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



package XrefParser::UniProtParser::Loader;

use strict;
use warnings;

use Carp;

# FIXME: for testing
use Data::Dumper;



sub new {
  my ( $proto, $arg_ref ) = @_;

  my $self = {
              'batch_size'  => $arg_ref->{'batch_size'} // 1,
              'dbh'         => $arg_ref->{'dbh'},
            };
  my $class = ref $proto || $proto;
  bless $self, $class;
  $self->_clear_send_buffer();

  return $self;
}


sub load {
  my ( $self, $transformed_data ) = @_;

  if ( ! defined $transformed_data ) {
    return;
  }

  # FIXME: pass $count from Transformer
  $self->_add_to_send_buffer( $transformed_data );

  if ( $self->{'send_backlog'} >= $self->{'batch_size'} ) {
    $self->flush();
  }

  return;
}


sub flush {
  my ( $self ) = @_;

  # FIXME: call BaseParser::upload_xref_object_graphs($buffer, $dbh)
  print "SENDING: " . Dumper( $self->{'send_buffer'} );

  $self->_clear_send_buffer();

  return;
}



sub _add_to_send_buffer {
  my ( $self, $entry ) = @_;

  push @{ $self->{'send_buffer'} }, $entry;
  # FIXME: update the counter CORRECTLY
  $self->{'send_backlog'}++;

  return;
}

sub _clear_send_buffer {
  my ( $self ) = @_;

  $self->{'send_buffer'}  = [];
  $self->{'send_backlog'} = 0;

  return;
}


1;
