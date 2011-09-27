package XrefMapper::uniparc;

use base qw(XrefMapper::db);

sub method {
  my ($self, $method) = @_;
  $self->{method} = $method if defined $method;
  return $self->{method};
}

1;