package Xref::Schema::ResultSet::Xref;
use strict;
use warnings;

use parent 'DBIx::Class::ResultSet';

sub check_direct_xref {
  my ($self,$params) = @_;

  my $hit = $self->find($params);
  # {
  #   accession => $params->{accession},
  #   label => $params->{display_label},
  #   description => $params->{description}
  # }
  return 1 if defined $hit;
  return;
}

1;