package Xref::Schema::ResultSet::DependentXref;
use strict;
use warnings;

use parent 'DBIx::Class::ResultSet';

sub fetch_dependent_xref {
  my ($self,$direct_accession,$dependent_accession) = @_;

  my $hit = $self->find({
      'master_xref.accession' => $direct_accession,
      'dependent_xref.accession' => $dependent_accession
    }, {
      join => [ 'master_xref','dependent_xref']
    });
  
  return $hit;
}

1;