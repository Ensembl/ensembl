=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

package Xref::Schema::ResultSet::Xref;
use strict;
use warnings;

use parent 'DBIx::Class::ResultSet';


=head2 check_direct_xref

Arg [1]    : Hashref of constraints for the xref, e.g. accession, label,
             info_type etc. Can be any column in the schema
Description: A query wrapper to reduce the call stack when checking a single
             Xref is in the database
Returntype : Boolean - True means the Xref was in the database

Example    : $db->schema->resultset('Xref')->check_direct_xref({
                accession => 'BOB'
              });

=cut

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
