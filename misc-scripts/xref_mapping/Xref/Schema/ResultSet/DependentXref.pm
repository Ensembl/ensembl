=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

package Xref::Schema::ResultSet::DependentXref;
use strict;
use warnings;

use parent 'DBIx::Class::ResultSet';


=head2 fetch_dependent_xref

Description: A canned query for fetching dependent xrefs by accession
             The result contains the result row from the query
Example    : my $row = $db->schema->resultset('DependentXref')
              ->fetch_dependent_xref("A0A075B6P5", "R-HSA-109582");
Returntype : Xref::Schema::Result::DependentXref

=cut

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
